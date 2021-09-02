#Internal function of speciation
fn.speciation <-
function (sites.list, n.indBYsite, n.new.sp) {
if (length(which(n.indBYsite >0))>=n.new.sp) { 
dat.abund <- table(as.character(sample(x = sites.list, 
prob = n.indBYsite, size = n.new.sp, replace = FALSE)))[sites.list]
names(dat.abund) <- sites.list
dat.abund[is.na(dat.abund)] <- 0
} else {
dat.abund <- table(as.character(sample(x = sites.list, 
prob = n.indBYsite, size = n.new.sp, replace = TRUE)))[sites.list]
names(dat.abund) <- sites.list
dat.abund[is.na(dat.abund)] <- 0
}
abundBYsite <- data.frame()[1:length(dat.abund),]
for (i in 1:length(dat.abund)) {
if (dat.abund[i] >0) {
new.spiBYsite <- dat.abund
new.spiBYsite[i] <- 1
new.spiBYsite[-i] <- 0
new.spiBYsite <- replicate(dat.abund[i], new.spiBYsite)
abundBYsite <- cbind(row.names = NULL, abundBYsite, new.spiBYsite)
}
}
return(abundBYsite)
}
#Internal function of habitat associations
fn.HA <-
function (pool.t0, E) {
d.temp.mean <- expand.grid(E, stringsAsFactors = FALSE, 
traits.mean = pool.t0$traits.mean)
d.temp.sd <- expand.grid(E, stringsAsFactors = FALSE, 
traits.sd = pool.t0$traits.sd)
d.temp <- cbind(d.temp.mean, d.temp.sd[, 2])
colnames(d.temp) <- c("E", "traits.mean", "traits.sd")
HA.siteBYind <- matrix(data = mapply(FUN = dnorm, 
x = d.temp$E, mean = d.temp$traits.mean, sd = d.temp$traits.sd), 
nrow = length(E), ncol = length(pool.t0$traits.mean), byrow = FALSE)
return(HA.siteBYind)
}
#Internal function of weighted lottery recruitment
fn.weighted.lottery <-
function (ind.list, recruitment.weights, J.remain) {
dat.abund <- table(as.character(sample(x = ind.list, 
prob = recruitment.weights, size = J.remain, replace = TRUE)))[ind.list]
names(dat.abund) <- ind.list
dat.abund[is.na(dat.abund)] <- 0
return(dat.abund)
}
#Function of the first time step
fn.first.timestep <-
function (Conditions.t0, E, nu, m, nb.mat, timestep) {
LC.t0.without.new.sp <- Conditions.t0$LC.t0
n.new.sp <- c()
for (i in 1:ncol(LC.t0.without.new.sp)) {
n.new.sp[i] <- sum(runif(sum(LC.t0.without.new.sp[,i])) < nu)
}
sites.list <- as.character(1:length(E))
n.sites <- length(E)
new.sp <- data.frame()[1:n.sites, ]
for (i in 1:ncol(LC.t0.without.new.sp)) {
if (n.new.sp[i] >0) {
new.spi <- fn.speciation(sites.list= sites.list, n.indBYsite = LC.t0.without.new.sp[,i], n.new.sp = n.new.sp[i])
names(new.spi) <- paste("t", timestep, "nsp.", names(LC.t0.without.new.sp[i]), "_", 1:n.new.sp[i], sep = "")
LC.t0.without.new.sp[,i] <- LC.t0.without.new.sp[,i]-rowSums(new.spi)
new.sp <- cbind(row.names = NULL, new.sp, new.spi)
}
}
LC.t0.RA <- as.matrix(LC.t0.without.new.sp/rowSums(LC.t0.without.new.sp))
I.RA <- nb.mat %*% LC.t0.RA
HA.siteBYind <- fn.HA(pool.t0= Conditions.t0$pool.t0, E)
HA.siteBYind <- HA.siteBYind/apply(HA.siteBYind,1,max)
LC <- HA.siteBYind * LC.t0.RA
LC <- LC/rowSums(LC)
WLR <- m * I.RA + (1 - m) * LC
WLR.list <- as.list(data.frame(t(WLR)))
J.remain.list <- as.list(rowSums(LC.t0.without.new.sp))
ind.list <- rownames(Conditions.t0$pool.t0)
LC.t1.without.new.sp <- data.frame(row.names = NULL, t(mapply(FUN = fn.weighted.lottery, recruitment.weights = WLR.list, J.remain=J.remain.list, 
MoreArgs = list(ind.list = ind.list))))
col_old <- colnames(LC.t1.without.new.sp) 
col_new <- gsub(pattern = "X",replacement = "", x  = col_old)
colnames(LC.t1.without.new.sp) <- col_new
LC.t1 <- cbind(row.names = NULL, LC.t1.without.new.sp, new.sp)
n.extinction <- length(which(colSums(LC.t1)==0, TRUE))
id.ancestors <- (colSums(LC.t1, na.rm=T) != 0)
LC.t1 <- LC.t1[,c(id.ancestors)]
names.new.sp <- names(new.sp)
pool.t1.new.sp <- data.frame(row.names = c(names.new.sp), species=c(names.new.sp))
id.new.sp <- gsub("^[^.]+.|_[^_]+$", "", names.new.sp)
select.traits.new.sp <- Conditions.t0$pool.t0[c(id.new.sp),c(2:4)]
pool.t1.new.sp <- cbind(pool.t1.new.sp, select.traits.new.sp)
pool.t1 <- rbind(Conditions.t0$pool.t0, pool.t1.new.sp)
pool.t1 <- pool.t1[c(id.ancestors),]
species.richness <- length(unique(pool.t1$species))
n.ancestors <- nrow(pool.t1)
n.groups <- length(unique(pool.t1$groups))
Conditions.t1 <- list(n.new.sp = sum(n.new.sp), n.extinction = n.extinction, species.richness = species.richness, n.ancestors = n.ancestors, n.groups = n.groups, LC.t0=LC.t1, pool.t0=pool.t1)
return(Conditions.t1)
}
#Calculate the species compositions of local communities and regional species pool forward in time
fn.forward.LC <- 
function(Conditions.t0, E, nu, m, nb.mat, n.timestep, keep = FALSE, stop = TRUE) {
Conditions.t.minus.1 <- Conditions.t0
n.new.sp <- c()
n.new.sp[1] <- Conditions.t0$n.new.sp
n.extinction <- c()
n.extinction[1] <- Conditions.t0$n.extinction 
species.richness <- c()
species.richness[1] <- Conditions.t0$species.richness 
n.ancestors <- c()
n.ancestors[1] <- Conditions.t0$n.ancestors
n.groups <- c()
n.groups[1] <- Conditions.t0$n.groups
if (keep) {
Conditions <- list()
Conditions[[1]] <- Conditions.t0
}
for (i in 2:n.timestep) {
Conditions.t <- fn.first.timestep(Conditions.t0 = Conditions.t.minus.1, E, nu, m, nb.mat, timestep=i)
n.new.sp[i] <- Conditions.t$n.new.sp
n.extinction[i] <- Conditions.t$n.extinction
species.richness[i] <- Conditions.t$species.richness
n.ancestors[i] <- Conditions.t$n.ancestors
n.groups[i] <- Conditions.t$n.groups
if (keep) {
Conditions[[i]] <- Conditions.t
}
print(paste("Timestep:", i))
if (stop) {
if (species.richness[i] >= n.ancestors[i]) { break }
}
Conditions.t.minus.1 <- Conditions.t
}
if (keep) 
return(list(Conditions = Conditions, n.new.sp=n.new.sp, n.extinction=n.extinction, species.richness=species.richness, n.ancestors=n.ancestors, n.groups=n.groups))
else return(list(Conditions = Conditions.t, n.new.sp=n.new.sp, n.extinction=n.extinction, species.richness=species.richness, n.ancestors=n.ancestors, n.groups=n.groups))    
}

