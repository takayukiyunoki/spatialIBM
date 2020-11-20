#Internal function to calculating the similarity of habitat associations between pairs of species
fn.dist.HA <- function(pool.t0) { 
pool.t0.HA.unique <- distinct(pool.t0[, 3:4])
d.temp.mean <- expand.grid(pool.t0.HA.unique$traits.mean, stringsAsFactors = FALSE, pool.t0.HA.unique$traits.mean)
d.temp.sd <- expand.grid(pool.t0.HA.unique$traits.sd, stringsAsFactors = FALSE, pool.t0.HA.unique$traits.sd)
d.temp <- cbind(d.temp.mean[, 1], d.temp.sd[, 1], d.temp.mean[, 2], d.temp.sd[, 2])
d.temp <- as.data.frame(d.temp)
for (i in 1:nrow(d.temp)) {
min.f1f2 <- function(x, mean1, mean2, sd1, sd2) {
f1 <- dnorm(x, mean=mean1, sd=sd1)
f2 <- dnorm(x, mean=mean2, sd=sd2)
pmin(f1, f2)
}
d.temp5 <- try(integrate(min.f1f2, -Inf, Inf, mean1=d.temp[i,1], mean2=d.temp[i,3], sd1=d.temp[i,2], sd2=d.temp[i,4])$value)
if (inherits(d.temp5, "try-error"))
{
d.temp5 <- integrate(min.f1f2, -Inf, Inf, mean1=d.temp[i,1], mean2=d.temp[i,3], sd1=d.temp[i,2], sd2=d.temp[i,4], rel.tol=1.e-10)$value
}
d.temp[i,5] <- d.temp5
}
pool.t0.HA <- distinct(pool.t0[, c(1,3:4)])
d.temp.mean2 <- expand.grid(pool.t0.HA$traits.mean, stringsAsFactors = FALSE, pool.t0.HA$traits.mean)
d.temp.sd2 <- expand.grid(pool.t0.HA$traits.sd, stringsAsFactors = FALSE, pool.t0.HA$traits.sd)
d.temp2 <- cbind(d.temp.mean2[, 1], d.temp.sd2[, 1], d.temp.mean2[, 2], d.temp.sd2[, 2])
d.temp2 <- as.data.frame(d.temp2)
d.temp3 <- left_join(d.temp2, d.temp, by = c("V1" = "V1", "V2" = "V2", "V3" = "V3", "V4" = "V4" ))
overlap <- matrix(data = d.temp3[,5], nrow = nrow(pool.t0.HA), ncol = nrow(pool.t0.HA), byrow = FALSE)
overlap <- overlap/max(overlap)
overlap[diag(overlap)] <- 1
return(overlap)
}
# Function to calculating the diversity statics of emergent biodiversity patterns  
fn.divpart.out <- function(pool.t0, LC.t0) {
overlap <- fn.dist.HA(pool.t0)
LC.t0 <- t(LC.t0)
MC <- MetaCommunity(Abundances = LC.t0, Weights = rep(1, ncol(LC.t0)))
dp.neutral <- DivProfile(seq(0, 2, 1), MC)
dp.functional <- DivProfile(seq(0, 2, 1), MC, Z = overlap)
divpart.out <- data.frame(q.order=c(0,1,2), alpha.neutral=dp.neutral$TotalAlphaDiversity, beta.neutral=dp.neutral$TotalBetaDiversity, gamma.neutral=dp.neutral$GammaDiversity, alpha.functional=dp.functional$TotalAlphaDiversity, beta.functional=dp.functional$TotalBetaDiversity, gamma.functional=dp.functional$GammaDiversity)
return(divpart.out)
}

