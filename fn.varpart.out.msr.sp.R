fn.varpart.out.msr.sp <- function(PCNM.pos, landscape, E.trans, pool.t0, LC.t0, listW) {
varpart.out <- data.frame(n.groups.out=NA, Df.env=NA, env.sign=NA, env.adjR2=NA, envspace.sign=NA, envspace.adjR2=NA, Df.pure.env=NA, pure.env.sign=NA, pure.env.adjR2=NA, Df.glob.PCNM=NA, glob.PCNM.sign=NA, glob.PCNM.adjR2=NA, Df.pure.glob.PCNM=NA, pure.glob.PCNM.sign=NA, pure.glob.PCNM.adjR2=NA, Df.all.glob=NA, all.glob.sign=NA, all.glob.adjR2=NA, Df.env.patch=NA, env.patch.sign=NA, env.patch.adjR2=NA, envspace.patch.sign=NA, envspace.patch.adjR2=NA, Df.pure.env.patch=NA, pure.env.patch.sign=NA, pure.env.patch.adjR2=NA, Df.PCNM.patch=NA, PCNM.patch.sign=NA, PCNM.patch.adjR2=NA, Df.pure.PCNM.patch=NA, pure.PCNM.patch.sign=NA, pure.PCNM.patch.adjR2=NA, Df.all.patch=NA, all.patch.sign=NA, all.patch.adjR2=NA)
varpart.out$n.groups.out <- length(unique(pool.t0$groups))
LC.h <- decostand(LC.t0, "hellinger")
LC.PCNM.rda <- rda(LC.h, PCNM.pos)
anova.LC.PCNM.rda <- anova.cca(LC.PCNM.rda)
LC.PCNM.R2a <- RsquareAdj(LC.PCNM.rda)$adj.r.squared
LC.PCNM.fwd <- forward.sel(LC.h, PCNM.pos, adjR2thresh=LC.PCNM.R2a, nperm=999)
PCNM.sign <- sort(LC.PCNM.fwd$order)
PCNM.red <- PCNM.pos[,c(PCNM.sign)]
PCNM.red <- as.data.frame(PCNM.red)
colnames(PCNM.red) <- LC.PCNM.fwd[order(LC.PCNM.fwd$order),]$variables
LC.PCNM.red.rda <- rda(LC.h, PCNM.red)  
anova.LC.PCNM.red.rda <- anova.cca(LC.PCNM.red.rda)
varpart.out$Df.glob.PCNM <- anova.LC.PCNM.red.rda$Df[1]
varpart.out$glob.PCNM.sign <- anova.LC.PCNM.red.rda$Pr[1]
varpart.out$glob.PCNM.adjR2<- RsquareAdj(LC.PCNM.red.rda)$adj.r.squared
env <- E.trans
LC.env.rda <- rda(LC.h ~., data=env)
anova.LC.env.rda <- anova.cca(LC.env.rda)
varpart.out$Df.env <- anova.LC.env.rda$Df[1]
R2.ab <- RsquareAdj(LC.env.rda)$r.squared
R2.ab.adj <- RsquareAdj(LC.env.rda)$adj.r.squared
R2.bc.adj <- RsquareAdj(LC.PCNM.red.rda)$adj.r.squared 
LC.pure.env.rda <- rda(LC.h, env, PCNM.red)
R2.a <- RsquareAdj(LC.pure.env.rda)$r.squared
R2.a.adj <- RsquareAdj(LC.pure.env.rda)$adj.r.squared
R2.b <- R2.ab- R2.a
R2.b.adj <- R2.ab.adj- R2.a.adj
alternative <- ifelse(R2.b < 0, "less", "greater")
E.ab <- c()
E.a <- c()
E.b <- c()
E.b.adj <- c()
MSR.ENV <- msr(env, listW, nrepet = 999, method = "singleton", simplify = FALSE)
for (k in 1:999) {
LC.varpart <- varpart(LC.h, MSR.ENV[[k]], PCNM.red)
E.ab <- c(E.ab, LC.varpart$part$fract$R.square[1])
E.a <- c(E.a, LC.varpart$part$fract$R.square[3]- LC.varpart$part$fract$R.square[2])
E.b <- c(E.b, LC.varpart$part$fract$R.square[3]- (LC.varpart$part$fract$R.square[3]- LC.varpart$part$fract$R.square[1])- (LC.varpart$part$fract$R.square[3]- LC.varpart$part$fract$R.square[2]))
E.b.adj <- c(E.b.adj, LC.varpart$part$indfract$Adj.R.square[2])
}
varpart.out$env.sign <- as.randtest(obs = R2.ab, sim = E.ab, alter = "greater")$pvalue
varpart.out$env.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))
global <- mem.select(env, listW, method = "global", MEM.autocor = "all", alpha = 0.05)$global.t
pval <- c(global$positive$pvalue, global$negative$pvalue)
if (! length(which(pval <= 0.05)) == 0) { 
varpart.out$envspace.sign <- as.randtest(obs = R2.b.adj, sim = E.b.adj, alter = alternative)$pvalue
}
varpart.out$envspace.adjR2 <- 1-(1- R2.b)/(1-mean(E.b))
varpart.out$pure.env.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))-(1-(1- R2.b)/(1-mean(E.b)))
anova.LC.pure.env.rda <- anova.cca(LC.pure.env.rda)
varpart.out$Df.pure.env <- anova.LC.pure.env.rda$Df[1]
varpart.out$pure.env.sign <- as.randtest(obs = R2.a, sim = E.a, alter = "greater")$pvalue
varpart.out$pure.glob.PCNM.adjR2 <- R2.bc.adj-(1-(1- R2.b)/(1-mean(E.b))) 
LC.pure.glob.PCNM.rda <- rda(LC.h, PCNM.red, env)
anova.LC.pure.glob.PCNM.rda <- anova.cca(LC.pure.glob.PCNM.rda)
varpart.out$Df.pure.glob.PCNM <- anova.LC.pure.glob.PCNM.rda$Df[1]
varpart.out$pure.glob.PCNM.sign <- anova.LC.pure.glob.PCNM.rda$Pr[1]
varpart.out$all.glob.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))+(R2.bc.adj-(1-(1- R2.b)/(1-mean(E.b))))
all.glob <- cbind(env, PCNM.red)
LC.all.glob.rda <- rda(LC.h, all.glob)
anova.LC.all.glob.rda <- anova.cca(LC.all.glob.rda)
varpart.out$Df.all.glob <- anova.LC.all.glob.rda$Df[1]
varpart.out$all.glob.sign <- anova.LC.all.glob.rda$Pr[1]
if (varpart.out$env.sign<=0.05 & length(unique(pool.t0$groups))>=2 |
anova.LC.env.rda$Pr[1] <=0.05 & anova.LC.pure.env.rda$Pr[1] <=0.05 & length(unique(pool.t0$groups))>=2) {
axes.test <- anova.cca(LC.env.rda, by="axis") 
choice <- length(which(!is.na(axes.test$Pr) & axes.test$Pr<=0.05))
LC.env.axes <- scores(LC.env.rda, choice=c(1:choice), display="lc", scaling=1)
LC.env.axes.KM.cascade <- try(cascadeKM(LC.env.axes, inf.gr=2, sup.gr= length(unique(pool.t0$groups)), iter=100, criterion="ssi"))
if (inherits(LC.env.axes.KM.cascade, "try-error"))
{
LC.env.axes.KM.cascade <- cascadeKM(LC.env.axes, inf.gr=2, sup.gr= nrow(distinct_all(as.data.frame(LC.env.axes))), iter=100, criterion="ssi")
}
if (ncol(LC.env.axes.KM.cascade$partition)==1) {
partition <- as.character(LC.env.axes.KM.cascade$partition)
} else {
partition <- as.character(LC.env.axes.KM.cascade$partition[,which(LC.env.axes.KM.cascade$results[2,] ==max(LC.env.axes.KM.cascade$results[2,] [sapply(apply(LC.env.axes.KM.cascade$partition, 2, table), min)!=1]))])
}
if (!length(partition)==nrow(LC.t0)) {
partition <- partition[c(1:nrow(LC.t0))]
}
env.patch <- as.data.frame(model.matrix(~partition)[,-1])
LC.env.patch.rda <- rda(LC.h, env.patch)
anova.LC.env.patch.rda <- anova.cca(LC.env.patch.rda)
varpart.out$Df.env.patch <- anova.LC.env.patch.rda$Df[1]
R2.ab <- RsquareAdj(LC.env.patch.rda)$r.squared
E.ab <- c()
MSR.ENV <- msr(env.patch, listW, nrepet = 999, method = "singleton", simplify = FALSE)
for (k in 1:999) {
E.ab <- c(E.ab, RsquareAdj(rda(LC.h, MSR.ENV[[k]]))$r.square)
}
env.patch.sign <- as.randtest(obs = R2.ab, sim = E.ab, alter = "greater")
if (env.patch.sign$pvalue>0.05) {
varpart.out$env.patch.sign <- env.patch.sign$pvalue
varpart.out$env.patch.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))
}
if (env.patch.sign$pvalue<=0.05) {
LC.groups <- data.frame(matrix(, nrow=length(landscape$E), ncol=0))
unique.groups <- unique(pool.t0$groups)
for(i in 1:length(unique.groups)) { 
if (nrow(pool.t0[which(pool.t0$groups==unique.groups[i]),]) >= 2) { 
LC.groups <- cbind(LC.groups, rowSums(LC.t0[,which(pool.t0$groups==unique.groups[i])]))
} else {
LC.groups <- cbind(LC.groups, LC.t0[,which(pool.t0$groups==unique.groups[i])])
}
}
colnames(LC.groups) <- paste("G_", unique.groups, sep = "")
env.patch.groups <- cbind(env.patch, LC.groups)
cor.env.patch.groups <- cor(env.patch.groups)
diag(cor.env.patch.groups) <- 0
patch.groups <- rownames(cor.env.patch.groups[apply(cor.env.patch.groups, 2, which.max), ])
patch.groups <- make.unique(patch.groups[1:ncol(env.patch)])
colnames(env.patch) <- patch.groups
if (ncol(env.patch)>=2) {
landscape.patch <- landscape[which(rowSums(env.patch)==0),]
} 
if (ncol(env.patch)==1) {
landscape.patch <- landscape[which(env.patch[,1]==0),]
}
dist.landscape.patch <- dist(landscape.patch[, c("x", "y")])
PCNM.patch <- PCNM(dist.landscape.patch)
select <- which(PCNM.patch$Moran_I$Positive==TRUE)
PCNM.patch.pos <- as.data.frame(PCNM.patch$vectors)[,select, drop=F]
colnames(PCNM.patch.pos) <- paste("others_PCNM", 1:ncol(PCNM.patch.pos), sep = "")
PCNM.patch.pos <- cbind(sites=landscape.patch[,1], PCNM.patch.pos)
for (i in 1:ncol(env.patch)) {
landscape.patchi <- landscape[which(env.patch[,i]==1),]
dist.landscape.patchi <- dist(landscape.patchi[, c("x", "y")])
PCNM.patchi <- PCNM(dist.landscape.patchi)
select <- which(PCNM.patchi$Moran_I$Positive==TRUE)
PCNM.patchi.pos <- as.data.frame(PCNM.patchi$vectors)[,select, drop=F]
colnames(PCNM.patchi.pos) <- paste(colnames(env.patch[i]), "_PCNM", 1:ncol(PCNM.patchi.pos), sep = "")
PCNM.patchi.pos <- cbind(sites=landscape.patchi[,1], PCNM.patchi.pos)
PCNM.patch.pos <- bind_rows(PCNM.patch.pos, PCNM.patchi.pos)
PCNM.patch.pos[is.na(PCNM.patch.pos)] <- 0
}
PCNM.patch.pos <- PCNM.patch.pos[order(PCNM.patch.pos$sites),]
PCNM.patch.pos <- as.data.frame(PCNM.patch.pos[,-1])
env.patch.pre <- predict(LC.env.patch.rda, type="response", newdata=env.patch, scaling=2)
residual <- LC.h - env.patch.pre
residual.PCNM.patch.rda <- rda(residual, PCNM.patch.pos)
anova.residual.PCNM.patch.rda <- anova.cca(residual.PCNM.patch.rda)
if (anova.residual.PCNM.patch.rda$Pr[1]<=0.05) {
residual.PCNM.patch.R2a <- RsquareAdj(residual.PCNM.patch.rda)$adj.r.squared
residual.PCNM.patch.fwd <- forward.sel(residual, PCNM.patch.pos, adjR2thresh=residual.PCNM.patch.R2a, nperm=999)
PCNM.patch.sign <- sort(residual.PCNM.patch.fwd$order)
PCNM.patch.red <- PCNM.patch.pos[,c(PCNM.patch.sign)]
PCNM.patch.red <- as.data.frame(PCNM.patch.red)
colnames(PCNM.patch.red) <- residual.PCNM.patch.fwd[order(residual.PCNM.patch.fwd$order),]$variables
LC.PCNM.patch.red.rda <- rda(LC.h, PCNM.patch.red)      
partition <- as.factor(partition)
anova.LC.PCNM.patch.red.rda <- permutest(LC.PCNM.patch.red.rda, permutations = how(nperm=999, within = Within(type="free"), plots = Plots(type="none", strata= partition)))
varpart.out$Df.PCNM.patch <- ncol(PCNM.patch.red)
varpart.out$PCNM.patch.sign <- as.randtest(obs = anova.LC.PCNM.patch.red.rda$F.0, sim = anova.LC.PCNM.patch.red.rda$F.perm, alter = "greater")$pvalue
num <- anova.LC.PCNM.patch.red.rda$F.perm*ncol(PCNM.patch.red)
den <- nrow(LC.h)-ncol(PCNM.patch.red)-1+num
R2.bc <- num/den
varpart.out$PCNM.patch.adjR2 <- 1-(1-RsquareAdj(LC.PCNM.patch.red.rda)$r.squared)/(1-mean(R2.bc))
R2.ab <- RsquareAdj(LC.env.patch.rda)$r.squared
R2.ab.adj <- RsquareAdj(LC.env.patch.rda)$adj.r.squared
R2.bc.adj <- 1-(1-RsquareAdj(LC.PCNM.patch.red.rda)$r.squared)/(1-mean(R2.bc))
LC.pure.env.patch.rda <- rda(LC.h, env.patch, PCNM.patch.red)
R2.a <- RsquareAdj(LC.pure.env.patch.rda)$r.squared
R2.a.adj <- RsquareAdj(LC.pure.env.patch.rda)$adj.r.squared
R2.b <- R2.ab- R2.a
R2.b.adj <- R2.ab.adj- R2.a.adj
alternative <- ifelse(R2.b < 0, "less", "greater")
E.ab <- c()
E.a <- c()
E.b <- c()
E.b.adj <- c()
MSR.ENV <- msr(env.patch, listW, nrepet = 999, method = "singleton", simplify = FALSE)
for (k in 1:999) {
LC.varpart <- varpart(LC.h, MSR.ENV[[k]], PCNM.patch.red)
E.ab <- c(E.ab, LC.varpart$part$fract$R.square[1])
E.a <- c(E.a, LC.varpart$part$fract$R.square[3]- LC.varpart$part$fract$R.square[2])
E.b <- c(E.b, LC.varpart$part$fract$R.square[3]- (LC.varpart$part$fract$R.square[3]- LC.varpart$part$fract$R.square[1])- (LC.varpart$part$fract$R.square[3]- LC.varpart$part$fract$R.square[2]))
E.b.adj <- c(E.b.adj, LC.varpart$part$indfract$Adj.R.square[2])
}
varpart.out$env.patch.sign <- as.randtest(obs = R2.ab, sim = E.ab, alter = "greater")$pvalue
varpart.out$env.patch.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))
global <- mem.select(env.patch, listW, method = "global", MEM.autocor = "all", alpha = 0.05)$global.t
pval <- c(global$positive$pvalue, global$negative$pvalue)
if (! length(which(pval <= 0.05)) == 0) { 
varpart.out$envspace.patch.sign <- as.randtest(obs = R2.b.adj, sim = E.b.adj, alter = alternative)$pvalue
}
varpart.out$envspace.patch.adjR2 <- 1-(1- R2.b)/(1-mean(E.b))
varpart.out$pure.env.patch.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))-(1-(1- R2.b)/(1-mean(E.b)))
anova.LC.pure.env.patch.rda <- anova.cca(LC.pure.env.patch.rda)
varpart.out$Df.pure.env.patch <- anova.LC.pure.env.patch.rda$Df[1]
varpart.out$pure.env.patch.sign <- as.randtest(obs = R2.a, sim = E.a, alter = "greater")$pvalue
varpart.out$pure.PCNM.patch.adjR2 <- R2.bc.adj-(1-(1- R2.b)/(1-mean(E.b))) 
LC.pure.patch.PCNM.rda <- rda(LC.h, PCNM.patch.red, env.patch)
anova.LC.pure.patch.PCNM.rda <- anova.cca(LC.pure.patch.PCNM.rda)
varpart.out$Df.pure.PCNM.patch <- anova.LC.pure.patch.PCNM.rda$Df[1]
varpart.out$pure.PCNM.patch.sign <- anova.LC.pure.patch.PCNM.rda$Pr[1]
varpart.out$all.patch.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))+(R2.bc.adj-(1-(1- R2.b)/(1-mean(E.b))))
all.patch <- cbind(env.patch, PCNM.patch.red)
LC.all.patch.rda <- try(rda(LC.h, all.patch))
if (inherits(LC.all.patch.rda, "try-error"))
{
all.patch <- cbind(PCNM.patch.red, env.patch)
LC.all.patch.rda <- rda(LC.h, all.patch)
}
anova.LC.all.patch.rda <- anova.cca(LC.all.patch.rda)
varpart.out$Df.all.patch <- anova.LC.all.patch.rda$Df[1]
varpart.out$all.patch.sign <- anova.LC.all.patch.rda$Pr[1]
}
if (anova.residual.PCNM.patch.rda$Pr[1]>0.05) {
LC.PCNM.patch.pos.rda <- rda(LC.h, PCNM.patch.pos)      
anova.LC.PCNM.patch.pos.rda <- anova.cca(LC.PCNM.patch.pos.rda)
varpart.out$Df.PCNM.patch <- anova.LC.PCNM.patch.pos.rda$Df[1]
varpart.out$PCNM.patch.sign <- anova.LC.PCNM.patch.pos.rda$Pr[1]
varpart.out$PCNM.patch.adjR2 <- RsquareAdj(LC.PCNM.patch.pos.rda)$adj.r.squared
varpart.out$Df.pure.PCNM.patch <- anova.LC.PCNM.patch.pos.rda$Df[1]
varpart.out$pure.PCNM.patch.sign <- anova.LC.PCNM.patch.pos.rda$Pr[1]
varpart.out$pure.PCNM.patch.adjR2 <- RsquareAdj(LC.PCNM.patch.pos.rda)$adj.r.squared
varpart.out$env.patch.sign <- env.patch.sign$pvalue
varpart.out$env.patch.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))
varpart.out$Df.pure.env.patch <- anova.LC.env.patch.rda$Df[1]
varpart.out$pure.env.patch.sign <- env.patch.sign$pvalue
varpart.out$pure.env.patch.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))
varpart.out$Df.all.patch <- anova.LC.env.patch.rda$Df[1]
varpart.out$all.patch.sign <- env.patch.sign$pvalue
varpart.out$all.patch.adjR2 <- 1-(1- R2.ab)/(1-mean(E.ab))
}
}
}
return(varpart.out)
}
