# Function to calculating the summary statics of variation partitioning for emergent biodiversity patterns
fn.varpart.out <- function(PCNM.pos, landscape, E.trans, pool.t0, LC.t0) {
varpart.out <- data.frame(n.groups.out=NA, Df.env=NA, env.sign=NA, env.adjR2=NA, env.cor.mean=NA, env.cor.sd=NA, env.correlog.mean=NA, env.correlog.sd=NA, env.mantel.cor=NA, env.mantel.sign=NA, Df.pure.env=NA, pure.env.sign=NA, pure.env.adjR2=NA, Df.glob.PCNM=NA, glob.PCNM.sign=NA, glob.PCNM.adjR2=NA, Df.pure.glob.PCNM=NA, pure.glob.PCNM.sign=NA, pure.glob.PCNM.adjR2=NA, neutral.glob.cor.mean=NA, neutral.glob.cor.sd=NA, neutral.glob.correlog.mean=NA, neutral.glob.correlog.sd=NA, neutral.glob.mantel.cor=NA, neutral.glob.mantel.sign=NA, Df.all.glob=NA, all.glob.sign=NA, all.glob.adjR2=NA, Df.env.patch=NA, env.patch.sign=NA, env.patch.adjR2=NA, Df.pure.env.patch=NA, pure.env.patch.sign=NA, pure.env.patch.adjR2=NA, Df.PCNM.patch=NA, PCNM.patch.sign=NA, PCNM.patch.adjR2=NA, Df.pure.PCNM.patch=NA, pure.PCNM.patch.sign=NA, pure.PCNM.patch.adjR2=NA, Df.all.patch=NA, all.patch.sign=NA, all.patch.adjR2=NA)
varpart.out$n.groups.out <- length(unique(pool.t0$groups))
LC.h <- decostand(LC.t0, "hellinger")
LC.PCNM.rda <- rda(LC.h, PCNM.pos)
anova.LC.PCNM.rda <- anova.cca(LC.PCNM.rda)
if (anova.LC.PCNM.rda$Pr[1]>0.05) {
varpart.out$Df.glob.PCNM <- anova.LC.PCNM.rda$Df[1]
varpart.out$glob.PCNM.sign <- anova.LC.PCNM.rda$Pr[1]
varpart.out$glob.PCNM.adjR2<- RsquareAdj(LC.PCNM.rda)$adj.r.squared
varpart.out$Df.pure.glob.PCNM <- anova.LC.PCNM.rda$Df[1]
varpart.out$pure.glob.PCNM.sign <- anova.LC.PCNM.rda$Pr[1]
varpart.out$pure.glob.PCNM.adjR2 <- RsquareAdj(LC.PCNM.rda)$adj.r.squared
} else {
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
} 
env <- E.trans
LC.env.rda <- rda(LC.h, env)
anova.LC.env.rda <- anova.cca(LC.env.rda)
if (anova.LC.env.rda$Pr[1]>0.05) {
varpart.out$Df.env <- anova.LC.env.rda$Df[1]
varpart.out$env.sign <- anova.LC.env.rda$Pr[1]
varpart.out$env.adjR2 <- RsquareAdj(LC.env.rda)$adj.r.squared
varpart.out$Df.pure.env <- anova.LC.env.rda$Df[1]
varpart.out$pure.env.sign <- anova.LC.env.rda$Pr[1]
varpart.out$pure.env.adjR2 <- RsquareAdj(LC.env.rda)$adj.r.squared
if (anova.LC.PCNM.rda$Pr[1]<=0.05) {
varpart.out$Df.pure.glob.PCNM <- anova.LC.PCNM.red.rda$Df[1]
varpart.out$pure.glob.PCNM.sign <- anova.LC.PCNM.red.rda$Pr[1]
varpart.out$pure.glob.PCNM.adjR2 <- RsquareAdj(LC.PCNM.red.rda)$adj.r.squared
varpart.out$Df.all.glob <- anova.LC.PCNM.red.rda$Df[1]
varpart.out$all.glob.sign <- anova.LC.PCNM.red.rda$Pr[1]
varpart.out$all.glob.adjR2 <- RsquareAdj(LC.PCNM.red.rda)$adj.r.squared
neutral.glob.pre <- predict(LC.PCNM.red.rda, type="response", newdata= PCNM.red, scaling=2)
cor.neutral.glob.pre <- cor(neutral.glob.pre)
correlogram <- correlog(landscape[,c("x", "y")], neutral.glob.pre[,1], method="Moran", nbclass = NULL)[,2]
for (i in 2:ncol(neutral.glob.pre)) 
correlogram <- cbind(correlogram, correlog(landscape[,c("x", "y")], neutral.glob.pre[,i], method="Moran", nbclass = NULL)[,2])
correlogram <- t(correlogram)
correlogram.MD <- dist(correlogram,method = "manhattan")
correlogram.MD <- as.matrix(correlogram.MD)
mantel.Pr <- mantel.test(cor.neutral.glob.pre, correlogram.MD,nperm = 999, alternative = "less")
varpart.out$neutral.glob.mantel.sign <- mantel.Pr$p 
correlation.lowtri <- cor.neutral.glob.pre[lower.tri(cor.neutral.glob.pre)]
varpart.out$neutral.glob.cor.mean <- mean(correlation.lowtri)
varpart.out$neutral.glob.cor.sd <-sd(correlation.lowtri)
correlogram.lowtri <- correlogram.MD[lower.tri(correlogram.MD)]
varpart.out$neutral.glob.correlog.mean <- mean(correlogram.lowtri)
varpart.out$neutral.glob.correlog.sd <-sd(correlogram.lowtri)
correlogram.MD <- as.dist(correlogram.MD)
cor.neutral.glob.pre <- as.dist(cor.neutral.glob.pre)
mantel.cor <- mantel(cor.neutral.glob.pre, correlogram.MD,method="pearson", permutations=999)
varpart.out$neutral.glob.mantel.cor  <- mantel.cor$statistic
} 
}
if (anova.LC.env.rda$Pr[1]<=0.05) {
LC.env.R2a <- RsquareAdj(LC.env.rda)$adj.r.squared
LC.env.fwd <- try(forward.sel(LC.h, env, adjR2thresh=LC.env.R2a, nperm=999))
if (!inherits(LC.env.fwd, "try-error"))
{
env.sign <- sort(LC.env.fwd$order)
env.red <- env[,c(env.sign)]
env.red <- as.data.frame(env.red)
colnames(env.red) <- LC.env.fwd[order(LC.env.fwd$order),]$variables
} 
if (inherits(LC.env.fwd, "try-error"))
{
env.red <- env
}
LC.env.red.rda <- rda(LC.h ~., data=env.red)
anova.LC.env.red.rda <- anova.cca(LC.env.red.rda)
varpart.out$Df.env <- anova.LC.env.red.rda$Df[1]
varpart.out$env.sign <- anova.LC.env.red.rda$Pr[1]
varpart.out$env.adjR2 <- RsquareAdj(LC.env.red.rda)$adj.r.squared
if (anova.LC.PCNM.rda$Pr[1]<=0.05) {
LC.varpart <- varpart(LC.h, env.red, PCNM.red)
varpart.out$pure.env.adjR2 <- LC.varpart$part$indfract$Adj.R.squared[1]
LC.pure.env.rda <- rda(LC.h, env.red, PCNM.red)
anova.LC.pure.env.rda <- anova.cca(LC.pure.env.rda)
varpart.out$Df.pure.env <- anova.LC.pure.env.rda$Df[1]
varpart.out$pure.env.sign <- anova.LC.pure.env.rda$Pr[1]
varpart.out$pure.glob.PCNM.adjR2 <- LC.varpart$part$indfract$Adj.R.squared[3]
LC.pure.glob.PCNM.rda <- rda(LC.h, PCNM.red, env.red)
anova.LC.pure.glob.PCNM.rda <- anova.cca(LC.pure.glob.PCNM.rda)
varpart.out$Df.pure.glob.PCNM <- anova.LC.pure.glob.PCNM.rda$Df[1]
varpart.out$pure.glob.PCNM.sign <- anova.LC.pure.glob.PCNM.rda$Pr[1]
varpart.out$all.glob.adjR2 <- LC.varpart$part$fract$Adj.R.squared[3]
all.glob <- cbind(env.red, PCNM.red)
LC.all.glob.rda <- rda(LC.h, all.glob)
anova.LC.all.glob.rda <- anova.cca(LC.all.glob.rda)
varpart.out$Df.all.glob <- anova.LC.all.glob.rda$Df[1]
varpart.out$all.glob.sign <- anova.LC.all.glob.rda$Pr[1]
if (anova.LC.pure.env.rda$Pr[1]<=0.05 & anova.LC.pure.glob.PCNM.rda$Pr[1] <=0.05) {
all.glob.pre <- predict(LC.all.glob.rda, type="response", newdata=all.glob, scaling=2)
env.pre <- predict(LC.env.red.rda, type="response", newdata=env.red, scaling=2)
neutral.glob.pre <- all.glob.pre - env.pre
env.pre <- as.data.frame(env.pre)
env.pre <- env.pre[, sapply(env.pre, var) !=0]
cor.env.pre <- cor(env.pre)
correlogram.env <- correlog(landscape[,c("x", "y")], env.pre[,1], method="Moran", nbclass = NULL)[,2]
for (i in 2:ncol(env.pre)) 
correlogram.env <- cbind(correlogram.env, correlog(landscape[,c("x", "y")], env.pre[,i], method="Moran", nbclass = NULL)[,2])
correlogram.env <- t(correlogram.env)
correlogram.MD.env <- dist(correlogram.env,method = "manhattan")
correlogram.MD.env <- as.matrix(correlogram.MD.env)
mantel.Pr.env <- mantel.test(cor.env.pre, correlogram.MD.env,nperm = 999, alternative = "less")
varpart.out$env.mantel.sign <- mantel.Pr.env$p
correlation.env.lowtri <- cor.env.pre[lower.tri(cor.env.pre)]
varpart.out$env.cor.mean <- mean(correlation.env.lowtri)
varpart.out$env.cor.sd <-sd(correlation.env.lowtri)
correlogram.env.lowtri <- correlogram.MD.env[lower.tri(correlogram.MD.env)]
varpart.out$env.correlog.mean <- mean(correlogram.env.lowtri)
varpart.out$env.correlog.sd <-sd(correlogram.env.lowtri)
correlogram.MD.env <- as.dist(correlogram.MD.env)
cor.env.pre <- as.dist(cor.env.pre)
mantel.cor.env <- mantel(cor.env.pre, correlogram.MD.env,method="pearson", permutations=999)
varpart.out$env.mantel.cor <- mantel.cor.env$statistic
cor.neutral.glob.pre <- cor(neutral.glob.pre)
correlogram <- correlog(landscape[,c("x", "y")], neutral.glob.pre[,1], method="Moran", nbclass = NULL)[,2]
for (i in 2:ncol(neutral.glob.pre)) 
correlogram <- cbind(correlogram, correlog(landscape[,c("x", "y")], neutral.glob.pre[,i], method="Moran", nbclass = NULL)[,2])
correlogram <- t(correlogram)
correlogram.MD <- dist(correlogram,method = "manhattan")
correlogram.MD <- as.matrix(correlogram.MD)
mantel.Pr <- mantel.test(cor.neutral.glob.pre, correlogram.MD,nperm = 999, alternative = "less")
varpart.out$neutral.glob.mantel.sign <- mantel.Pr$p 
correlation.lowtri <- cor.neutral.glob.pre[lower.tri(cor.neutral.glob.pre)]
varpart.out$neutral.glob.cor.mean <- mean(correlation.lowtri)
varpart.out$neutral.glob.cor.sd <-sd(correlation.lowtri)
correlogram.lowtri <- correlogram.MD[lower.tri(correlogram.MD)]
varpart.out$neutral.glob.correlog.mean <- mean(correlogram.lowtri)
varpart.out$neutral.glob.correlog.sd <-sd(correlogram.lowtri)
correlogram.MD <- as.dist(correlogram.MD)
cor.neutral.glob.pre <- as.dist(cor.neutral.glob.pre)
mantel.cor <- mantel(cor.neutral.glob.pre, correlogram.MD,method="pearson", permutations=999)
varpart.out$neutral.glob.mantel.cor  <- mantel.cor$statistic
}
if (anova.LC.pure.env.rda$Pr[1]>0.05 & anova.LC.pure.glob.PCNM.rda$Pr[1] <=0.05) {
neutral.glob.pre <- predict(LC.PCNM.red.rda, type="response", newdata= PCNM.red, scaling=2)
cor.neutral.glob.pre <- cor(neutral.glob.pre)
correlogram <- correlog(landscape[,c("x", "y")], neutral.glob.pre[,1], method="Moran", nbclass = NULL)[,2]
for (i in 2:ncol(neutral.glob.pre)) 
correlogram <- cbind(correlogram, correlog(landscape[,c("x", "y")], neutral.glob.pre[,i], method="Moran", nbclass = NULL)[,2])
correlogram <- t(correlogram)
correlogram.MD <- dist(correlogram,method = "manhattan")
correlogram.MD <- as.matrix(correlogram.MD)
mantel.Pr <- mantel.test(cor.neutral.glob.pre, correlogram.MD,nperm = 999, alternative = "less")
varpart.out$neutral.glob.mantel.sign <- mantel.Pr$p 
correlation.lowtri <- cor.neutral.glob.pre[lower.tri(cor.neutral.glob.pre)]
varpart.out$neutral.glob.cor.mean <- mean(correlation.lowtri)
varpart.out$neutral.glob.cor.sd <-sd(correlation.lowtri)
correlogram.lowtri <- correlogram.MD[lower.tri(correlogram.MD)]
varpart.out$neutral.glob.correlog.mean <- mean(correlogram.lowtri)
varpart.out$neutral.glob.correlog.sd <-sd(correlogram.lowtri)
correlogram.MD <- as.dist(correlogram.MD)
cor.neutral.glob.pre <- as.dist(cor.neutral.glob.pre)
mantel.cor <- mantel(cor.neutral.glob.pre, correlogram.MD,method="pearson", permutations=999)
varpart.out$neutral.glob.mantel.cor  <- mantel.cor$statistic
} 
} 
if (anova.LC.PCNM.rda$Pr[1]>0.05) {
varpart.out$Df.pure.env <- anova.LC.env.red.rda$Df[1]
varpart.out$pure.env.sign <- anova.LC.env.red.rda$Pr[1]
varpart.out$pure.env.adjR2 <- RsquareAdj(LC.env.red.rda)$adj.r.squared
varpart.out$Df.all.glob <- anova.LC.env.red.rda$Df[1]
varpart.out$all.glob.sign <- anova.LC.env.red.rda$Pr[1]
varpart.out$all.glob.adjR2 <- RsquareAdj(LC.env.red.rda)$adj.r.squared
env.pre <- predict(LC.env.red.rda, type="response", newdata=env.red, scaling=2)
env.pre <- as.data.frame(env.pre)
env.pre <- env.pre[, sapply(env.pre, var) !=0]
cor.env.pre <- cor(env.pre)
correlogram.env <- correlog(landscape[,c("x", "y")], env.pre[,1], method="Moran", nbclass = NULL)[,2]
for (i in 2:ncol(env.pre)) 
correlogram.env <- cbind(correlogram.env, correlog(landscape[,c("x", "y")], env.pre[,i], method="Moran", nbclass = NULL)[,2])
correlogram.env <- t(correlogram.env)
correlogram.MD.env <- dist(correlogram.env,method = "manhattan")
correlogram.MD.env <- as.matrix(correlogram.MD.env)
mantel.Pr.env <- mantel.test(cor.env.pre, correlogram.MD.env,nperm = 999, alternative = "less")
varpart.out$env.mantel.sign <- mantel.Pr.env$p
correlation.env.lowtri <- cor.env.pre[lower.tri(cor.env.pre)]
varpart.out$env.cor.mean <- mean(correlation.env.lowtri)
varpart.out$env.cor.sd <-sd(correlation.env.lowtri)
correlogram.env.lowtri <- correlogram.MD.env[lower.tri(correlogram.MD.env)]
varpart.out$env.correlog.mean <- mean(correlogram.env.lowtri)
varpart.out$env.correlog.sd <-sd(correlogram.env.lowtri)
correlogram.MD.env <- as.dist(correlogram.MD.env)
cor.env.pre <- as.dist(cor.env.pre)
mantel.cor.env <- mantel(cor.env.pre, correlogram.MD.env,method="pearson", permutations=999)
varpart.out$env.mantel.cor <- mantel.cor.env$statistic
}
}
if (varpart.out$pure.env.sign<=0.05 & length(unique(pool.t0$groups))>=2) {
axes.test <- anova.cca(LC.env.red.rda, by="axis") 
choice <- length(which(!is.na(axes.test$Pr) & axes.test$Pr<=0.05))
LC.env.axes <- scores(LC.env.red.rda, choice=c(1:choice), display="lc", scaling=1)
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
if (anova.LC.env.patch.rda$Pr[1]<=0.05 & ncol(env.patch)>=2) {
LC.env.patch.R2a <- RsquareAdj(LC.env.patch.rda)$adj.r.squared
LC.env.patch.fwd <- try(forward.sel(LC.h, env.patch, adjR2thresh=LC.env.patch.R2a, nperm=999))
if (!inherits(LC.env.patch.fwd, "try-error"))
{
env.patch.sign <- sort(LC.env.patch.fwd$order)
env.patch.red <- env.patch[,c(env.patch.sign)]
env.patch.red <- as.data.frame(env.patch.red)
colnames(env.patch.red) <- LC.env.patch.fwd[order(LC.env.patch.fwd$order),]$variables
} 
if (inherits(LC.env.patch.fwd, "try-error"))
{
env.patch.red <- env.patch
}
LC.env.patch.red.rda <- rda(LC.h, env.patch.red)
anova.LC.env.patch.red.rda <- anova.cca(LC.env.patch.red.rda)
varpart.out$Df.env.patch <- anova.LC.env.patch.red.rda$Df[1]
varpart.out$env.patch.sign <- anova.LC.env.patch.red.rda$Pr[1]
varpart.out$env.patch.adjR2 <- RsquareAdj(LC.env.patch.red.rda)$adj.r.squared
}
if (anova.LC.env.patch.rda$Pr[1]<=0.05 & ncol(env.patch)==1) {
env.patch.red <- env.patch
LC.env.patch.red.rda <- rda(LC.h, env.patch.red)
anova.LC.env.patch.red.rda <- anova.cca(LC.env.patch.red.rda)
varpart.out$Df.env.patch <- anova.LC.env.patch.red.rda$Df[1]
varpart.out$env.patch.sign <- anova.LC.env.patch.red.rda$Pr[1]
varpart.out$env.patch.adjR2 <- RsquareAdj(LC.env.patch.red.rda)$adj.r.squared
}
if (!anova.LC.env.patch.rda$Pr[1]<=0.05) {
varpart.out$Df.env.patch <- anova.LC.env.patch.rda$Df[1]
varpart.out$env.patch.sign <- anova.LC.env.patch.rda$Pr[1]
varpart.out$env.patch.adjR2 <- RsquareAdj(LC.env.patch.rda)$adj.r.squared
}
if (anova.LC.env.patch.rda$Pr[1]<=0.05) {
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
env.patch.groups <- cbind(env.patch.red, LC.groups)
cor.env.patch.groups <- cor(env.patch.groups)
diag(cor.env.patch.groups) <- 0
patch.groups <- rownames(cor.env.patch.groups[apply(cor.env.patch.groups, 2, which.max), ])
patch.groups <- make.unique(patch.groups[1:ncol(env.patch.red)])
colnames(env.patch.red) <- patch.groups
if (ncol(env.patch.red)>=2) {
landscape.patch <- landscape[which(rowSums(env.patch.red)==0),]
} 
if (ncol(env.patch.red)==1) {
landscape.patch <- landscape[which(env.patch.red[,1]==0),]
}
dist.landscape.patch <- dist(landscape.patch[, c("x", "y")])
PCNM.patch <- PCNM(dist.landscape.patch)
select <- which(PCNM.patch$Moran_I$Positive==TRUE)
PCNM.patch.pos <- as.data.frame(PCNM.patch$vectors)[,select, drop=F]
colnames(PCNM.patch.pos) <- paste("others_PCNM", 1:ncol(PCNM.patch.pos), sep = "")
PCNM.patch.pos <- cbind(sites=landscape.patch[,1], PCNM.patch.pos)
for (i in 1:ncol(env.patch.red)) {
landscape.patchi <- landscape[which(env.patch.red[,i]==1),]
dist.landscape.patchi <- dist(landscape.patchi[, c("x", "y")])
PCNM.patchi <- PCNM(dist.landscape.patchi)
select <- which(PCNM.patchi$Moran_I$Positive==TRUE)
PCNM.patchi.pos <- as.data.frame(PCNM.patchi$vectors)[,select, drop=F]
colnames(PCNM.patchi.pos) <- paste(colnames(env.patch.red[i]), "_PCNM", 1:ncol(PCNM.patchi.pos), sep = "")
PCNM.patchi.pos <- cbind(sites=landscape.patchi[,1], PCNM.patchi.pos)
PCNM.patch.pos <- bind_rows(PCNM.patch.pos, PCNM.patchi.pos)
PCNM.patch.pos[is.na(PCNM.patch.pos)] <- 0
}
PCNM.patch.pos <- PCNM.patch.pos[order(PCNM.patch.pos$sites),]
PCNM.patch.pos <- as.data.frame(PCNM.patch.pos[,-1])
env.patch.pre <- predict(LC.env.patch.red.rda, type="response", newdata=env.patch.red, scaling=2)
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
anova.LC.PCNM.patch.red.rda <- anova.cca(LC.PCNM.patch.red.rda)
varpart.out$Df.PCNM.patch <- anova.LC.PCNM.patch.red.rda$Df[1]
varpart.out$PCNM.patch.sign <- anova.LC.PCNM.patch.red.rda$Pr[1]
varpart.out$PCNM.patch.adjR2 <- RsquareAdj(LC.PCNM.patch.red.rda)$adj.r.squared
LC.patch.varpart <- varpart(LC.h, env.patch.red, PCNM.patch.red)
varpart.out$pure.env.patch.adjR2 <- LC.patch.varpart$part$indfract$Adj.R.squared[1]
LC.pure.env.patch.rda <- rda(LC.h, env.patch.red, PCNM.patch.red)
anova.LC.pure.env.patch.rda <- anova.cca(LC.pure.env.patch.rda)
varpart.out$Df.pure.env.patch <- anova.LC.pure.env.patch.rda$Df[1]
varpart.out$pure.env.patch.sign <- anova.LC.pure.env.patch.rda$Pr[1]
varpart.out$pure.PCNM.patch.adjR2 <- LC.patch.varpart$part$indfract$Adj.R.squared[3]
LC.pure.PCNM.patch.rda <- rda(LC.h, PCNM.patch.red, env.patch.red)
anova.LC.pure.PCNM.patch.rda <- anova.cca(LC.pure.PCNM.patch.rda)
varpart.out$Df.pure.PCNM.patch <- anova.LC.pure.PCNM.patch.rda$Df[1]
varpart.out$pure.PCNM.patch.sign <- anova.LC.pure.PCNM.patch.rda$Pr[1]
varpart.out$all.patch.adjR2 <- LC.patch.varpart$part$fract$Adj.R.squared[3]
all.patch <- cbind(env.patch.red, PCNM.patch.red)
LC.all.patch.rda <- try(rda(LC.h, all.patch))
if (inherits(LC.all.patch.rda, "try-error"))
{
all.patch <- cbind(PCNM.patch.red, env.patch.red)
LC.all.patch.rda <- rda(LC.h, all.patch)
}
anova.LC.all.patch.rda <- anova.cca(LC.all.patch.rda)
varpart.out$Df.all.patch <- anova.LC.all.patch.rda$Df[1]
varpart.out$all.patch.sign <- anova.LC.all.patch.rda$Pr[1]
all.patch.pre <- predict(LC.all.patch.rda, type="response", newdata=all.patch, scaling=2)
neutral.patch.pre <- all.patch.pre - env.patch.pre
env.patch.groups <- c("others", colnames(env.patch.red))
others <- rep(1, nrow(landscape))-rowSums(env.patch.red)
others.env.patch.red <- cbind(others, env.patch.red)
PCNMbygroups <- list()
varpart.out.groups <- data.frame(matrix(, nrow=1, ncol=0))
for (i in 1:length(env.patch.groups)) {
varpart.out.groupsi <- data.frame(matrix(, nrow=1, ncol=0))
PCNMbygroups[[i]] <- as.data.frame(PCNM.patch.red[,grepl(paste("^",env.patch.groups[i],"$", sep=""), gsub("_PCNM.*","",names(PCNM.patch.red)))])
if (ncol(PCNMbygroups[[i]]) >= 1) { 
neutral.patch.prei <- neutral.patch.pre[which(others.env.patch.red[,i]==1),]
landscapei <- landscape[which(others.env.patch.red[,i]==1),]
cor.neutral.patch.prei <- cor(neutral.patch.prei)
correlogrami <- try(correlog(landscapei[,c("x", "y")], neutral.patch.prei[,1], method="Moran", nbclass = 3)[,2])
if (!inherits(correlogrami, "try-error"))
{
for (j in 2:ncol(neutral.patch.prei)) 
correlogrami <- try(cbind(correlogrami, correlog(landscapei[,c("x", "y")], neutral.patch.prei[,j], method="Moran", nbclass = 3)[,2])) 
varpart.out.groupsi$groups <- env.patch.groups[i]
varpart.out.groupsi$PCNM.correlog.n.c <- 3
} 
if (inherits(correlogrami, "try-error"))
{
correlogrami <- correlog(landscapei[,c("x", "y")], neutral.patch.prei[,1], method="Moran", nbclass = 2)[,2]
for (j in 2:ncol(neutral.patch.prei)) 
correlogrami <- cbind(correlogrami, correlog(landscapei[,c("x", "y")], neutral.patch.prei[,j], method="Moran", nbclass = 2)[,2])
varpart.out.groupsi$groups <- env.patch.groups[i]
varpart.out.groupsi$PCNM.correlog.n.c <- 2
} 
correlogrami <- t(correlogrami)
correlogram.MDi <- dist(correlogrami,method = "manhattan")
correlation.lowtrii <- cor.neutral.patch.prei[lower.tri(cor.neutral.patch.prei)]
correlogram.lowtrii <- as.matrix(correlogram.MDi)[lower.tri(as.matrix(correlogram.MDi))]
varpart.out.groupsi$PCNM.cor.mean <- mean(correlation.lowtrii)
varpart.out.groupsi$PCNM.cor.sd <- sd(correlation.lowtrii)
varpart.out.groupsi$PCNM.correlog.mean <- mean(correlogram.lowtrii)
varpart.out.groupsi$PCNM.correlog.sd <- sd(correlogram.lowtrii)  
cor.neutral.patch.prei <- as.dist(cor.neutral.patch.prei)
mantel.cori <- mantel(cor.neutral.patch.prei, correlogram.MDi,method="pearson", permutations=999)
varpart.out.groupsi$PCNM.mantel.cor  <- mantel.cori$statistic
correlogram.MDi <- as.matrix(correlogram.MDi)
cor.neutral.patch.prei <- as.matrix(cor.neutral.patch.prei)
mantel.Pri <- mantel.test(cor.neutral.patch.prei, correlogram.MDi,nperm = 999, alternative = "less")
varpart.out.groupsi$PCNM.mantel.sign <- mantel.Pri$p
varpart.out.groups <- cbind(varpart.out.groups, varpart.out.groupsi)
}
}
} 
if (anova.residual.PCNM.patch.rda$Pr[1]>0.05) {
LC.PCNM.patch.pos.rda <- rda(LC.h, PCNM.patch.pos)	
anova.LC.PCNM.patch.pos.rda <- anova.cca(LC.PCNM.patch.pos.rda)
varpart.out$Df.PCNM.patch <- anova.LC.PCNM.patch.pos.rda$Df[1]
varpart.out$PCNM.patch.sign <- anova.LC.PCNM.patch.pos.rda$Pr[1]
varpart.out$PCNM.patch.adjR2 <- RsquareAdj(LC.PCNM.patch.pos.rda)$adj.r.squared
LC.patch.varpart <- varpart(LC.h, env.patch.red, PCNM.patch.pos)
varpart.out$pure.env.patch.adjR2 <- LC.patch.varpart$part$indfract$Adj.R.squared[1]
LC.pure.env.patch.rda <- rda(LC.h, env.patch.red, PCNM.patch.pos)
anova.LC.pure.env.patch.rda <- anova.cca(LC.pure.env.patch.rda)
varpart.out$Df.pure.env.patch <- anova.LC.pure.env.patch.rda$Df[1]
varpart.out$pure.env.patch.sign <- anova.LC.pure.env.patch.rda$Pr[1]
varpart.out$pure.PCNM.patch.adjR2 <- LC.patch.varpart$part$indfract$Adj.R.squared[3]
LC.pure.PCNM.patch.rda <- rda(LC.h, PCNM.patch.pos, env.patch.red)
anova.LC.pure.PCNM.patch.rda <- anova.cca(LC.pure.PCNM.patch.rda)
varpart.out$Df.pure.PCNM.patch <- anova.LC.pure.PCNM.patch.rda$Df[1]
varpart.out$pure.PCNM.patch.sign <- anova.LC.pure.PCNM.patch.rda$Pr[1]
varpart.out$all.patch.adjR2 <- LC.patch.varpart$part$fract$Adj.R.squared[3]
all.patch <- cbind(env.patch.red, PCNM.patch.pos)
LC.all.patch.rda <- try(rda(LC.h, all.patch))
if (inherits(LC.all.patch.rda, "try-error"))
{
all.patch <- cbind(PCNM.patch.pos, env.patch.red)
LC.all.patch.rda <- rda(LC.h, all.patch)
}
anova.LC.all.patch.rda <- anova.cca(LC.all.patch.rda)
varpart.out$Df.all.patch <- anova.LC.all.patch.rda$Df[1]
varpart.out$all.patch.sign <- anova.LC.all.patch.rda$Pr[1]
}
}
}
if (!varpart.out$pure.env.sign<=0.05 | !length(unique(pool.t0$groups))>=2) {
return(varpart.out)
} else if (!anova.LC.env.patch.rda$Pr[1]<=0.05) {
return(varpart.out)
} else if (!anova.residual.PCNM.patch.rda$Pr[1]<=0.05) {
return(varpart.out)
} else {
varpart.out <- cbind(varpart.out, varpart.out.groups)
return(varpart.out)
} 
}
