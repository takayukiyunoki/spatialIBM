# For large system and random
n.sites = 1600
landscape <- data.frame(sites = 1:n.sites, expand.grid(x = 1:sqrt(n.sites), y = 1:sqrt(n.sites)))
landscape.rand <- read.csv("landscape.rand.csv", header=TRUE, row.names=1)
rand <- matrix(landscape.rand$E, byrow=FALSE, nrow=20, ncol=20)
rand2 <- matrix(rep(rand,2), ncol=ncol(rand), byrow=TRUE)
rand4 <- matrix(rep(rand2,2), ncol=ncol(rand2)*2, byrow=TRUE)
landscape$E <- as.vector(rand4)
write.csv(landscape, file = "landscape.rand.large.csv")
landscape.rand.large <- landscape
par(mfrow = c(4,1))
sr.value(landscape.rand.large[, c("x", "y")], landscape.rand.large$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
landscape.center <- landscape.rand.large[landscape.rand.large$x >= 16 & landscape.rand.large$x <= 25 & landscape.rand.large$y >= 16 & landscape.rand.large$y <= 25, ]
E.trans <- poly(landscape.center$E, 2)
colnames(E.trans) <- c("E1", "E2")
rownames(E.trans) <- rownames(landscape.center) 
E.trans <- as.data.frame(E.trans)
write.csv(E.trans, file = "E.trans.rand.large.center.csv")
sr.value(landscape.center[, c("x", "y")], landscape.center$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E1, method = "greylevel", csize = 0.35, sub = "E1",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E2, method = "greylevel", csize = 0.35, sub = "E2",
csub=2)


