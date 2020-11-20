# For large system and two waves
n.sites = 1600
landscape <- data.frame(sites = 1:n.sites, expand.grid(x = 1:sqrt(n.sites), y = 1:sqrt(n.sites)))
landscape.wave <- read.csv("landscape.wave.csv", header=TRUE, row.names=1)
wave <- matrix(landscape.wave$E, byrow=FALSE, nrow=20, ncol=20)
wave2 <- matrix(rep(wave,2), ncol=ncol(wave), byrow=TRUE)
wave2 <- matrix(rep(wave2,2), ncol=ncol(wave2)*2, byrow=TRUE)
landscape$E <- as.vector(wave2)
write.csv(landscape, file = "landscape.wave2.csv")
landscape.wave2 <- landscape
par(mfrow = c(4,1))
sr.value(landscape.wave2[, c("x", "y")], landscape.wave2$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
landscape.center <- landscape.wave2[landscape.wave2$x >= 16 & landscape.wave2$x <= 25 & landscape.wave2$y >= 16 & landscape.wave2$y <= 25, ]
E.trans <- poly(landscape.center$E, 2)
colnames(E.trans) <- c("E1", "E2")
rownames(E.trans) <- rownames(landscape.center) 
E.trans <- as.data.frame(E.trans)
write.csv(E.trans, file = "E.trans.wave2.center.csv")
sr.value(landscape.center[, c("x", "y")], landscape.center$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E1, method = "greylevel", csize = 0.35, sub = "E1",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E2, method = "greylevel", csize = 0.35, sub = "E2",
csub=2)


