# For large system and 64 humps
n.sites = 1600
landscape <- data.frame(sites = 1:n.sites, expand.grid(x = 1:sqrt(n.sites), y = 1:sqrt(n.sites)))
landscape.hump16 <- read.csv("landscape.hump16.csv", header=TRUE, row.names=1)
hump16 <- matrix(landscape.hump16$E, byrow=FALSE, nrow=20, ncol=20)
hump32 <- matrix(rep(hump16,2), ncol=ncol(hump16), byrow=TRUE)
hump64 <- matrix(rep(hump32,2), ncol=ncol(hump32)*2, byrow=TRUE)
landscape$E <- as.vector(hump64)
write.csv(landscape, file = "landscape.hump64.csv")
landscape.hump64 <- landscape
par(mfrow = c(4,1))
sr.value(landscape.hump64[, c("x", "y")], landscape.hump64$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
landscape.center <- landscape.hump64[landscape.hump64$x >= 16 & landscape.hump64$x <= 25 & landscape.hump64$y >= 16 & landscape.hump64$y <= 25, ]
E.trans <- poly(landscape.center$E, 2)
colnames(E.trans) <- c("E1", "E2")
rownames(E.trans) <- rownames(landscape.center) 
E.trans <- as.data.frame(E.trans)
write.csv(E.trans, file = "E.trans.hump64.center.csv")
sr.value(landscape.center[, c("x", "y")], landscape.center$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E1, method = "greylevel", csize = 0.35, sub = "E1",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E2, method = "greylevel", csize = 0.35, sub = "E2",
csub=2)

