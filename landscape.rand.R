# For small system and random
n.sites = 400
landscape <- data.frame(sites = 1:n.sites, expand.grid(x = 1:sqrt(n.sites), y = 1:sqrt(n.sites)))
set.seed(1) # For replicability
landscape$E <- sample(rep(1:10, each = n.sites/10))
landscape$E <- landscape$E/10
write.csv(landscape, file = "landscape.rand.csv")
landscape.rand <- landscape
par(mfrow = c(4,1))
sr.value(landscape.rand[, c("x", "y")], landscape.rand$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
landscape.center <- landscape.rand[landscape.rand$x >= 6 & landscape.rand$x <= 15 & landscape.rand$y >= 6 & landscape.rand$y <= 15, ]
E.trans <- poly(landscape.center$E, 2)
colnames(E.trans) <- c("E1", "E2")
rownames(E.trans) <- rownames(landscape.center) 
E.trans <- as.data.frame(E.trans)
write.csv(E.trans, file = "E.trans.rand.center.csv")
sr.value(landscape.center[, c("x", "y")], landscape.center$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E1, method = "greylevel", csize = 0.35, sub = "E1",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E2, method = "greylevel", csize = 0.35, sub = "E2",
csub=2)
 
