# For small system and one wave
n.sites = 400
landscape <- data.frame(sites = 1:n.sites, expand.grid(x = 1:sqrt(n.sites), y = 1:sqrt(n.sites)))
landscape$E <- rep(c(6,7,8,9,10,10,9,8,7,6,5,4,3,2,1,1,2,3,4,5), each = 20)
landscape$E <- landscape$E/10
write.csv(landscape, file = "landscape.wave.csv")
landscape.wave <- landscape
par(mfrow = c(4,1))
sr.value(landscape.wave[, c("x", "y")], landscape.wave$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
landscape.center <- landscape.wave[landscape.wave$x >= 6 & landscape.wave$x <= 15 & landscape.wave$y >= 6 & landscape.wave$y <= 15, ]
E.trans <- poly(landscape.center$E, 2)
colnames(E.trans) <- c("E1", "E2")
rownames(E.trans) <- rownames(landscape.center) 
E.trans <- as.data.frame(E.trans)
write.csv(E.trans, file = "E.trans.wave.center.csv")
sr.value(landscape.center[, c("x", "y")], landscape.center$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E1, method = "greylevel", csize = 0.35, sub = "E1",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E2, method = "greylevel", csize = 0.35, sub = "E2",
csub=2)

