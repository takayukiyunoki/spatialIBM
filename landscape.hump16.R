# For small system and 16 humps
n.sites = 400
landscape <- data.frame(sites = 1:n.sites, expand.grid(x = 1:sqrt(n.sites), y = 1:sqrt(n.sites)))
hump4 <- matrix(NA, nrow = sqrt(n.sites/4)/2, ncol = sqrt(n.sites/4)/2)
for (i in 1: sqrt(n.sites/4)/2) for (j in 1: sqrt(n.sites/4)/2) 
hump4[i, j] <- dist(rbind(c(i, j), rep((sqrt(n.sites/4)/2 + 1)/2, 2)))
hranks <- rank(hump4)
breaks <- rep((n.sites/40)/4, 10)
breaks <- cumsum(breaks)
hump4[hranks <= breaks[1]] <- 1
for (i in 2:10) 
hump4[hranks > breaks[i - 1] & hranks <= breaks[i]] <- i
hump4 <- cbind(hump4, hump4)
hump4 <- rbind(hump4, hump4)
hump8 <- matrix(rep(hump4,2), ncol=ncol(hump4), byrow=TRUE)
hump16 <- matrix(rep(hump8,2), ncol=ncol(hump8)*2, byrow=TRUE)
hump16 <- as.vector(hump16)/10
landscape$E <- hump16
write.csv(landscape, file = "landscape.hump16.csv")
landscape.hump16 <- landscape
par(mfrow = c(4,1))
sr.value(landscape.hump16[, c("x", "y")], landscape.hump16$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
landscape.center <- landscape.hump16[landscape.hump16$x >= 6 & landscape.hump16$x <= 15 & landscape.hump16$y >= 6 & landscape.hump16$y <= 15, ]
E.trans <- poly(landscape.center$E, 2)
colnames(E.trans) <- c("E1", "E2")
rownames(E.trans) <- rownames(landscape.center) 
E.trans <- as.data.frame(E.trans)
write.csv(E.trans, file = "E.trans.hump16.center.csv")
sr.value(landscape.center[, c("x", "y")], landscape.center$E, method = "greylevel", csize = 0.35, sub = "E",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E1, method = "greylevel", csize = 0.35, sub = "E1",
csub=2)
sr.value(landscape.center[, c("x", "y")], E.trans$E2, method = "greylevel", csize = 0.35, sub = "E2",
csub=2)

