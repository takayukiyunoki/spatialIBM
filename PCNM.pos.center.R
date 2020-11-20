# For small system
n.sites = 400
landscape <- data.frame(sites = 1:n.sites, expand.grid(x = 1:sqrt(n.sites), y = 1:sqrt(n.sites)))
landscape.center <- landscape[landscape$x >= 6 & landscape$x <= 15 & landscape$y >= 6 & landscape$y <= 15, ]
dist.landscape <- dist(landscape.center[, c("x", "y")])
pcnm <- pcnm(dist.landscape)
PCNM.overall <- PCNM(dist.landscape)
length(which(pcnm$values>0))
select <- which(PCNM.overall$Moran_I$Positive==TRUE)
length(select)
par(mfrow = c(3,2))
some.pcnm <- c(1, 2, 5, 10, 20, 50)
for(i in 1:length(some.pcnm)) {
sr.value(landscape.center[, c("x", "y")], PCNM.overall$vectors[, some.pcnm[i]], method = "greylevel", csize = 0.35, sub = some.pcnm[i], csub=2)
}
PCNM.pos.center <- as.data.frame(PCNM.overall$vectors)[,select]
write.csv(PCNM.pos.center, file = "PCNM.pos.center.csv")
