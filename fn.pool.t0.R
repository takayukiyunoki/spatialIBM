#Internal function to generating the population sizes of guilds
fn.g.sizes <- function(JM, g) {
if (g ==1) {
g.sizes <- JM
} else {
vec <- sort(sample(x=1:(JM-1), size=g-1, replace = FALSE))
vec[g] <- JM
g.sizes <- c()
g.sizes[1] <- vec[1]
for (i in 2:g) {
g.sizes[i] <- vec[i]- vec[i-1]
}
}
return(g.sizes)
}
#Function to generating the initial conditions of regional species pool
fn.pool.t0 <- function(JM, g, diversity, n.samples.g.sizes=5) {
traits.mean <- runif(g, 0, 1)
traits.sd <- runif(g, 0, 10)
if (n.samples.g.sizes > 1) {
pool.t0 <- list()
set.seed(1) 
g.sizes <- replicate(n.samples.g.sizes, {fn.g.sizes(JM, g)})
for (i in 1:n.samples.g.sizes) {
if (diversity > 0) {
pool.t0[[i]] <- data.frame(row.names = 1:JM, species=1:JM, groups = rep(1:g, g.sizes[,i]), traits.mean = rep(traits.mean, g.sizes[,i]), traits.sd = rep(traits.sd, g.sizes[,i]))
} else {
pool.t0[[i]] <- data.frame(row.names = 1:JM, species = rep(1:g, g.sizes[,i]), groups = rep(1:g, g.sizes[,i]), traits.mean = rep(traits.mean, g.sizes[,i]), traits.sd = rep(traits.sd, g.sizes[,i]))
}
}
} else {
g.sizes <- fn.g.sizes(JM, g)
if (diversity > 0) {
pool.t0 <- data.frame(row.names = 1:JM, species=1:JM, groups = rep(1:g, g.sizes), traits.mean = rep(traits.mean, g.sizes), traits.sd = rep(traits.sd, g.sizes))
} else {
pool.t0 <- data.frame(row.names = 1:JM, species = rep(1:g, g.sizes), groups = rep(1:g, g.sizes), traits.mean = rep(traits.mean, g.sizes), traits.sd = rep(traits.sd, g.sizes))
}
}
return(pool.t0)
}
