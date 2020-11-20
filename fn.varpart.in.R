# Function to calculating the summary statics of variation partitioning for inicial conditions
fn.varpart.in <- function(PCNM.pos, E.trans, pool, pool.t0) {
varpart.in <- data.frame(matrix(, nrow=1, ncol=0))
varpart.in$n.groups.in <- length(unique(pool.t0$groups))
pool$species <- pool.t0$species
LC <- dcast(pool, locations~species, length) [,-1] 
LC.h <- decostand(LC, "hellinger")
LC.PCNM.rda <- rda(LC.h, PCNM.pos)
anova.LC.PCNM.rda <- anova.cca(LC.PCNM.rda)
varpart.in$space.sign <- anova.LC.PCNM.rda$Pr[1]
LC.env.rda <- rda(LC.h, E.trans)
anova.LC.env.rda <- anova.cca(LC.env.rda)
varpart.in$env.sign <- anova.LC.env.rda$Pr[1]
return(varpart.in)
}
