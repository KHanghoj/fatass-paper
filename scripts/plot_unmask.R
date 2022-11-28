
require(ggplot2)
args <- commandArgs(trailing=T)
eigenvecs <- args[1]
labels <- args[2]
out <- args[3]
df <- read.table(eigenvecs)
pops <- as.character(read.table(labels)[,1])
plotdf <- cbind(df[,1:3], pops)
bitmap(out, res=500)
ggplot(plotdf, aes(V1,V2, col=pops)) + geom_point() + theme_bw()
dev.off()
