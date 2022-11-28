
require(ggplot2)
args <- commandArgs(trailing=T)
eigenvecs <- args[1]
labels <- args[2]
pc1 <- args[3]
pc2 <- args[4]
out <- args[5]
df <- read.table(eigenvecs)
colnames(df) <- paste0("PC", 1:ncol(df))
labels <- as.character(read.table(labels)[,1])
pops <- substr(labels, 1, nchar(labels)-2)
ks <- substr(labels, nchar(labels)-1, nchar(labels))
plotdf <- cbind(df[,c(pc1,pc2)], pops, ks)
bitmap(out, res=500, width=7, height=5)
ggplot(plotdf, aes(.data[[pc1]], .data[[pc2]], col = pops, shape = ks)) + 
  geom_point() + 
  theme_bw()
dev.off()
