
require(ggplot2)
args <- commandArgs(trailing=T)
eigenvecs <- args[1]
labels <- args[2]
names <- args[3]
pc1 <- args[4]
pc2 <- args[5]
out <- args[6]
df <- read.table(eigenvecs)
colnames(df) <- paste0("PC", 1:ncol(df))
labels <- as.character(read.table(labels)[,1])
pops <- substr(labels, 1, nchar(labels)-2)
ks <- substr(labels, nchar(labels)-1, nchar(labels))

names <- as.character(read.table(names)[,1])
names <- substr(names, 1, nchar(names)-2)
namesAdm <- names[duplicated(names)]
Admixed <- rep("Unadmixed", length(names))
Admixed[names %in% namesAdm] <- "Admixed"

plotdf <- cbind(df[,c(pc1,pc2)], p=pops, k=ks, Admixed=factor(Admixed))

df_unadm <- subset(plotdf, Admixed == "Unadmixed")
df_adm <- subset(plotdf, Admixed == "Admixed")

# bitmap(out, res=500, width=7, height=5)
p <- ggplot() +
  scale_size_manual(values = c("Unadmixed" = .75, "Admixed" = 1.25)) +
  scale_color_brewer(palette = "Paired") +
  # scale_alpha_manual(values = c("Unadmixed" = .2, "Admixed" = .5)) +
  guides(col = guide_legend(ncol = 2, byrow = TRUE),
         shape = guide_legend(ncol = 2, byrow = TRUE)) +
  geom_point(data = df_adm, aes(.data[[pc1]], .data[[pc2]],
                      col = p, shape = k,
                      size = Admixed)) +
  geom_point(data = df_unadm, aes(.data[[pc1]], .data[[pc2]],
                      col = p, shape = k,
                      size = Admixed)) +
  theme_bw()
# # bitmap(out, res=500, width=7, height=5)
# p<-ggplot(plotdf, aes(.data[[pc1]], .data[[pc2]], 
#                       col = p, shape = k, 
#                       size=Admixed)) + 
#                       # alpha=Admixed)) + 
#   scale_size_manual(values = c("Unadmixed" = .5, "Admixed" = 1)) +
#   scale_color_brewer(palette = "Paired") +
#   # scale_alpha_manual(values = c("Unadmixed" = .2, "Admixed" = .5)) +
#   guides(col = guide_legend(ncol = 2, byrow = TRUE),
#          shape = guide_legend(ncol = 2, byrow = TRUE)) +
#   geom_point() +
#   theme_bw()
# dev.off()
ggsave(out, p, width = 8, height = 5, dpi = 350)