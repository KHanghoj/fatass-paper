require(ggplot2)
require(reshape2)
df <- read.table("test.txt.log")
d2 <- df / rowSums(df)
bitmap("plot.png", res = 400, width = 6, height = 2)
ggplot(melt(d2), aes(value, fill = variable)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Distribution for each haplotype with min 200 loglike units improvement") +
  labs(x = "freq")
dev.off()