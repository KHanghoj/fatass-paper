require(ggplot2)
require(reshape2)
args <- commandArgs(trailing = TRUE)
df <- read.table(args[1])
pop <- rep(scan(args[2], what = "thefuck"), each = 2)
d2 <- df / rowSums(df)
a <- seq(0, 1000)
good_sub <- a[bitwAnd(a, a - 1) == 0 & a != 1]
colnames(d2) <- paste0("sub", good_sub[seq_len(ncol(d2))])
bitmap(args[3], res = 400, width = 9, height = 3)
ggplot(melt(d2), aes(value, fill = variable)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Distribution for each haplotype") +
  labs(x = "freq")
dev.off()

d3 <- data.frame(pop = pop, d2)
d3$group <- seq_len(nrow(d3))

bitmap(args[4], res = 400, width = 10, height = 8)
ggplot(melt(d3, id.vars = c("pop", "group")), aes(variable, value)) +
  geom_point() +
  geom_line(aes(group = group, col = pop)) +
  facet_wrap(~pop) +
  theme_bw() +
  labs(x = NULL, y = "Freq of windows")
dev.off()