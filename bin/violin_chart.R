#!/usr/bin/env Rscript

library(ggplot2)
require(readr)

names <- c('deg', 'Bases', 'ratio')
df1 <- read_tsv("input.tsv", col_names=names)

plot <- ggplot(df1, aes(x=Bases, y=ratio, fill=Bases)) +
  scale_fill_brewer(palette="Set2") +
  geom_violin(alpha=0.2, position = position_dodge(width = .75)) +
  coord_cartesian(ylim=c(0.2,0.8)) +
  #stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
  #stat_summary(fun.y=median, geom="point", size=2, color="red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color='black') +
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.2, show.legend = FALSE) +
  facet_grid(cols = vars(deg),scales = "free", space = "free") +
  ylab(  c("Base fraction for each affected position")  ) +
  xlab(  c("Observed types of base combination - labeled by their respective IUPAC code")  ) +
  theme_classic() +
  theme(legend.position="bottom") 



svg("frequency.svg", height=5, width=10)  # , height=as.numeric(length(combined))*2
plot
dev.off()
