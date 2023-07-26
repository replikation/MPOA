#!/usr/bin/env Rscript

library(ggplot2)
require(readr)

names <- c('deg', 'base', 'ratio')
df1 <- read_tsv("input.tsv", col_names=names)

plot <- ggplot(df1, aes(x=deg, y=ratio, fill=deg)) +
  scale_fill_brewer(palette="Blues") +
  geom_violin(alpha=0.4, position = position_dodge(width = .75),show.legend = FALSE) +
  coord_cartesian(ylim=c(0,1)) +
  #stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
  #stat_summary(fun.y=median, geom="point", size=2, color="red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color='black') +
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1, show.legend = FALSE) +
  ylab(  c("Base ratio")  ) +
  xlab(  c("masked Base")  ) +
  theme_classic()



svg("chart.svg")  # , height=as.numeric(length(combined))*2
plot
dev.off()
