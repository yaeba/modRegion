#!/usr/bin/env Rscript
"
Contains functions to plot per-read modification data
"

suppressMessages(library(ggplot2))


plot_smoothed_read <- function(df, title="Title") {
  contig <- unique(df$seqname)
  p <- ggplot(df) +
    geom_smooth(method="loess",
                aes(x=pos, y=prob_meth, group=read_id, color=sample),
                se=FALSE) +
    ggplot2::ylim(0, 1) +
    geom_point(aes(x=pos), y=0, pch='|') +
    theme_bw() +
    ggtitle(title) +
    xlab(paste(contig, "pos", sep=", "))
  p
}

plot_read <- function(df, title="Title") {
  contig <- unique(df$seqname)
  p <- ggplot(df, aes(x=pos, group=read_id, y=1, color=prob_meth)) + 
    geom_line(position=ggstance::position_dodgev(height=0.2), size=3) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0.5, space="Lab") +
    facet_grid(rows=vars(sample)) +
    theme_bw() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggtitle(title) +
    xlab(paste(contig, "pos", sep=", "))
  p
}