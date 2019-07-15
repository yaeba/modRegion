#!/usr/bin/env Rscript
"
Contains functions to plot per-read modification data
"

suppressMessages(library(ggplot2))


# plot_smoothed_read <- function(df) {
#   contig <- unique(df$seqname)
#   p <- ggplot(df) +
#     geom_smooth(method="loess",
#                 aes(x=pos, y=prob_meth, group=read_id, color=sample),
#                 se=FALSE) +
#     ggplot2::ylim(0, 1) +
#     geom_point(aes(x=pos), y=0, pch='|') +
#     theme_bw() +
#     ggtitle(title) +
#     xlab(paste(contig, "pos", sep=", "))
#   p
# }


plot_smoothed_read <- function(df, title) {
  contig <- unique(df$seqname)
  smoothed <- df %>%
    tidyr::nest(-read_id) %>%
    dplyr::mutate(
      m = purrr::map(data, loess, formula=prob_meth ~ pos,
                     method.args=list(span=0.1 + 8 * 10^(-11) * (max(10^5 - (max(pos) - min(pos)), 0)^2))),
      fitted = purrr::map(m, `[[`, "fitted")
    ) %>%
    select(-m) %>%
    tidyr::unnest()
  
  p <- ggplot(smoothed) +
    geom_line(aes(x=pos, y=fitted, group=read_id, color=sample), size=1) +
    ggplot2::ylim(0, 1) +
    geom_point(aes(x=pos), y=0, pch='|') +
    theme_bw() +
    ggtitle(title) +
    xlab(paste(contig, "pos", sep=", "))
  p
}


plot_read <- function(df, title) {
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

plot_aggregated_read <- function(df, title) {
  contig <- unique(df$seqname)
  p <- df %>%
    group_by(sample, pos) %>%
    summarise(prob_meth = mean(prob_meth)) %>%
    ungroup() %>%
    ggplot(aes(x=pos, y=prob_meth)) +
    geom_point(alpha=0.1) +
    geom_smooth(method="loess", aes(color=sample)) +
    coord_cartesian(ylim=c(0, 1)) +
    geom_point(aes(x=pos), y=0, pch='|') +
    theme_bw() +
    ggtitle(title) +
    xlab(paste(contig, "pos", se=", "))
  p
}
