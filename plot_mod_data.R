#!/usr/bin/env Rscript
"
Contains functions to plot per-read modification data
"

suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))


plot_smoothed_read <- function(df, title=NULL) {
  contig <- unique(df$seqname)
  p <- ggplot(df) +
    geom_smooth(method="loess",
                aes(x=pos, y=prob_meth, group=read_id, color=sample),
                se=FALSE, span=1) +
    ggplot2::ylim(0, 1) +
    geom_point(aes(x=pos), y=0, pch='|') +
    theme_bw() +
    ggtitle(title) +
    xlab(paste0("Contig: ", contig, ", pos"))
  p
}


# plot_smoothed_read <- function(df, title) {
#   contig <- unique(df$seqname)
#   smoothed <- df %>%
#     tidyr::nest(-read_id) %>%
#     dplyr::mutate(
#       m = purrr::map(data, loess, formula=prob_meth ~ pos,
#                      method.args=list(span=0.1 + 8 * 10^(-11) * (max(10^5 - (max(pos) - min(pos)), 0)^2))),
#       fitted = purrr::map(m, `[[`, "fitted")
#     ) %>%
#     select(-m) %>%
#     tidyr::unnest()
#   
#   p <- ggplot(smoothed) +
#     geom_line(aes(x=pos, y=fitted, group=read_id, color=sample), size=1) +
#     ggplot2::ylim(0, 1) +
#     geom_point(aes(x=pos), y=0, pch='|') +
#     theme_bw() +
#     ggtitle(title) +
#     xlab(paste(contig, "pos", sep=", "))
#   p
# }


plot_read <- function(df, title=NULL) {
  contig <- unique(df$seqname)
  
  p <- ggplot(df, aes(x=pos, group=read_id, y=1, color=prob_meth)) + 
    geom_line(position=ggstance::position_dodgev(height=0.2), size=3) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0.5, space="Lab") +
    facet_grid(rows=vars(sample), switch='y') +
    theme_bw() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(x=paste0("Contig: ", contig, ", pos"), y="read", title=title)
  p
}

plot_aggregated_read <- function(df, title=NULL) {
  contig <- unique(df$seqname)
  p <- df %>%
    group_by(sample, pos) %>%
    summarise(prob_meth = mean(prob_meth)) %>%
    ungroup() %>%
    ggplot(aes(x=pos, y=prob_meth, color=sample)) +
    geom_point(alpha=0.1) +
    geom_smooth(method="loess") +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    ggtitle(title) +
    xlab(paste0("Contig: ", contig, ", pos"))
  p
}


plot_annotation <- function(gtf, seqname, genome_start, genome_end) {
  
  filtered <- gtf[seqnames(gtf) == seqname & start(gtf) >= genome_start & end(gtf) <= genome_end] %>%
    as_tibble() %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  
  split(filtered, filtered$gene_name) %>%
    GRangesList() %>%
    ggplot() +
    geom_alignment(label=TRUE)
}


plot_tracks <- function(df, gtf=NULL, title="TITLE", highlights=NULL) {
  suppressMessages(library(ggbio))
  
  contig <- unique(df$seqname)
  genome_start <- min(df$pos)
  genome_end <- max(df$pos)
  
  
  plots <- list()
  
  plots$Reads <- plot_read(df)
  
  plots$Smoothed <- plot_smoothed_read(df)
  
  plots$Mean <- plot_aggregated_read(df)
  
  if (!is.null(gtf)) {
    plots$Genes <- gtf %>%
      subset(type %in% "exon") %>%
      plot_annotation(contig, genome_start, genome_end)
  }
  
  p <- tracks(plots, main=title) +
    theme_tracks_sunset()
  
  invisible(p)
}
