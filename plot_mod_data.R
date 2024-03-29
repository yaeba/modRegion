#!/usr/bin/env Rscript
"
Contains functions to plot per-read modification data
"

suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))


plot_smoothed_read <- function(df, title=NULL) {
  contig <- unique(df$seqname)
  p <- df %>%
    mutate(read_id = paste(read_id, sample, sep='_')) %>%
    ggplot() +
    geom_smooth(method="loess",
                aes(x=pos, y=prob_mod, group=read_id, color=sample),
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
#       m = purrr::map(data, loess, formula=prob_mod ~ pos,
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
  
  p <- ggplot(df, aes(x=pos, group=read_id, y=1, color=prob_mod)) + 
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
    summarise(prob_mod = mean(prob_mod)) %>%
    ungroup() %>%
    ggplot(aes(x=pos, y=prob_mod, color=sample)) +
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
  
  p <- split(filtered, filtered$gene_name) %>%
    GRangesList() %>%
    ggplot()
  
  if (!isEmpty(filtered)) {
    p <- p + geom_alignment(label=TRUE)
  } else {
    p <- p + geom_point() + ggplot2::xlim(genome_start, genome_end)
  }
  p
}


plot_tracks <- function(df, gtf=NULL, title="TITLE", highlights=NULL, max_reads=50) {
  suppressMessages(library(ggbio))
  
#  if (dim(df)[1] > max_reads) {
#    n <- max_reads %/% length(unique(df$sample))
#    ## subsample from df
#    df <- df %>%
#      group_by(sample) %>%
#      subset(read_id %in% sample(unique(.$read_id), n)) %>%
#      ungroup()
#    
#  }
  
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
  
  p
}
