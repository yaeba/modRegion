#!/usr/bin/env Rscript
"
Contains functions that accepts data output from base modification callers, selects only 
required rows/columns and returns an indexed dataframe for lazy access.
"

suppressMessages(library(vroom))
suppressMessages(library(dplyr))

ORDER <- c("seqname", "pos", "read_id", "strand", "log_lik_ratio", "prob_meth")

lazy_load <- function(filename, cols, coltypes, colnames=NULL) {
  if (is.null(colnames)) {
    vroom(filename, col_select=cols, col_types=coltypes)
  } else {
    vroom(filename, col_select=cols, col_types=coltypes, col_names=colnames)
  }
}

maybe_filter <- function(df, contig, start, end) {
  df <- df %>%
    filter(
      if (is.null(contig)) TRUE else seqname == contig,
      if (is.null(start)) TRUE else pos >= start,
      if (is.null(end)) TRUE else pos <= end
    )
  df
}

load_tombo <- function(filename, contig=NULL, start=NULL, end=NULL) {
  
  cols <- c("chrm", "pos", "read_id", "strand", "stat")
  coltypes <- list(chrm = 'c', pos = 'i', read_id = 'c', stat = 'd',
                   strand = col_factor(levels=c('+', '-')))
  
  data <- lazy_load(filename, cols, coltypes) %>%
    dplyr::rename(seqname = chrm,
                  log_lik_ratio = stat) %>%
    maybe_filter(contig, start, end) %>%
    mutate(prob_meth = 1 / (1 + exp(log_lik_ratio)),
           pos = ifelse(strand == "-", pos - 1, pos)) %>%
    select(ORDER)
  
  data
}

load_nanopolish <- function(filename, contig=NULL, start=NULL, end=NULL, motif="CG") {
  
  cols <- c("chromosome", "start", "read_name", "strand", "log_lik_ratio", 
            "num_motifs", "sequence")
  coltypes <- list(chromosome = 'c', start = 'i', read_name = 'c', 
                   log_lik_ratio = 'd', num_motifs = 'i', sequence = 'c',
                   strand = col_factor(levels=c('+', '-')))
  
  data <- lazy_load(filename, cols, coltypes) %>%
    dplyr::rename(seqname = chromosome,
                  pos = start,
                  read_id = read_name) %>%
    maybe_filter(contig, start, end) %>%
    mutate(log_lik_ratio = -1 * log_lik_ratio)
  
  ## Assign same log-lik-ratio to all bases covered in motif
  offset <- unlist(gregexpr(pattern=motif, data$sequence)) - 1
  data <- data %>%
    uncount(num_motifs) %>%
    mutate(pos = pos - 5 + offset) %>%
    select(-sequence)
  
  ## Upon inspecting, a single read sometimes produces more than one statistic at a position
  ## Assuming statistic in the middle contribute the most, use weights from binomial distribution
  data <- data %>%
    group_by(seqname, pos, read_id, strand) %>%
    summarise(log_lik_ratio = sum(
      dbinom(0:(n()-1), n()-1, 0.5) * log_lik_ratio
    )) %>%
    ungroup() %>%
    mutate(prob_meth = 1 / (1 + exp(log_lik_ratio))) %>%
    select(ORDER)
  
  data
}


load_deepsignal <- function(filename, contig=NULL, start=NULL, end=NULL) {
  
  colnames <- c("seqname", "pos", "strand", "x1", 
            "read_id", "x2", "x3", "prob_meth", "x4", "x5")
  cols <- c("seqname", "pos", "read_id", "strand", "prob_meth")
  coltypes <- list(seqname = 'c', pos = 'i', read_id = 'c', prob_meth = 'd',
                   strand = col_factor(levels=c('+', '-')))
  
  data <- lazy_load(filename, cols, coltypes, colnames=colnames) %>%
    maybe_filter(contig, start, end) %>%
    mutate(log_lik_ratio = log((1 - prob_meth) / prob_meth)) %>%
    select(ORDER)
  
  data
}


#fileA <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/large_nanopolish.tsv.gz"
# fileB <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/small_tombo.tsv"
# fileC <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/small_deepsignal.tsv"
# a <- load_tombo(fileB, 'tig00000080|arrow|arrow', 460000, 470000)
#b <- load_nanopolish(fileA, "OCVW01000001.1", 26000, 28000)
# c <- load_deepsignal(fileC, "NC_001144.5", 460000, 470000)
