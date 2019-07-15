#!/usr/bin/env Rscript
"
Contains functions that accepts data output from base modification callers, selects only 
required rows/columns and returns an indexed dataframe for lazy access.
"

suppressMessages(library(vroom))
suppressMessages(library(dplyr))

ORDER <- c("seqname", "pos", "read_id", "strand", "log_lik_ratio", "prob_meth")


parse_genomic_region <- function(string) {
  splitted <- unlist(strsplit(string, ':'))
  reg_pos <- as.numeric(gsub(',', '', unlist(strsplit(splitted[2], '-'))))
  
  if (length(splitted) != 2 || length(reg_pos) != 2) {
    stop(paste("Cannot parse region", string))
    quit(save='n')
  }
  
  region <- list()
  region$contig <- splitted[1]
  region$start <- reg_pos[1]
  region$end <- reg_pos[2]
  
  region
}

parse_regions <- function(raw_regions) {
  lapply(raw_regions, parse_genomic_region)
}



lazy_load <- function(filename, cols, coltypes, colnames=NULL) {
  if (is.null(colnames)) {
    vroom(filename, col_select=cols, col_types=coltypes)
  } else {
    vroom(filename, col_select=cols, col_types=coltypes, col_names=colnames)
  }
}



filter_region <- function(df, regions) {
  if (is_empty(regions)) {
    df
  } else {
    ## filter to those in region(s)
    in_region <- function(s, p, r) {s == r$contig & p >= r$start & p <= r$end}

    df %>%
      filter(
        rowSums(
          sapply(regions, function(r, s, p) {in_region(s, p, r)}, s=.$seqname, p=.$pos)
        ) > 0L
      )
  }
}

load_tombo <- function(filename, raw_regions=NULL) {
  regions <- parse_regions(raw_regions)
  
  cols <- c("chrm", "pos", "read_id", "strand", "stat")
  coltypes <- list(chrm = 'c', pos = 'i', read_id = 'c', stat = 'd',
                   strand = col_factor(levels=c('+', '-')))
  
  data <- lazy_load(filename, cols, coltypes) %>%
    dplyr::rename(seqname = chrm,
                  log_lik_ratio = stat) %>%
    filter_region(regions) %>%
    mutate(prob_meth = 1 / (1 + exp(log_lik_ratio)),
           pos = ifelse(strand == "-", pos - 1, pos)) %>%
    select(ORDER)
  
  data
}

load_nanopolish <- function(filename, raw_regions=NULL, motif="CG") {
  regions <- parse_regions(raw_regions)
  
  cols <- c("chromosome", "start", "read_name", "strand", "log_lik_ratio", 
            "num_motifs", "sequence")
  coltypes <- list(chromosome = 'c', start = 'i', read_name = 'c', 
                   log_lik_ratio = 'd', num_motifs = 'i', sequence = 'c',
                   strand = col_factor(levels=c('+', '-')))
  
  data <- lazy_load(filename, cols, coltypes) %>%
    dplyr::rename(seqname = chromosome,
                  pos = start,
                  read_id = read_name) %>%
    filter_region(regions) %>%
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


load_deepsignal <- function(filename, raw_regions=NULL) {
  regions <- parse_regions(raw_regions)
  
  colnames <- c("seqname", "pos", "strand", "x1", 
            "read_id", "x2", "x3", "prob_meth", "x4", "x5")
  cols <- c("seqname", "pos", "read_id", "strand", "prob_meth")
  coltypes <- list(seqname = 'c', pos = 'i', read_id = 'c', prob_meth = 'd',
                   strand = col_factor(levels=c('+', '-')))
  
  data <- lazy_load(filename, cols, coltypes, colnames=colnames) %>%
    filter_region(regions) %>%
    mutate(log_lik_ratio = log((1 - prob_meth) / prob_meth)) %>%
    select(ORDER)
  
  data
}


# fileA <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/small_nanopolish.tsv.gz"
# fileB <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/small_tombo.tsv"
# fileC <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/small_deepsignal.tsv"
# a <- load_tombo(fileB, 'OCVW01001791.1:40000-50000')
# b <- load_nanopolish(fileA, "OCVW01000001.1:26000-30000")
# c <- load_deepsignal(fileC, "NC_001144.5:460000-470000")
