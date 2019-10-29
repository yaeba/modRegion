#!/usr/bin/env Rscript
## Script to load aggregated nanopolish methylation data (seqname, pos, strand, llr, prob) 
## and reshape to (chr, pos, forward, forward_prob, reverse_prob)
## Usage: ./reshape_nanopolish_strand.R  <path>  <pattern>  <output.rds>

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


args = commandArgs(trailingOnly=TRUE)


load_aggregated <- function(infile) {
  cat(paste("Loading", infile, "\n"), file=stdout())
  readRDS(infile) %>%
    dplyr::rename(chr = seqname,
    		  prob_meth = prob_mod)
}

reshape_strand <- function(df) {
  cat("Reshaping\n", file=stdout())
  df %>%
    select(-log_lik_ratio) %>%
    group_by(chr, pos) %>%
    spread(key=strand, value=prob_meth) %>%
    dplyr::rename(forward_prob = '+', reverse_prob = '-') %>%
    ungroup() %>%
    arrange(chr, pos)
}




main <- function(args) {
  if (length(args) < 3) {
    cat("Wrong arguments supplied\n", file=stderr())
    quit(save='no')
  }
  
  path <- args[1]
  pattern <- args[2]
  output <- args[3]

  if (!grepl(".rds$", output)) {
    output <- paste0(output, ".rds")
  }

  files <- list.files(path, pattern, full.names=TRUE)

  cat(paste("Found", paste(files, collapse=' '), '\n'), file=stdout())
  
  files %>%
    lapply(function(x) {
      load_aggregated(x) %>%
      reshape_strand() 
    }) %>%
    rbindlist() %>%
    saveRDS(output)
}


main(args)