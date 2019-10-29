#!/usr/bin/env Rscript
## Script to load all nanopolish chunk files and aggregate to (seqname, pos, strand) level
## Usage: ./aggregated_nanopolish.R  <path>  <pattern>  <output.rds>

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
source("../load_mod_data.R", chdir=T)

args = commandArgs(trailingOnly=TRUE)


load_partial_nanopolish <- function(infile) {
  cat(paste("Loading", infile, "\n"), file=stdout())
  load_file.nanopolish(infile) %>%
    preprocess.nanopolish() %>%
    group_by(seqname, pos, strand) %>%
    summarise(sum_llr = sum(log_lik_ratio),
              n = n()) %>%
    ungroup()
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
    lapply(load_partial_nanopolish) %>%
    rbindlist() %>%
    group_by(seqname, pos, strand) %>%
    summarise(log_lik_ratio = sum(sum_llr) / sum(n),
              prob_mod = 1 / (1 + exp(log_lik_ratio))) %>%
    ungroup() %>%
    saveRDS(output)
}


main(args)