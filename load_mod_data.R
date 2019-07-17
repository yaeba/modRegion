#!/usr/bin/env Rscript
"
Contains functions that accepts data output from base modification callers, selects only 
required rows/columns and returns an indexed dataframe for lazy access.
"

suppressMessages(library(vroom))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))

source("mod_nanopolish.R")
source("mod_tombo.R")
source("mod_deepsignal.R")

ORDER_EXTRACT <- 
  c("seqname", "pos", "read_id", "strand", "log_lik_ratio", "prob_meth")


parse_genomic_region <- function(string) {
  splitted <- unlist(strsplit(string, ':'))
  reg_pos <- as.numeric(gsub(',', '', unlist(strsplit(splitted[2], '-'))))
  
  if (length(splitted) != 2 || length(reg_pos) != 2) {
    stop(paste("Cannot parse region", string))
    quit(save='n')
  }
  
  region <- list()
  region$contig <- splitted[1]
  
  ## region strings are 1-based
  region$start <- reg_pos[1] - 1
  region$end <- reg_pos[2] - 1
  
  region
}

parse_regions <- function(raw_regions) {
  lapply(raw_regions, parse_genomic_region)
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


load_mod_data <- function(filename, caller, order=ORDER_EXTRACT, raw_regions=NULL, motif="CG") {
  data_file <- add_class(filename, caller)
  
  # load
  data <- load_file(data_file)
  
  # filter
  parsed_regions <- parse_regions(raw_regions)
  data <- filter_region(data, parsed_regions)
  
  # preprocess
  data <- add_class(data, caller)
  data <- preprocess(data, motif)

  data %>%
    select(order)
}

add_class <- function(x, classname) structure(x, class=append(class(x), classname))

load_file <- function(filename) UseMethod("load_file")

preprocess <- function(df, ...) UseMethod("preprocess")


# fileA <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/small_nanopolish.tsv.gz"
# fileB <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/small_tombo.tsv"
# fileC <- "/stornext/HPCScratch/home/tay.x/scripts/notebooks/small_deepsignal.tsv"
# a <- load_tombo(fileB, 'OCVW01001791.1:40000-50000')
# a <- load_file.tombo(fileB)
# b <- load_nanopolish(fileA, "OCVW01000001.1:26000-30000")
# b <- load_file.nanopolish(fileA)
# c <- load_deepsignal(fileC, "NC_001144.5:460000-470000")
# c <- load_file.deepsignal(fileC)
