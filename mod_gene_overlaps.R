#!/usr/bin/env Rscript
"
Contains functions to read and find overlaps between modification sites and
gene data in GTF file
"

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

ORDER_OVERLAP <-
  c("seqname", "pos", "read_id", "strand", "log_lik_ratio", "prob_meth",
    "gene_name", "gene_biotype", "gene_strand", "gene_start", "gene_end")

load_genes <- function(gtf, gene_names=NULL, gene_biotypes=NULL, overhang=2000) {

  genes <- gtf %>%
    as_tibble() %>%
    dplyr::rename(seqname = seqnames,
                  gene_strand = strand) %>%
    filter(type == "gene",
           if (is.null(gene_names)) TRUE else gene_name %in% gene_names,
           if (is.null(gene_biotypes)) TRUE else gene_biotype %in% gene_biotypes) %>%
    group_by(gene_name, seqname, gene_biotype, gene_strand) %>%
    summarise(gene_start = min(start) - 1,
              gene_end = max(end) - 1) %>%
    mutate(start = gene_start - overhang,
           end = gene_end + overhang) %>%
    
    # tss = ifelse(gene_strand == "+", start, end),
    # tts = ifelse(gene_strand == "+", end, start)) %>%
    as.data.table()
  
  genes
}


find_overlaps <- function(mod_df, genes) {
  setkey(genes, seqname, start, end)
  
  ## not optimised for large dataset in terms of memory
  mod_table <- mod_df %>%
    dplyr::rename(start = pos) %>%
    as.data.table()
  mod_table[, end := start]
  
  overlaps <- foverlaps(mod_table, genes, nomatch=0, type='within')
  setnames(overlaps,
           old="i.start",
           new="pos")
  
  overlaps %>%
    as_tibble() %>%
    select(-i.end, -start, -end)
}


mod_gene_overlaps <- function(filename, caller, genes, order=ORDER_OVERLAP,
                              raw_regions=NULL, motif="CG") {
  source("load_mod_data.R")
  
  ## load modification data
  mod_df <- filename %>%
    add_class(caller) %>%
    load_file()
  
  ## filter if regions are present
  if (!is.null(raw_regions)) {
    mod_df <- mod_df %>%
      filter_region(raw_regions %>%
                      parse_regions())
  }
  
  ## find overlaps with genes
  mod_overlaps <- find_overlaps(mod_df, genes)
  
  ## preprocess remaining columns
  mod_overlaps %>%
    add_class(caller) %>%
    preprocess(motif) %>%
    select(order)
  
}
