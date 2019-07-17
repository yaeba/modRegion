#!/usr/bin/env Rscript


'Main Program.

Usage:
	modRegion.R extract (-r=<region> | -f=<regions_txt>)... (-o=<output> | -p=<plot> [--title=<title>])... [--nanopolish <[label:]fileA>]... [--tombo <[label:]fileB>]... [--deepsignal <[label:]fileC>]...
	modRegion.R overlap (-o=<output> | -p=<plot>)... [--overhang <overhang>] [--gene_name <gene_name>]... [--gene_biotype <gene_biotype>]... [-r=<region> | -f=<regions_txt>]... [--nanopolish <[label:]fileA>]... [--tombo <[label:]fileB>]... [--deepsignal <[label:]fileC>]... <gtf_file>
	modRegion.R --help

Options:
	-h --help  Show this screen.
	-r --region=<region>  Genomic region in UCSC/IGV/samtools format, eg "chr9:3500-4500", can be given multiple times.
	-f --regions_file=<regions_txt>  File containing genomic regions (one per line).
	--nanopolish <[label:]fileA>  Nanopolish methylation data and (optional) label name.
	--tombo <[label:]fileB>  Tombo methylation data and (optional) label name.
	--deepsignal <[label:]fileC>  DeepSignal methylation data and (optional) label name.
	-o --output=<output>  Save dataframe to .tsv or .tsv.gz file.
	-p --plot=<plot>  Save plot(s) to pdf.
	--gene_name <gene_name>  Name of the gene, can be given multiple times. 
	--gene_biotype <gene_biotype>  Biotype of gene, eg "protein_coding", can be given multiple times. [default: protein_coding]
	--overhang <overhang>  Up to how many bases up and downstream from gene. [default: 2000]
	--title=<title>  Title of the plot.
	
Arguments:
	gtf_file  Gene Transfer Format that contains gene informations
' -> doc

suppressMessages(library(docopt))
suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer))
suppressMessages(library(data.table))

source("load_mod_data.R")
source("plot_mod_data.R")
source("mod_gene_overlaps.R")

arguments <- docopt(doc)

print(arguments)

CALLERS <- c("nanopolish", "tombo", "deepsignal")

# main.R overlap [-r <region] (-o <output> | -p <plot>)... --gene <gtf_file> [--nanopolish <[label:]fileA>]... [--tombo <[label:]fileB>]... [--deepsignal <[label:]fileC]...


verbose <- function(...) cat(sprintf(...), sep='', file=stderr())


parse_label_file <- function(string) {
  count <- attr(parse_label_file, "running_count")
  
  if(is.null(count)) {
     count <- 0
  }
  count <- count + 1

  # arg could be "file.tsv" or "label:file.tsv"
  splitted <- unlist(strsplit(string, ':'))
  filename <- splitted
  if (length(splitted) == 1) {
    label <- paste0("sample_", count)
  } else {
    label <- splitted[1]
    filename <- splitted[2]
  }
  
  attr(parse_label_file, "running_count") <<- count
  
  
  list(label=label, filename=filename)
}


write_output <- function(df, outfile) {
  df %>%
    mutate(log_lik_ratio = format(log_lik_ratio, digits=4, format='f', scientific=FALSE),
           prob_meth = format(prob_meth, digits=4, format='f', scientific=FALSE)) %>%
    vroom_write(outfile)
}


plot_region <- function(df, raw_regions, outfile, width=10, height=10, title="Methylation Pattern") {
  options(warn=-1)
  
  pdf(outfile)
  
  for (r in raw_regions) {
    reg <- parse_regions(r)
    reg_title <- paste(title, r, sep=", ")
    reg_df <- filter_region(df, reg)
    print(plot_read(reg_df, reg_title))
    print(plot_smoothed_read(reg_df, reg_title))
    print(plot_aggregated_read(reg_df, reg_title))
    
  }
  # p <- df %>%
  #   plot_read()
  # ggsave(outfile, plot=p, width=width, height=height)
  invisible(dev.off())
}


extract_regions <- function(arguments, raw_regions) {
  dfs <- list()
  
  ## read all input files
  for (caller in CALLERS) {
    for (arg in arguments[[caller]]) {
      parsed <- parse_label_file(arg)
      verbose("Loading %s from %s..\n", parsed$label, parsed$filename)
      
      dfs[[parsed$label]] <- load_mod_data(parsed$filename, caller, raw_regions=raw_regions)
    }
  }
  
  ## merge the queried dataframes
  bind_rows(dfs, .id="sample")
}



overlap_genes <- function(arguments, raw_regions) {
  dfs <- list()
  
  verbose("Loding genes data from %s..\n", arguments$gtf_file)
  ## read genes data
  genes <- load_gtf(arguments$gtf_file, 
                    arguments$gene_name, arguments$gene_biotype,
                    as.numeric(arguments$overhang))
  
  ## find overlaps for all input files
  for (caller in CALLERS) {
    for (arg in arguments[[caller]]) {
      parsed <- parse_label_file(arg)
      verbose("Loading %s from %s..\n", parsed$label, parsed$filename)
      
      dfs[[parsed$label]] <- mod_gene_overlaps(parsed$filename, caller, genes,
                                               raw_regions=raw_regions)
    }
  }
  
  ## merge all overlap dataframes
  bind_rows(dfs, .id="sample")
  
  
  # b <- q[seqnames(q)=='10' & start(q) >  69398773 & end(q) <        70027438]
  # 
  # q <- a %>% as_tibble() %>% filter(type %in% c('exon')) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  # b <- q[seqnames(q)=='10' & start(q) >=  69422683 & end(q) <=        69566361]
  #         
  # grlist <- split(b, b$gene_name) %>% GRangesList()
  # p <- ggplot(grlist) + geom_alignment(label=TRUE)
  
}

main <- function(arguments) {

  ## parse all the genomic regions
  from_file <- NULL
  if (!is.null(arguments$regions_file)) {
    from_file <- read.table(arguments$regions_file, sep='\n') %>%
      pull('V1') %>%
      as.vector()
  }
  raw_regions <- c(from_file, arguments$region)
  
  verbose("Region(s): ")
  verbose(paste(raw_regions, collapse=', '))
  verbose('\n')
  
  df <- list()
  
  if (arguments$extract) {
    df <- extract_regions(arguments, raw_regions)
  }
  
  if (arguments$overlap) {
    
    
    # ## not yet implemented
    # verbose("Not yet implemented\n")
    # quit(save='no')
    
    df <- overlap_genes(arguments, raw_regions)
  }
  
  
  ## output merged dataframe
  if (!is.null(arguments$output)) {
    outfile <- arguments$output
    verbose("Writing to %s..\n", outfile)
    write_output(df, outfile)
  }
  
  ## plot and save to pdf
  if (!is.null(arguments$plot) && !arguments$overlap) {
    outfile <- arguments$plot
    verbose("Plotting to %s..\n", outfile)
    plot_region(df, raw_regions, outfile)
  }
  
}

main(arguments)
