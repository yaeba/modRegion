#!/usr/bin/env Rscript


'Main Program.

Usage:
	modRegion.R extract (-r=<region> | -f=<regions_txt>)... (-o=<output> | -p=<plot> [--title=<title>])... [--nanopolish <[label:]fileA>]... [--tombo <[label:]fileB>]... [--deepsignal <[label:]fileC]...
	modRegion.R overlap [--nanopolish <[label:]fileA>]... [--tombo <[label:]fileB>]... [--deepsignal <[label:]fileC]...

Options:
	-h --help  Show this screen.
	-r --region=<region>  Genomic region in UCSC/IGV/samtools format, eg "chr9:3500-4500", can be given multiple times.
	-f --regions_file=<regions_txt>  File containing genomic regions (one per line).
	--nanopolish <[label:]fileA>  Nanopolish methylation data and (optional) label name.
	--tombo <[label:]fileB>  Tombo methylation data and (optional) label name.
	--deepsignal <[label:]fileC>  DeepSignal methylation data and (optional) label name.
	--gene <gtf_file>  Gene information in Gene Transfer Format.
	-o --output=<output>  Save dataframe to tsv.
	-p --plot=<plot>  Save plot(s) to pdf.
	--title=<title>  Title of the plot.


' -> doc

suppressMessages(library(docopt))
suppressMessages(library(tidyverse))

source("load_mod_data.R")
source("plot_mod_data.R")

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
    ## not yet implemented
    verbose("Not yet implemented\n")
    quit(save='no')
  }
  
  
  ## output merged dataframe
  if (!is.null(arguments$output)) {
    outfile <- arguments$output
    verbose("Writing to %s..\n", outfile)
    write_output(df, outfile)
  }
  
  ## plot and save to pdf
  if (!is.null(arguments$plot)) {
    outfile <- arguments$plot
    verbose("Plotting to %s..\n", outfile)
    plot_region(df, raw_regions, outfile)
  }
  
}

main(arguments)