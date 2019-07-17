suppressMessages(require(tidyverse))

load_file.nanopolish <- function(filename) {
  cols <- c("chromosome", "start", "read_name", "strand", "log_lik_ratio", 
            "num_motifs", "sequence")
  coltypes <- list(chromosome = 'c', start = 'i', read_name = 'c', 
                   log_lik_ratio = 'd', num_motifs = 'i', sequence = 'c',
                   strand = col_factor(levels=c('+', '-')))
  
  vroom(filename, col_select=cols, col_types=coltypes) %>%
    dplyr::rename(seqname = chromosome,
                  pos = start,
                  read_id = read_name)
}

preprocess.nanopolish <- function(df, motif="CG") {
  ## Assign same log-lik-ratio to all bases covered in motif
  offset <- unlist(gregexpr(pattern=motif, df$sequence)) - 1
  data <- df %>%
    uncount(num_motifs) %>%
    mutate(pos = pos - 5 + offset,
           log_lik_ratio = -1 * log_lik_ratio) %>%
    select(-sequence)
  
  ## Upon inspecting, a single read sometimes produces more than one statistic at a position
  ## Assuming statistic in the middle contribute the most, use weights from binomial distribution
  data <- data %>%
    group_by(seqname, pos, read_id, strand) %>%
    summarise(log_lik_ratio = sum(
      dbinom(0:(n()-1), n()-1, 0.5) * log_lik_ratio
    )) %>%
    ungroup() %>%
    mutate(pos = pos + 1,
           prob_meth = 1 / (1 + exp(log_lik_ratio)))
    
  data
}