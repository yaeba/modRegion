suppressMessages(require(tidyverse))


load_file.deepsignal <- function(filename) {
  colnames <- c("seqname", "pos", "strand", "x1", 
                "read_id", "x2", "x3", "prob_mod", "x4", "x5")
  cols <- c("seqname", "pos", "read_id", "strand", "prob_mod")
  coltypes <- list(seqname = 'c', pos = 'i', read_id = 'c', prob_mod = 'd',
                   strand = col_factor(levels=c('+', '-')))
  
  vroom(filename, delim='\t',
        col_select=cols, col_types=coltypes, col_names=colnames)
}


preprocess.deepsignal <- function(df, ...) {
  data <- df %>%
    mutate(pos = pos + 1,
           log_lik_ratio = log((1 - prob_mod) / prob_mod))
  
  data
}