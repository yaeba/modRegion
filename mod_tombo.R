suppressMessages(require(tidyverse))

load_file.tombo <- function(filename) {
  cols <- c("chrm", "pos", "read_id", "strand", "stat")
  coltypes <- list(chrm = 'c', pos = 'i', read_id = 'c', stat = 'd',
                   strand = col_factor(levels=c('+', '-')))
  
  vroom(filename, col_select=cols, col_types=coltypes) %>%
    dplyr::rename(seqname = chrm,
                  log_lik_ratio = stat)
}

preprocess.tombo <- function(df, order, ...) {
  data <- df %>%
    mutate(pos = pos + 1,
           pos = ifelse(strand == "-", pos - 1, pos),
           prob_meth = 1 / (1 + exp(log_lik_ratio))) %>%
    select(order)
  
  data
}