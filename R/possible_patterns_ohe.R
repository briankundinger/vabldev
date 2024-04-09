possible_patterns_ohe <- function(levels){
  possible_values <- lapply(levels, function(x){
    c(0, seq_len(x))
  })
  possible_patterns <- data.frame(do.call(expand.grid, possible_values))

  data.frame(t(apply(possible_patterns, 1, function(x){
    fs_to_ohe(x, levels)
  })))
}
