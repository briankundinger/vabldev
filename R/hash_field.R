hash_field <- function(L_f, k, Lf_vec){
  level_seq <- seq_len(L_f)
  as.numeric(level_seq > 0) * 2 ^ ((level_seq) + (as.numeric(k > 1)  * Lf_vec[k]))
}
