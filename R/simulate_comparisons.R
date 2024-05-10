#' Simulate comparison Vectors
#'
#' @export

simulate_comparisons <- function(m, u, levels, n1, n2, overlap,
                                 previous_matches = 0){
  field_marker <- unlist(lapply(1:length(levels), function(x){
    rep(x, levels[x])
  }))
  FF <- length(levels)

  N <- n1 * n2
  ids <- expand.grid(1:n1, 1:n2)
  indicators <- matrix(NA, nrow = N, ncol = length(levels))

  df2matches <- seq_len(overlap)
  df1matches <- df2matches + previous_matches

  Z_true <- rep(0, n2)
  Z_true[df2matches] <- df1matches

  match_index <- which(ids[,1]  == (ids[,2]+ previous_matches))[seq_len(overlap)]

  m_list <- split(m, field_marker)
  u_list <- split(u, field_marker)

  gamma_match <- sapply(m_list, function(x){
    sample(seq_along(x) - 1, overlap, replace = T, x) + 1
  })

  gamma_nonmatch <- sapply(u_list, function(x){
    sample(seq_along(x) - 1, N - overlap, replace = T, x) + 1
  })

  if(overlap == 0){
    indicators <- gamma_nonmatch
  } else {
  indicators[match_index,] <- gamma_match
  indicators[-match_index,] <- gamma_nonmatch
  }

  ohe <- vector("list", FF)
  for(f in 1:FF){
    ohe[[f]] <- matrix(0, nrow = n1 * n2, ncol = levels[f])
    for(ell in 1:levels[f]){
      ohe[[f]][indicators[, f] == ell, ell] <- 1
    }
  }
  gamma <- do.call(cbind, ohe)


  # ohe <- purrr::map2(data.frame(indicators), levels, ~fs_to_ohe(.x, .y)) %>%
  #   do.call(cbind, .)

  list(comparisons = gamma,
       n1 = n1,
       n2 = n2,
       nDisagLevs = levels,
       Z_true = Z_true)
}
