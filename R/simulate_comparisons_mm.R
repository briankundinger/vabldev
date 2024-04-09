#' Simulate comparison Vectors
#'
#' @export
#'

simulate_comparisons_mm <- function(m, u, levels, n1, n2, overlap){
  field_marker <- unlist(lapply(1:length(levels), function(x){
    rep(x, levels[x])
  }))
  FF <- length(levels)
  N <- n1 * n2
  ids <- expand.grid(1:n1, 1:n2)
  indicators <- matrix(NA, nrow = N, ncol = length(levels))



  df2_matches <- lapply(overlap, seq_len)
  matches_in_2 <- seq_len(overlap[1])
  match_breaks <- cumsum(c(0, overlap[ -length(overlap)]))
  df1_matches <- purrr::map2(df2_matches, match_breaks, ~.x + .y)

  Z_temp <- data.frame(target = unlist(df1_matches),
                       base = unlist(df2_matches))

  match_index <- apply(Z_temp, 1, function(x){
    which(ids[, 1] == unlist(x)[1] & ids[, 2] == unlist(x)[2])
  })

  Z_true <- Z_temp %>%
    setNames(c("id_1", "id_2")) %>%
    tidyr::complete(id_2 = 1:n2) %>%
    relocate(id_2, .after = id_1)


  # %>%
  #   nest_by(id_2)


  m_list <- split(m, field_marker)
  u_list <- split(u, field_marker)

  gamma_match <- sapply(m_list, function(x){
    sample(seq_along(x) - 1, sum(overlap), replace = T, x) + 1
  })

  gamma_nonmatch <- sapply(u_list, function(x){
    sample(seq_along(x) - 1, N - sum(overlap), replace = T, x) + 1
  })

  if(sum(overlap) == 0){
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


  list(comparisons = gamma,
       n1 = n1,
       n2 = n2,
       nDisagLevs = levels,
       Z_true = Z_true,
       Z_true_long = Z_temp)
}
