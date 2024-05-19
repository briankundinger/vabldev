#' @export
#'
#'
#'
#'

get_mms <- function(out, hash, transitivity = T){

  identify_conflicts <- function(set, mms){
    common_entities <- sapply(seq_along(mms), function(j){
      intersect(set, mms[[j]]) %>%
        length()
    })
    max(0, which(common_entities > 0 & common_entities < length(set)))
  }

  n1 <- hash$n1
  n2 <- hash$n2
  Z_samps <- out$Z
  samps <- ncol(Z_samps)
  # if(typeof(Z_samps) == "list"){
  #   samps <- length(Z_samps)
  #   View(Z_samps)
  #
  # } else {
  #
  # samps <- ncol(Z_samps)
  # }

  # thing <- array(0, dim = c(2, 3, 4)
  # abind(thing, along = #)

  mms <- vector("list", n2)
  mms_probs <- vector("double", n2)

  for(j in 1:n2){
    samples <- Z_samps[j, , ] %>%
      apply(., 1, sort, simplify = F)

    unique_sets <- samples %>%
      unique()

    prob <- match(samples, unique_sets) %>%
      table(.)/ samps

    mms[[j]] <- unique_sets[[which.max(prob)]]
    mms_probs[j] <- max(prob)
  }



  if(transitivity == TRUE){
    unique_mms <- unique(mms)
    unique_mms_map <- match(mms, unique_mms)
    conflicts <- sapply(unique_mms, identify_conflicts, unique_mms)

    if(any(conflicts >0)){
      conflicts_df <- data.frame(code_1 = conflicts) %>%
        mutate(code_2 = row_number()) %>%
        filter(code_1 > 0)

      for(i in seq_len(nrow(conflicts_df))){
        higher_prob <- data.frame(set = unique_mms_map,
                                  probs = mms_probs) %>%
          mutate(record = row_number()) %>%
          filter(set %in% conflicts_df[i, ]) %>%
          mutate(size = sapply(mms[record], length)) %>%
          mutate(total_prob = probs * size) %>%
          filter(total_prob == max(total_prob)) %>%
          select(set) %>%
          pull()

        lower_prob <-  conflicts_df[i, ][conflicts_df[i, ] != higher_prob]
        mms[which(unique_mms_map ==  lower_prob)] <- 0
      }
    }
  }

  Z_hat <- lapply(1:n2, function(j){
    data.frame(target_id = mms[[j]],
               base_id = j,
               prob = mms_probs[j])
  }) %>%
    do.call(rbind, .) %>%
    filter(target_id != 0)

  return(list(Z_hat = Z_hat[, 1:2],
              prob = Z_hat$prob))
}
