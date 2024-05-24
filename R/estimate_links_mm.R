#' @export
#'


# estimate_links_mm <- function(out, hash, lFNM=1, lFM1=1, lFM2=2, lR=Inf,
#                               resolve = T){
#   #
#   # This is a complete copy of "linkrecords" from BRL, only modified
#   # so that it passes on the posterior link probabilities
#   #
#   #
#   #
#
#   # control the input
#   #if(!is.matrix(Z_samps)) stop("Z_samps should be a matrix")
#   #n2 <- nrow(Z_samps)
#   # make sure the labels in Z_samps are within the expected range
#   #if(max(Z_samps) > n1 + n2) stop("Labels in Z_samps exceed n1+n2")
#   # - positive losses
#   C0 <- (lFNM > 0) & (lFM1 > 0) & (lFM2 > 0) & (lR > 0)
#   # - conditions of Theorem 1 of Sadinle (2017)
#   C1 <- (lR == Inf) & (lFNM <= lFM1) & (lFNM + lFM1 <= lFM2)
#   # - conditions of Theorem 2 of Sadinle (2017)
#   C2 <- ((lFM2 >= lFM1) & (lFM1 >= 2*lR)) | ((lFM1 >= lFNM) & (lFM2 >= lFM1 + lFNM))
#   # - conditions of Theorem 3 of Sadinle (2017)
#   C3 <- (lFM2 >= lFM1) & (lFM1 >= 2*lR) & (lFNM >= 2*lR)
#   # check we can handle the specified losses
#   if(!C0) stop("Losses need to be positive")
#   if(!any(c(C1,C2,C3))) stop("Invalid configuration of losses")
#
#   n1 <- hash$n1
#   n2 <- hash$n2
#   Z_samps <- out$Z
#
#   #threshold <- lFM1/(lFM1+lFNM)
#   threshold <- 1/2
#   temp <- Z_samps %>%
#     do.call(rbind, .) %>%
#     group_by(id_1, id_2) %>%
#     count() %>%
#     mutate(prob = n / length(Z_samps)) %>%
#     filter(prob > threshold) %>%
#     ungroup()
#
#   Z_hat <- temp %>%
#     select(id_1, id_2)
#
#   prob <- temp %>%
#     select(prob)
#
#   double_matches <- Z_hat$id_1[duplicated(Z_hat$id_1)]
#
#   if(resolve == TRUE & length(double_matches) > 0){
#     if (lR == Inf){
#       to_resolve <- unlist(lapply(double_matches, function(x){
#         df2_options <- which(Z_hat$id_1 == x)
#         df2_probs <- prob[df2_options, ]
#         non_matches <- df2_options[-which.max(df2_probs)]
#         non_matches
#       }))
#       Z_hat <- Z_hat[-to_resolve, ]
#       prob <- prob[-to_resolve, ]
#     }
#   }
#
#   return(list(Z_hat = Z_hat,
#               prob = prob))
# }

estimate_links_mm <- function(out, hash, lFNM=1, lFM1=1, lFM2=2, lR=Inf,
                              resolve = T, transitivity = F){

  if(transitivity == T) {
    resolve <- F
  }

  threshold <- 1/2

  # control the input
  #if(!is.matrix(Z_samps)) stop("Z_samps should be a matrix")
  #n2 <- nrow(Z_samps)
  # make sure the labels in Z_samps are within the expected range
  #if(max(Z_samps) > n1 + n2) stop("Labels in Z_samps exceed n1+n2")
  # - positive losses
  C0 <- (lFNM > 0) & (lFM1 > 0) & (lFM2 > 0) & (lR > 0)
  # - conditions of Theorem 1 of Sadinle (2017)
  C1 <- (lR == Inf) & (lFNM <= lFM1) & (lFNM + lFM1 <= lFM2)
  # - conditions of Theorem 2 of Sadinle (2017)
  C2 <- ((lFM2 >= lFM1) & (lFM1 >= 2*lR)) | ((lFM1 >= lFNM) & (lFM2 >= lFM1 + lFNM))
  # - conditions of Theorem 3 of Sadinle (2017)
  C3 <- (lFM2 >= lFM1) & (lFM1 >= 2*lR) & (lFNM >= 2*lR)
  # check we can handle the specified losses
  if(!C0) stop("Losses need to be positive")
  if(!any(c(C1,C2,C3))) stop("Invalid configuration of losses")

  n1 <- hash$n1
  n2 <- hash$n2
  Z_samps <- out$Z
  #
  #   if(typeof(Z_samps) == "list"){
  #     samps <- length(Z_samps)
  #     all_draws <- do.call(cbind, Z_samps)
  #     probs <- apply(all_draws, 1, function(x){
  #       table(x)/samps
  #     }, simplify = F)
  #
  #     Z_hat <- purrr::imap(probs, ~data.frame(.x, base_id = .y)) %>%
  #       do.call(rbind, .) %>%
  #       rename(target_id = x,
  #              prob = Freq) %>%
  #       mutate(target_id = as.numeric(as.character(target_id))) %>%
  #       relocate(base_id, .before = prob) %>%
  #       filter(prob > threshold) %>%
  #       filter(target_id != 0)
  #
  #     probs_matches <- Z_hat$prob
  #     # Z_hat <- Z_hat %>%
  #     #   select(prob)
  #   } else {

  samps <- ncol(Z_samps)
  probs <- apply(Z_samps, 1, function(x){
    table(x)/samps
  }, simplify = F)

  probs_matches <- lapply(probs, function(x){
    matches <- x[which(x >.5)]
    matches <- matches[names(matches) != "0"]
    matches
  }) %>%
    unlist()

  target_id <- sapply(probs, function(x){
    matches <- names(which(x >.5))
    matches <- matches[matches != "0"]
    matches
  }, simplify = F)

  n_matches <- sapply(target_id, length)
  base_id <- sapply(1:n2, function(j){
    rep(j, n_matches[j])
  }) %>%
    unlist()

  target_id <- target_id %>%
    unlist() %>%
    as.numeric()

  Z_hat <- data.frame(target_id = target_id,
                      base_id = base_id,
                      prob = probs_matches)

  double_matches <- Z_hat$target_id[duplicated(Z_hat$target_id)]

  if(resolve == TRUE & length(double_matches) > 0){
    if (lR == Inf){
      to_resolve <- unlist(lapply(double_matches, function(x){
        index <- which(Z_hat$target_id == x)
        base_probs <- probs_matches[index]
        non_matches <- index[-which.max(base_probs)]
        non_matches
      }))
      Z_hat <- Z_hat[-to_resolve, ]
      probs_matches <- probs_matches[-to_resolve]
    }
  }


  identify_conflicts <- function(set, mms){
    common_entities <- sapply(seq_along(mms), function(j){
      intersect(set, mms[[j]]) %>%
        length()
    })
    max(0, which(common_entities > 0 & common_entities < length(set)))
  }

  if(transitivity == TRUE){

    mms_df <- Z_hat %>%
      group_split(base_id, .keep = T)
    # mms <- Z_hat %>%
    #   group_split(base_id) %>%
    #   lapply(., `[[`, "target_id")
    #
    # mms_prob <- Z_hat %>%
    #   group_split(base_id) %>%
    #   lapply(., `[[`, "prob")
    # unique_mms <- unique(mms)
    # unique_mms_map <- match(mms, unique_mms)
    # conflicts <- sapply(unique_mms, identify_conflicts, unique_mms)

    # mms_df[[201]] <- data.frame(target_id = 25,
    #                             base_id = 201,
    #                             prob = .9)

    mms <- mms_df %>%
      lapply(., `[[`, "target_id")

    mms_prob <- mms_df %>%
      lapply(., `[[`, "prob")


    unique_mms <- unique(mms)
    unique_mms_map <- match(mms, unique_mms)

    set_id_df <- lapply(seq_along(mms), function(j){
      data.frame(mms_df[[j]], set_id = unique_mms_map[j])
    }) %>%
      do.call(rbind, .)

    conflicts <- lapply(unique_mms, identify_conflicts, unique_mms)
    conflicts_df <- lapply(seq_along(conflicts), function(x){
      data.frame(set_1 = conflicts[[x]], set_2 = x)
    }) %>%
      do.call(rbind, .) %>%
      filter(set_1 >0)

    #if(any(conflicts >0)){
    if(nrow(conflicts_df > 0)){
      # conflicts_df <- data.frame(code_1 = conflicts) %>%
      #   mutate(code_2 = row_number()) %>%
      #   filter(code_1 > 0)

      for(i in seq_len(nrow(conflicts_df))){
        higher_prob <- set_id_df %>%
          filter(set_id %in% conflicts_df[i, ]) %>%
          group_by(base_id) %>%
          mutate(total_prob = sum(prob)) %>%
          ungroup() %>%
          filter(total_prob == max(total_prob)) %>%
          select(set_id) %>%
          pull() %>%
          unique()

        lower_prob <-  conflicts_df[i, ][conflicts_df[i, ] != higher_prob]

        set_id_df <- set_id_df %>%
          filter(!(set_id %in% lower_prob))
        # mms[[which(unique_mms_map ==  lower_prob)]] <- 0
        # mms_prob[[which(unique_mms_map ==  lower_prob)]] <- 0
      }
    }

    Z_hat <- set_id_df %>%
      select(target_id, base_id)

    probs_matches <- set_id_df %>%
      select(prob)
  }

  # Z_hat <- lapply(1:n2, function(j){
  #   data.frame(target_id = mms[[j]],
  #              base_id = j,
  #              prob = mms_prob[j])
  # }) %>%
  #   do.call(rbind, .) %>%
  #   filter(target_id != 0)



  return(list(Z_hat = Z_hat,
              prob = probs_matches))


}



