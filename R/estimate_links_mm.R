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
                              resolve = T){
  #
  # This is a complete copy of "linkrecords" from BRL, only modified
  # so that it passes on the posterior link probabilities
  #
  #
  #

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


  samps <- length(Z_samps)
  all_draws <- do.call(cbind, Z_samps)
  probs <- apply(all_draws, 1, function(x){
    table(x)/samps
  }, simplify = F)

  # lapply(probs, function(x){
  #   names(which(x > threshold)) %>%
  #     as.numeric() %>%
  #     as.
  # })


  # Z_hat <- purrr::map(probs, ~names(which(.x > threshold))) %>%
  #   purrr::imap(., ~data.frame(target_id = .x, base_id = .y)) %>%
  #   do.call(rbind, .) %>%
  #   filter(target_id > 0)

  Z_hat <- purrr::imap(probs, ~data.frame(.x, base_id = .y)) %>%
    do.call(rbind, .) %>%
    rename(target_id = x,
           prob = Freq) %>%
    relocate(base_id, .before = prob) %>%
    #mutate(target_id = as.integer(target_id)) %>%
    filter(prob > threshold) %>%
    filter(target_id != 0)

  # prob_no_link <- sapply(probs, function(x){
  #   sum(x[names(x) == 0])
  # })
  # threshold <- 1/2
  # temp <- Z_samps %>%
  #   do.call(rbind, .) %>%
  #   group_by(id_1, id_2) %>%
  #   count() %>%
  #   mutate(prob = n / length(Z_samps)) %>%
  #   filter(prob > threshold) %>%
  #   ungroup()
  #
  # Z_hat <- temp %>%
  #   select(id_1, id_2)
  #
  # prob <- temp %>%
  #   select(prob)
  #
  # double_matches <- Z_hat$id_1[duplicated(Z_hat$id_1)]
  #
  # if(resolve == TRUE & length(double_matches) > 0){
  #   if (lR == Inf){
  #     to_resolve <- unlist(lapply(double_matches, function(x){
  #       df2_options <- which(Z_hat$id_1 == x)
  #       df2_probs <- prob[df2_options, ]
  #       non_matches <- df2_options[-which.max(df2_probs)]
  #       non_matches
  #     }))
  #     Z_hat <- Z_hat[-to_resolve, ]
  #     prob <- prob[-to_resolve, ]
  #   }
  # }

  # return(list(Z_hat = Z_hat,
  #             prob = prob))
return(list(Z_hat = Z_hat))
}
