#' @export
#'
estimate_links_fl <- function(out, hash, l_FNM=1, l_FM1=1, l_FM2=2, l_R=Inf,
                              nonmatch_label = "zero", resolve = T, flags_only = F){
  # "out" can be the output from either fabl or vabl.

  n1 <- hash$n1
  n2 <- hash$n2
  fs_probs <- out$fs_probs
  pattern_weights <- out$pattern_weights

  if(flags_only == T){
    flag_probs <- possible_records <- lapply(1:n2, function(j){
      record <- hash$flags[[j]]$eligible_records
      prob <- fs_probs[hash$flags[[j]]$eligible_patterns]
      data.frame(id_1 = record,
                 id_2 = j,
                 prob = prob)
    }) %>%
      do.call(rbind, .) %>%
      filter(prob > .5)

    return(flag_probs)
  }

  ids <- expand.grid(1:n1, 1:n2)
  all_probs <- data.frame(id_1 = ids[, 1],
                   id_2 = ids[, 2],
                   pattern = as.integer(hash$hash_id),
                   fs_prob = fs_probs[hash$hash_id],
                   weight = pattern_weights[hash$hash_id]) %>%
    group_by(id_2) %>%
    #mutate(vabl_est = weight / (sum(weight) + 1))
    mutate(vabl_est = weight / (sum(weight) + exp(digamma(out$b_lambda/n1)) + log(n1)))

  total_match_prob <- all_probs %>%
    group_by(id_2) %>%
    summarize(matching_prob = sum(vabl_est)) %>%
    pull()



  prob_no_link <- 1 - total_match_prob

  no_link_df <- data.frame(id_1 = 0, id_2 = 1:n2,
                           pattern = 0,
                           weight = 1,
                           vabl_est = prob_no_link)

  complete_df <- rbind(all_probs, no_link_df) %>%
    arrange(id_2, id_1)

  Zhat_df <- complete_df %>%
    group_by(id_2) %>%
    filter(vabl_est == max(vabl_est)) %>%
    filter(!duplicated(id_2)) # eliminates ties in vabl_est, ensures nrow == n2

  best_match <- Zhat_df$id_1
  prob_best_match <- Zhat_df$vabl_est
  link_indicator <- best_match > 0

  if(l_R == Inf){# if not using reject option and conditions of Theorem 1

    if (nonmatch_label == "n_1 + j"){
      Z_hat <- (n1+1):(n1+n2)
    }
    if (nonmatch_label == "zero"){
      Z_hat <- rep(0, n2)
    }
    threshold <- l_FM1/(l_FM1+l_FNM) +
      (l_FM2-l_FM1-l_FNM)*(1 - prob_no_link - prob_best_match)/(l_FM1+l_FNM)

    Z_hat[link_indicator & (prob_best_match > threshold)] <-
      best_match[link_indicator & (prob_best_match > threshold)]

  }

  if(l_R < Inf){
      Z_hat <- rep(-1, n2) # represents the reject option
      threshold <- 1 - l_R/l_FM1 +
        (l_FM2-l_FM1)*(1 - prob_no_link - prob_best_match) / l_FM1

      Z_hat[link_indicator & (prob_best_match > threshold) ] <-
        best_match[link_indicator & (prob_best_match > threshold) ]
      nonlink_indicator <- prob_no_link > 1 - l_R / l_FNM

      if (nonmatch_label == "n_1 + j"){
        Z_hat[nonlink_indicator] <- ((n1+1):(n1+n2))[nonlink_indicator]
      }
      if (nonmatch_label == "zero"){
        Z_hat[nonlink_indicator] <- 0
      }
  }

    # TODO: Write code for Theorem 2


  # Enforce one-to-one matching
  if (resolve == T){
    double_matches <- Z_hat[duplicated(Z_hat) & Z_hat > 0]
    if (l_R == Inf){
      to_resolve <- unlist(lapply(double_matches, function(x){
        df1_options <- which(Z_hat == x)
        df1_probs <- prob_best_match[df1_options]
        non_matches <- df1_options[-which.max(df1_probs)]
        non_matches
      }))
      Z_hat[to_resolve] <- 0
    } else {
      to_resolve <- unlist(lapply(double_matches, function(x){
        df1_options <- which(Z_hat == x)
        df1_options
      }))
      Z_hat[to_resolve] <- -1
    }
  }

  field_marker <- hash$field_marker

  m <- split(out$a, hash$field_marker) %>%
    lapply(., function(x){
      x/sum(x)
    }) %>%
    unlist()

  m_p <- sweep(hash$ohe, 2, log(m), "*") %>%
    rowSums() %>%
    exp()

  u <- split(out$b, hash$field_marker) %>%
    lapply(., function(x){
      x/sum(x)
    }) %>%
    unlist()

  u_p <- sweep(hash$ohe, 2, log(u), "*") %>%
    rowSums() %>%
    exp()

  w_p <- m_p / u_p

  E_m <- sum(w_p * m_p)
  KL <- sum(log(w_p) * m_p)

  params_fields <- list(m = m,
                        u = u)

  params_patterns <- list(m_p = m_p,
                          u_p = u_p,
                          w_p = w_p,
                          E_m = E_m,
                          KL = KL)

  fs <- jaro(complete_df)






  return(list(Z_hat = Z_hat,
              prob = prob_best_match,
              all_probs = complete_df,
              params_fields = params_fields,
              params_patterns = params_patterns,
              Z_hat_jaro = fs$Z_hat,
              fs_linkages = fs$fs_no_jaro))
}
