#' @export
#'
estimate_links<- function(out, hash, l_FNM=1, l_FM1=1, l_FM2=2, l_R=Inf,
                          nonmatch_label = "zero", resolve = T){
  # "out" can be the output from either fabl or vabl.

  n1 <- hash$n1
  n2 <- hash$n2
  params_fields <- NULL
  params_patterns <- NULL
  complete_df <- NULL


  if(names(out)[1] == "Z"){
    Z_samps <- out$Z

    #Z_samps[Z_samps > n1+1] <- n1+1

    samps <- ncol(Z_samps)
    probs <- apply(Z_samps, 1, function(x){
      table(x)/samps
    }, simplify = F)
    prob_no_link <- sapply(probs, function(x){
      1 - sum(x[names(x) != n1 + 1])
    })
    # prob_no_link <- sapply(probs, function(x){
    #   1 - sum(x[names(x) != 0])
    # })
    Z_hat <- rep(0, n2)
    best_match <- sapply(probs, function(x){
      names(which.max(x))
    }) %>%
      as.numeric()
    prob_best_match <- sapply(probs, function(x){
      max(x)
    })
    link_indicator <- best_match < n1 + 1
    all_probs <- NULL
  }

  if(names(out)[1] == "pattern_weights"){

    n2 <- hash$n2
    pattern_probs <- lapply(1:n2, function(j){
      out$pattern_weights/out$C[j]
    })

    possible_records <- lapply(1:n2, function(j){
      record <- c(hash$flags[[j]]$eligible_records, 0)
      prob <- c(pattern_probs[[j]][hash$flags[[j]]$eligible_patterns],
                exp(digamma(out$b_pi)) / out$C[j]) %>%
        unname()

      data.frame(record, prob)
    })

    max_prob <- lapply(possible_records, function(x){
      x[which.max(x$prob), ]
    }) %>%
      do.call(rbind, .)

    best_match <- max_prob$record
    prob_best_match <- max_prob$prob
    prob_no_link <- out$b_pi/out$C
    link_indicator <- best_match > 0

    ids <- expand.grid(1:n1, 1:n2)
    complete_df <- data.frame(id_1 = ids[, 1],
                     id_2 = ids[, 2],
                     pattern = as.integer(hash$hash_id),
                     weight = out$pattern_weights[hash$hash_id]) %>%
      group_by(id_2) %>%
      rbind(., data.frame(id_1 = 0,
                          id_2 = 1:n2,
                          pattern = 0,
                          weight = out$b_pi)) %>%
      arrange(id_2, id_1) %>%
      mutate(prob = weight / (sum(weight)))
  }

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

  if(names(out)[1] == "pattern_weights"){
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
  }

  return(list(Z_hat = Z_hat,
              prob = prob_best_match,
              all_probs = complete_df,
              params_fields = params_fields,
              params_patterns = params_patterns))
}
