#' @export
#'
fabl_mm <- function(hash, m_prior = 1, u_prior = 1,
                               alpha = 1, beta = 1, S = 1000, burn = round(S * .1),
                               show_progress = T, max_K = Inf){
  # Implements bipartite record linkage with BK Sampling Mechanism
  #
  # Arguments
  # comparisons = list calculated from from BRL::compareRecords
  # m.prior = prior distribution for m parameters
  # u.prior= prior distribution for u parameters
  # alpha = first parameter of prior for linkage probability
  # beta = second parameter of prior for linkage probability
  # S = number of Gibbs iterations
  # burn = number of iterations to be discarded as burn-in
  # show_progress = set to false to show simulation progress

  n1 <- as.double(hash$n1)
  n2 <- as.double(hash$n2)
  field_marker <- hash$field_marker

  unique_patterns <- hash$ohe
  pattern_counts <- hash$total_counts
  P <- nrow(unique_patterns)
  hash_count_list <- hash$hash_count_list
  hash_to_file_1 <-hash$hash_to_file_1


  candidates <- 0:P
  #Z_compact <- vector("list", S)
  highest_K <- 1
  Z_samps <- array(NA, dim = c(n2, S, highest_K))
  Z_samps[, , 1] <- 0

  m_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  u_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  n_possible_list <- list()
  L_list <- list()
  eta_samps <- list()
  L <- 0

  m <- u <- rep(0, length(field_marker))
  matches <- rep(0,P)
  eta <- rbeta(1, alpha, beta)
  n_k <- c(n2, 0)
  #set.seed(1)

  # Gibbs
  for(s in 1:S){


    AZ <- sweep(unique_patterns, MARGIN = 1, STAT = matches, FUN = "*") %>%
      colSums() %>%
      unname()

    nonmatches <- pattern_counts - matches

    BZ <- sweep(unique_patterns, MARGIN = 1, STAT = nonmatches, FUN = "*") %>%
      colSums() %>%
      unname()


    m_post <- m_prior + AZ
    u_post <- u_prior + BZ

    m_post <- split(m_post, field_marker)
    m <- as.vector(unlist(sapply(m_post, function(x){
      prob <- MCMCpack::rdirichlet(1, x)
      prob/sum(prob)
    })))

    u_post <- split(u_post, field_marker)
    u <- as.vector(unlist(sapply(u_post, function(x){
      prob <- MCMCpack::rdirichlet(1, x)
      prob/sum(prob)
    })))


    ratio <- (log(m) - log(u)) %>%
      rep(., P) %>%
      matrix(., nrow = P, byrow = TRUE)

    unique_weights <- exp(rowSums(ratio * unique_patterns, na.rm = TRUE))
    hash_count_list <- hash$hash_count_list
    hash_to_file_1 <- hash$hash_to_file_1

    #pi_vec <- c()
    n_possible_vec <- n2
    matchable <- 1:n2
    Z_pattern <- matrix(0, nrow = n2, ncol = 1)

    for(k in seq_along(n_k)[-1]){
      eta[k-1] <- rbeta(1, n_k[k] + alpha, n_k[k-1] - n_k[k] + beta)
    }
    k <-  1

    while(TRUE){
      # if(s == 1){
      #   n_last_iter = 0
      # } else if(length(n_possible_list[[s - 1]]) < k){
      #   n_last_iter = 0
      # } else {
      #   n_last_iter <- sum(Z_samps[matchable, s - 1, k] > 0, na.rm = T)
      # }

      Z_pattern <- cbind(Z_pattern, rep(NA, n2))

      #n_possible <- n_possible_vec[k]

      if(k > length(eta)){
        eta_k <- rbeta(1, alpha, length(matchable) + beta)
      } else {
        eta_k <- eta[k]
      }


      hash_weights <- lapply(hash_count_list, function(x){
        x * unique_weights
      })

      for(j in matchable){
        Z_pattern[j, k] <- sample(candidates, 1,
                                  prob = c(1 - eta_k, hash_weights[[j]] * eta_k / (n1 - (k - 1))))
        if(Z_pattern[j, k] > 0){
          hash_count_list[[j]][Z_pattern[j, k]] <- hash_count_list[[j]][Z_pattern[j, k]] - 1
          index <- ceiling(runif(1) * length(hash_to_file_1[[j]][[Z_pattern[j, k]]]))
          record <- hash_to_file_1[[j]][[Z_pattern[j, k]]][index]
          hash_to_file_1[[j]][[Z_pattern[j, k]]] <- hash_to_file_1[[j]][[Z_pattern[j, k]]][hash_to_file_1[[j]][[Z_pattern[j, k]]] != record]
          Z_samps[j, s, k] <- record
        }
      }

      matchable <- seq_len(n2)[is.element(Z_pattern[, k] > 0, T)]
      n_possible_vec <- c(n_possible_vec, length(matchable))

      #pi_vec <- c(pi_vec, pi)

      k <-  k + 1

      if(k > highest_K){
        Z_samps <- abind::abind(Z_samps, array(NA, dim = c(n2, S, 1)), along = 3)
        highest_K <- highest_K + 1
      }

      if(length(matchable) == 0){
        break
      }

      if(k > max_K){
        break
      }
    }

    match_indicator <- Z_pattern > 0
    n_k <- c(n2, colSums(match_indicator, na.rm = T))

    #n_possible_list[[s]] <- n_possible_vec[-1]
    eta_samps[[s]] <- eta

    matches <- factor(Z_pattern, levels = 0:P) %>%
      table() %>%
      .[-1] %>%
      unname()

    #Z.SAMPS[,s] <- Z
    m_samps[,s] <- m
    u_samps[,s] <- u

    if(show_progress){
      if (s %% (S / 100) == 0) {
        flush.console()
        cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete", sep = ""))
      }
    }
  }


  Z_samps <- Z_samps[, -(1:burn), ]
  m_samps <- m_samps[ ,-(1:burn)]
  u_samps <- u_samps[ ,-(1:burn)]
  eta_samps <- eta_samps[-(1:burn)]

  list(Z = Z_samps,
       m = m_samps,
       u = u_samps,
       eta = eta_samps)
}
