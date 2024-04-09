#' @export
#'
BRL_hash <- function(hash, m_prior = 1, u_prior = 1,
                          alpha = 1, beta = 1, S = 1000, burn = round(S * .1),
                          show_progress = T, seed = 0, reject_iter = round(hash$n2/length(hash$total_counts)),
                          mode = "rejection"){
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

  set.seed(seed)

  n1 <- hash$n1
  n2 <- hash$n2
  field_marker <- hash$field_marker

  unique_patterns <- hash$ohe
  pattern_counts <- hash$total_counts
  P <- nrow(unique_patterns)
  hash_count_list <- hash$hash_count_list
  hash_to_file_1 <-hash$hash_to_file_1
  pair_to_pattern <- hash$pair_to_pattern

  candidates_P <- 0:P
  Z_samps <- matrix(NA, nrow = n2, ncol = S)
  m_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  u_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  L_samps <- vector(length = S)

  # Initialization
  Z <- rep(0, n2)
  Z_pattern <- rep(0, n2)
  Z_inv <- rep(0, n1)
  L <- 0
  m <- u <- rep(0, length(field_marker))
  matches <- rep(0, P)

  # Gibbs
  for(s in 1:S){

    # Update m and u
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

    # Calculate weights
    ratio <- (log(m) - log(u)) %>%
      rep(., P) %>%
      matrix(., nrow = P, byrow = TRUE)

    unique_weights <- exp(rowSums(ratio * unique_patterns, na.rm = TRUE))

    for(j in 1:n2){

      if(Z[j] > 0){
        L <- L - 1
        Z_inv[Z[j]] <- 0
      }
      Z[j] <- 0

      empty_weight <- (n1 - L) * (n2 - L - 1 + beta) / (L + alpha)

      if(mode == "base"){
        available <- which(Z_inv == 0)
        temp_weights <- c(empty_weight, unique_weights[pair_to_pattern[[j]][available]])

        Z[j] <- sample(c(0, available), 1, prob = temp_weights)
        if(Z[j] > 0){
          Z_pattern[j] <- pair_to_pattern[[j]][Z[j]]
        }
        else{
          Z_pattern[j] <- 0
        }
      }
      else if(mode == "efficient"){
        n_current <- hash_count_list[[j]]
        for(k in 1:n2){
          if(Z[k] > 0){
            ind <- pair_to_pattern[[j]][Z[k]]
            n_current[ind] <- n_current[ind] - 1
          }
        }
        temp_weights <- n_current * unique_weights
        probs <- c(empty_weight, temp_weights)
        pattern <- sample(candidates_P, 1, prob = probs)
        if(pattern == 0){
          Z[j] <- 0
          Z_pattern[j] <- 0
        }
        else{
          flag_2 <- 1
          npj <- hash_count_list[[j]][pattern]
          while(flag_2 == 1){
            index <- ceiling(runif(1) * npj)
            i <- hash_to_file_1[[j]][[pattern]][index]
            if(Z_inv[i] == 0){
              Z[j] <- i
              Z_pattern[j] <- pattern
              flag_2 <- 0
            }
          }
        }
      }
      else if(mode == "rejection"){
        hash_weights <- hash_count_list[[j]] * unique_weights
        probs <- c(empty_weight, hash_weights)
        flag <- 1
        iter <- 0
        while(flag == 1){
          pattern <- sample(candidates_P, 1, prob = probs)
          if(pattern == 0){
            Z[j] <- 0
            Z_pattern[j] <- 0
            flag <- 0
          }
          else{
            index <- ceiling(runif(1) * hash_count_list[[j]][pattern])
            i <- hash_to_file_1[[j]][[pattern]][index]
            if(Z_inv[i] == 0){
              Z[j] <- i
              Z_pattern[j] <- pattern
              flag <- 0
            }
          }
          iter <- iter + 1
          if(iter == reject_iter){

            # O(n1)
            available <- which(Z_inv == 0)
            temp_weights <- c(empty_weight, unique_weights[pair_to_pattern[[j]][available]])

            Z[j] <- sample(c(0, available), 1, prob = temp_weights)
            if(Z[j] > 0){
              Z_pattern[j] <- pair_to_pattern[[j]][Z[j]]
            }
            else{
              Z_pattern[j] <- 0
            }
            flag <- 0
          }
        }
      }

      if(Z[j] > 0){
        L <- L + 1
        Z_inv[Z[j]] <- 1
      }
    }
    hash_matches <- factor(Z_pattern, levels = 0:P)
    df <- data.frame(hash_matches)
    matches <- df %>%
      group_by(hash_matches, .drop = F) %>%
      count() %>%
      filter(hash_matches != 0) %>%
      pull()

    Z_samps[,s] <- Z
    m_samps[,s] <- m
    u_samps[,s] <- u
    L_samps[s] <- L

    if(show_progress){
      if (s %% (S / 100) == 0) {
        flush.console()
        cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete",
                        sep = ""))
      }
    }
  }

  Z_samps[Z_samps == 0] <- n1 + 1

  list(Z = Z_samps,
       m = m_samps,
       u = u_samps,
       overlap = L_samps)

}
