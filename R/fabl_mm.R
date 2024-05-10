#' @export
#'
fabl_mm <- function(hash, m_prior = 1, u_prior = 1,
                               alpha = 1, beta = 1, S = 1000, burn = round(S * .1),
                               show_progress = T, max_K = Inf, tau = 0){
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

  # nonmatch_df <- data.frame(id_2 = 1:n2,
  #                           pattern = NA,
  #                           count = 1)

  # weight_df <- lapply(1:n2, function(x){
  #   data.frame(id_2 = x,
  #              pattern = factor(1:P),
  #              count = hash_count_list[[x]])
  # }) %>%
  #   do.call(rbind, .) %>%
  #   mutate(rn = row_number())
  #
  # %>%
  #   bind_rows(nonmatch_df) %>%
  #   mutate(rn = row_number())

  candidates <- 0:P
  #Z_compact <- vector("list", S)
  if(max_K == Inf){
    Z_samps <- vector("list", S)
  } else {
    Z_samps <- array(NA, dim = c(n2, S, max_K))
    Z_samps[, , 1] <- 0
  }

  #Z.SAMPS <- matrix(NA, nrow = n2, ncol = S)
  m_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  u_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  n_possible_list <- list()
  L_list <- list()
  pi_samps <- list()
  #Z <- rep(n1+1, n2)
  #Z <- matrix(n1 + 1, nrow = n2, ncol = 1)
  L <- 0

  m <- u <- rep(0, length(field_marker))
  matches <- rep(0,P)
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
    k <-  1
    pi_vec <- c()
    n_possible_vec <- n2
    matchable <- 1:n2
    #Z_pattern <- vector()
    if(max_K == Inf){
      Z_samps[[s]] <- vector()
      Z_pattern <- matrix(0, nrow = n2, ncol = 1)
    } else {
      Z_pattern <- matrix(0, nrow = n2, ncol = max_K)
    }
    #removed_set <- c()


    while(TRUE){
      if(s == 1){
        n_last_iter = 0
      } else if(length(n_possible_list[[s - 1]]) < k){
        n_last_iter = 0
      } else {
        n_last_iter <- n_possible_list[[s - 1]][k]
      }

      # Z_pattern <- cbind(Z_pattern, rep(NA, n2))
      # Z_samps[[s]] <- cbind(Z_samps[[s]], rep(NA, n2))

      # if(k == 1){
      # Z_pattern <- cbind(Z_pattern, rep(0, n2))
      # if(max_K == Inf){
      # Z_samps[[s]] <- cbind(Z_samps[[s]], rep(0, n2))
      # }
      # } else {
      #   Z_pattern <- cbind(Z_pattern, rep(NA, n2))
      #   if(max_K == Inf){
      #     Z_samps[[s]] <- cbind(Z_samps[[s]], rep(NA, n2))
      #   }
      # }
      #
      if(max_K == Inf){
        if(k == 1){
          Z_samps[[s]] <- cbind(Z_samps[[s]], rep(0, n2))
        } else {
          Z_pattern <- cbind(Z_pattern, rep(NA, n2))
          Z_samps[[s]] <- cbind(Z_samps[[s]], rep(0, n2))
        }
      }

      n_possible <- n_possible_vec[k]
      n_last_iter <- min(n_possible, n_last_iter)

      pi <- rbeta(1, n_last_iter + alpha, n_possible - n_last_iter + k^tau)

      #pi <- rbeta(1, n_last_iter + alpha, n_possible - n_last_iter + beta_k)
      # pi <- rbeta(1, n_last_iter + alpha, n_possible - n_last_iter + k^3)
      #pi <- rbeta(1, n_last_iter + alpha, n2 - n_last_iter + beta)

      hash_weights <- lapply(hash_count_list, function(x){
        x * unique_weights
      })

      # weight_split <- weight_df %>%
      #   group_split(id_2)

      for(j in matchable){
        Z_pattern[j, k] <- sample(candidates, 1,
                       prob = c(1 - pi, hash_weights[[j]] * pi / (n1 - (k - 1))))
        if(Z_pattern[j, k] > 0){
          hash_count_list[[j]][Z_pattern[j, k]] <- hash_count_list[[j]][Z_pattern[j, k]] - 1
          index <- ceiling(runif(1) * length(hash_to_file_1[[j]][[Z_pattern[j, k]]]))
          record <- hash_to_file_1[[j]][[Z_pattern[j, k]]][index]
          hash_to_file_1[[j]][[Z_pattern[j, k]]] <- hash_to_file_1[[j]][[Z_pattern[j, k]]][hash_to_file_1[[j]][[Z_pattern[j, k]]] != record]
          #L[k] <- L + 1
          if(max_K == Inf){
            Z_samps[[s]][j, k] <- record
          } else {
            Z_samps[j, s, k] <- record
          }
        }
      }

      matchable <- seq_len(n2)[is.element(Z_pattern[, k] > 0, T)]
      n_possible_vec <- c(n_possible_vec, length(matchable))

      pi_vec <- c(pi_vec, pi)

      k <-  k + 1

      if(length(matchable) == 0){
        break
      }

      if(k > max_K){
        break
      }

    }

    n_possible_list[[s]] <- n_possible_vec[-1]
    pi_samps[[s]] <- pi_vec

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



  if(max_K == Inf){
    Z_samps <- Z_samps[-(1:burn)]
  } else {
    Z_samps <- Z_samps[, -(1:burn), ]
  }
  m_samps <- m_samps[ ,-(1:burn)]
  u_samps <- u_samps[ ,-(1:burn)]
  pi_samps <- pi_samps[-(1:burn)]

  # Format pi_samps

  list(Z = Z_samps,
       m = m_samps,
       u = u_samps,
       pi = pi_samps)

}
