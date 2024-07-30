#' @export
#'
fabl_mm_partition <- function(hash, m_prior = 1, u_prior = 1,
                              alpha = 1, beta = 1, S = 1000, burn = round(S * .1),
                              show_progress = T, max_K = Inf, rejection_sampler = F,
                              rejection_iter = 10){
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
  highest_K <- 1
  Z <- matrix(0, nrow = n2, ncol = 1)
  Z_pattern <- matrix(0, nrow = n2, ncol = 1)

  Z_inv <- matrix(0, nrow = n1, ncol = 1)
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
  n_possible <- n2
  n_last_iter <- 0
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

    #matchable <- 1:n2
    #Z_pattern <- matrix(0, nrow = n2, ncol = highest_K)

    eta <- sapply(1:highest_K, function(k){
      rbeta(1,  n_last_iter[k] + alpha, n_possible[k] - n_last_iter[k] + beta)
    })

    for(j in 1:n2){
      k = 1
      #iter <- 0
      flag <- 1
      if(any(Z_pattern[j, ] > 0)){
        n_last_iter[Z_pattern[j, ] > 0] <- n_last_iter[Z_pattern[j, ] > 0] - 1
      }
      Z_inv[Z_inv == j] <- 0
      Z[j, ] <- NA
      Z[j, 1] <- 0
      Z_pattern[j, ] <- 0


      while(flag == 1){
        hash_weights <- hash_count_list[[j]] * unique_weights
        probs <- c(1 - eta[k], (eta[k] / (n1 - k + 1)) *  hash_weights)

        Z_pattern[j, k] <- sample(candidates, 1, prob = probs)
        if(Z_pattern[j, k] > 0){
          hash_count_list[[j]][Z_pattern[j, k]] <- hash_count_list[[j]][Z_pattern[j, k]] - 1
          index <- ceiling(runif(1) * length(hash_to_file_1[[j]][[Z_pattern[j, k]]]))
          record <- hash_to_file_1[[j]][[Z_pattern[j, k]]][index]
          hash_to_file_1[[j]][[Z_pattern[j, k]]] <- hash_to_file_1[[j]][[Z_pattern[j, k]]][hash_to_file_1[[j]][[Z_pattern[j, k]]] != record]
          #Z_samps[j, s, k] <- record
          Z[j, k] <- record
          Z_inv[record, which(Z_inv[record, ] == 0)[1]] <- j
          n_last_iter[k] <- n_last_iter[k] + 1
          n_possible <- c(n1, n_last_iter[-highest_K])
        }

        matching_complete <- Z_pattern[j, k] == 0

          if(rejection_sampler == F & matching_complete){
            break
          }

          if(matching_complete & rejection_sampler == T & k > 2){
            clean_set <- Z[j, ][!is.na(Z[j, ])]

            inverse <- Z_inv[clean_set, 1:(k-1)]

            clusters <- inverse %>%
              apply(., 1, sort, simplify = F) %>%
              lapply(., function(x){
                x[x != 0]
              }) %>%
              unique()
            if(length(clusters) == 1){
              break
            }

            n_last_iter[1:k] <- n_last_iter[1:k] - 1
            n_possible <- c(n1, n_last_iter[-highest_K])
            k <- 0
            hash_count_list <- hash$hash_count_list
            hash_to_file_1 <- hash$hash_to_file_1
            Z[j, ] <- NA
            Z[j, 1] <- 0
            Z_inv[Z_inv == j] <- 0
            iter <- iter + 1

            if(iter == rejection_iter){
              k = 1
              # Z[j, ] <- 0
              # Z_inv[Z_inv == j] <- 0
              while(TRUE){
                available <- Z_inv %>%
                  rowSums() == 0
                records_available <- which(available)
                temp_weights <- c(1 - eta[k], (eta[k] / sum(available)) * unique_weights[hash$pair_to_pattern[[j]][available]])
                record <- sample(c(0, records_available), 1, prob = temp_weights)
                Z[j, k] <- record
                if(record == 0){
                  break
                }
                Z_inv[record, which(Z_inv[record, ] == 0)[1]] <- j
                k <- k + 1
                if(k > max_K){
                  break
                }
              }
            }
          }





        k <- k + 1
        if(k > max_K){
          break
        }

        if(k > highest_K){
          highest_K <- highest_K + 1
          n_possible[k] <- n1
          n_last_iter[k] <- 0
          eta[k] <- rbeta(1,  n_last_iter[k] + alpha, n_possible[k] - n_last_iter[k] + beta)
          Z_pattern <- cbind(Z_pattern, 0)
          Z_inv <- cbind(Z_inv, 0)
          Z <- cbind(Z, NA)
          Z_samps <- abind::abind(Z_samps, array(NA, dim = c(n2, S, 1)), along = 3)
        }
      }
    }

    #colSums(Z > 0)
    #n_possible <- c(n1, colSums(Z_pattern > 0)[-highest_K]


    #n_possible_list[[s]] <- n_possible_vec[-1]
    eta_samps[[s]] <- eta

    matches <- factor(Z_pattern, levels = 0:P) %>%
      table() %>%
      .[-1] %>%
      unname()

    Z_samps[, s, ] <- Z
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
  pi_samps <- pi_samps[-(1:burn)]

  list(Z = Z_samps,
       m = m_samps,
       u = u_samps,
       pi = pi_samps)

}
