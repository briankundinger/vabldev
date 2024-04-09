#' @export
#'

fabl <- function(hash, m_prior = 1, u_prior = 1, alpha = 1, beta = 1,
                 S = 1000, burn = round(S * .1),
                 show_progress = T){
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

  n1 <- hash$n1
  n2 <- hash$n2
  field_marker <- hash$field_marker

  unique_patterns <- hash$ohe
  pattern_counts <- hash$total_counts
  P <- nrow(unique_patterns)
  hash_count_list <- hash$hash_count_list
  hash_to_file_1 <-hash$hash_to_file_1

  candidates_P <- 0:P
  Z_samps <- matrix(0, nrow = n2, ncol = S)
  m_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  u_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  L_samps <- vector(length = S)
  pi_samps <- vector(length = S)

  # Initialize
  Z <- rep(0, n2)
  L <- 0
  m <- u <- rep(0, length(field_marker))
  matches <- rep(0,P)

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

    hash_weights <- lapply(hash_count_list, function(x){
      x * unique_weights
    })

    pi <- rbeta(1, L + alpha, n2 - L + beta)

    for(j in 1:n2){
      if(Z[j] > 0){
        L <- L - 1
      }
      Z[j] <- sample(candidates_P, 1,
                     prob = c(1 - pi, hash_weights[[j]] * pi / n1))
      if(Z[j] > 0){
        index <- ceiling(runif(1) * length(hash_to_file_1[[j]][[Z[j]]]))
        Z_samps[j, s] <- hash_to_file_1[[j]][[Z[j]]][index]
        L <- L + 1
      }
    }
    hash_matches <- factor(Z, levels = 0:P)
    df <- data.frame(hash_matches)
    matches <- df %>%
      group_by(hash_matches, .drop = F) %>%
      count() %>%
      filter(hash_matches != 0) %>%
      pull()

    m_samps[,s] <- m
    u_samps[,s] <- u
    L_samps[s] <- L
    pi_samps[s] <- pi

    if(show_progress){
      if (s %% (S / 100) == 0) {
        flush.console()
        cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete", sep = ""))
      }
    }
  }

  Z_samps[Z_samps == 0] <- n1 + 1

  list(Z = Z_samps[, -(1:burn)],
       m = m_samps[, -(1:burn)],
       u = u_samps[, -(1:burn)],
       overlap = L_samps[-(1:burn)],
       pi = pi_samps[-(1:burn)])

}
