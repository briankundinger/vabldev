#' @export
#'
fabl_mm <- function(hash, m_prior = 1, u_prior = 1,
                               alpha = 1, beta = 1, S = 1000, burn = round(S * .1),
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

  nonmatch_df <- data.frame(id_2 = 1:n2,
                            pattern = NA,
                            count = 1)

  weight_df <- lapply(1:n2, function(x){
    data.frame(id_2 = x,
               pattern = factor(1:P),
               count = hash_count_list[[x]])
  }) %>%
    do.call(rbind, .) %>%
    mutate(rn = row_number())
  #
  # %>%
  #   bind_rows(nonmatch_df) %>%
  #   mutate(rn = row_number())

  candidates_P <- 1:(P+1)
  #Z_compact <- vector("list", S)
  Z_list <- vector("list", S)

  #Z.SAMPS <- matrix(NA, nrow = n2, ncol = S)
  m_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  u_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  n_possible_list <- list()
  PI.SAMPS <- list()
  Z.temp <- rep(0, n1*n2)
  Z <- rep(n1+1, n2)
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

    # counts_list <- hash$pattern_hash_count_listord
    unique_weights <- exp(rowSums(ratio * unique_patterns, na.rm = TRUE))
    # counts_long <- hash$pattern_hash_count_listord %>%
    #   unlist
    # weights_long <- exp(rowSums(ratio * unique_patterns, na.rm = TRUE)) %>%
    #   rep(., n2)

    weight_df$count <- hash$hash_count_list %>%
      unlist()

    k  = 1
    #n_sampled_vec <- c()
    pi_vec <- c()
    n_possible_vec <- n2
    #n_possible <- n2
    matchable <- 1:n2
    removed_set <- c()


    while(TRUE){
      if(s == 1){
        n_prior = 0
      } else if(length(n_possible_list[[s - 1]]) < k){
        n_prior = 0
      } else {
        n_prior <- n_possible_list[[s - 1]][k]
      }

      n_possible <- n_possible_vec[k]
      n_prior <- min(n_possible, n_prior)

      pi <- rbeta(1, n_prior + alpha, n_possible - n_prior + beta)
      offset <- (n1 - k + 1) * (1 - pi) / pi

      weight_split <- weight_df %>%
        group_split(id_2)

      sampled_rows <- sapply(matchable, function(x){
        sample(c(0, weight_split[[x]]$rn),
               1,
               prob = c(1- pi, weight_split[[x]]$count * unique_weights * pi / (n1 - k)))
      }) %>%
        unname()

      removed <- sampled_rows[sampled_rows > 0]

      # sampled <- weight_df %>%
      #   filter(id_2 %in% matchable) %>%
      #   group_by(id_2) %>%
      #   sample_n(1, weight = weights) %>%
      #   ungroup() %>%
      #   tidyr::complete(id_2 = 1:n2)

      matchable <- weight_df %>%
        filter(rn %in% removed) %>%
        select(id_2) %>%
        pull()

      # removed <- weight_df %>%
      #   filter(rn %in% matches) %>%
      #   select(pattern) %>%
      #   pull() %>%
      #   as.integer

      weight_df$count[removed] <- weight_df$count[removed] - 1


      removed_set <- c(removed_set, removed)
      n_possible_vec <- c(n_possible_vec, length(matchable))

      #n_sampled_vec <- c(n_sampled_vec, length(matchable))
      pi_vec <- c(pi_vec, pi)

      k = k + 1

      if(length(matchable) == 0){
        break
      }

    }

    n_possible_list[[s]] <- n_possible_vec[-1]
    PI.SAMPS[[s]] <- pi_vec

    Z_list[[s]] <- weight_df[removed_set, ] %>%
      select(pattern, id_2) %>%
      arrange(id_2, pattern)

    # Z_compact[[s]] <- weight_df %>%
    #   filter(rn %in% removed_set) %>%
    #   tidyr::complete(id_2 = 1:n2) %>%
    #   select(pattern, id_2) %>%
    #   nest_by(id_2, .key = "pattern") %>%
    #   ungroup() %>%
    #   select(pattern)

    matches <- weight_df$pattern[removed_set] %>%
      data.frame(pattern = .) %>%
      group_by(pattern, .drop = F) %>%
      count() %>%
      pull()

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

  final_gibbs <- lapply(Z_list, function(y){
    temp <- y %>%
      group_by(id_2, pattern) %>%
      count() %>%
      ungroup() %>%
      group_split(id_2)

    rec_1 <- lapply(temp, function(w){
      w %>%
        group_split(pattern) %>%
        lapply(., function(z){
          sample_with_1(x = hash_to_file_1[[z$id_2]][[z$pattern]],
                        size = z$n)
        }) %>%
        unlist()
    }) %>%
      unlist()

    data.frame(id_1 = rec_1,
               id_2 = y$id_2)

    # %>%
    #     tidyr::complete(id_2 = 1:n2) %>%
    #     select(rec_1, id_2) %>%
    #     nest_by(id_2, .key = "id_1") %>%
    #     ungroup() %>%
    #     select(id_1)
  })
  # %>%
  #   do.call(cbind, .)



  #final_gibbs <- final_gibbs[ ,-(1:burn)]

  Z.SAMPS <- lapply((burn + 1):S, function(x){
    final_gibbs[[x]]
  })
  m_samps <- m_samps[ ,-(1:burn)]
  u_samps <- u_samps[ ,-(1:burn)]

  # Format PI.SAMPS

  list(Z = Z.SAMPS,
       #Z_list = Z_list,
       m = m_samps,
       u = u_samps)

}
