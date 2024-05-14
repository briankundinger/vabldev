#' @export
#'

variational_fastlink <- function(hash, threshold = 1e-6, tmax = 200, fixed_iterations = NULL, a_init = NULL,
                           b_init = TRUE, check_every = 10, store_every = check_every){

  ohe <- hash$ohe
  P <- dim(ohe)[1]
  total_counts <- hash$total_counts #N_p
  hash_count_list <- hash$hash_count_list
  field_marker <- hash$field_marker
  n1 <- as.double(hash$n1)
  n2 <- as.double(hash$n2)

  # Priors
  alpha <- rep(1, length(field_marker))
  beta <- rep(1, length(field_marker))
  alpha_lambda <- 1
  beta_lambda <- 1

  # Initialize
  if(is.null(a_init)){
    a <- rep(1, length(field_marker))
  } else {
    a <- a_init
  }
  if(b_init == T){
    b <- hash$ohe %>%
      sweep(., 1, hash$total_counts, "*") %>%
      colSums() + beta
  } else {
    b = rep(1, length(field_marker))
  }
    a_lambda <- 1
    b_lambda <- n1*n2 - 1 #initialize towards nonmatching, faster convergence
    #b_lambda <- 1

  t <- 1
  ratio <- 1
  elbo_seq <- vector()
  a_lambda_vec <- rep(NA, tmax)

  while(t <= tmax){
    a_sum <- a %>%
      split(., field_marker) %>%
      sapply(., sum) %>%
      digamma(.) %>%
      .[field_marker]

    a_chunk <- digamma(a) - a_sum

    b_sum <- b %>%
      split(., field_marker) %>%
      sapply(., sum) %>%
      digamma(.) %>%
      .[field_marker]
    b_chunk <- digamma(b) - b_sum

    m_p <- ohe %>%
      sweep(., 2, a_chunk, "*") %>%
      rowSums()

    u_p <- ohe %>%
      sweep(., 2, b_chunk, "*") %>%
      rowSums()

    # w_p
    weights = m_p - u_p

    # phi_single
    phi_tilde <- exp(digamma(a_lambda)  + weights)
    phi <- exp(log(phi_tilde) - log(phi_tilde + exp(digamma(b_lambda))))

    matching_weight_by_pattern <- hash$total_counts * phi
    total_match <- sum(matching_weight_by_pattern)

    AZ <- ohe %>%
      sweep(., 1, matching_weight_by_pattern, "*") %>%
      colSums()

    BZ <- ohe %>%
      sweep(., 1, total_counts - (matching_weight_by_pattern), "*") %>%
      colSums()

    a <- alpha + AZ
    b <- beta+ BZ

    a_lambda <- alpha_lambda + total_match
    b_lambda <- beta_lambda + (n1 * n2) - total_match

    a_lambda_vec[t] <- a_lambda


    t <- t + 1
    if(t > tmax){
      #print("Max iterations have passed before convergence")
      break
    }

    # if(!is.null(fixed_iterations)){
    #   if(t == fixed_iterations){
    #     break
    #   }
    # }
#print(t)
  }

  list(fs_probs = phi,
       pattern_weights = phi_tilde,
       a = a,
       b = b,
       a_lambda = a_lambda,
       b_lambda = b_lambda,
       elbo_seq = elbo_seq,
       t = t,
       a_lambda_vec = a_lambda_vec)

}

