#' SVI
#'
#' @export
#'
#'

svabl <- function(hash, threshold = 1e-6, tmax = 1000, fixed_iterations = NULL,
                          b_init = TRUE, B = min(1000, hash$n2),
                          k = 1, tau = 1, seed = 0,
                          check_every = 10, store_every = check_every){

  set.seed(seed)

  ohe <- hash$ohe
  P <- dim(ohe)[1]
  total_counts <- hash$total_counts #N_p
  hash_count_list <- hash$hash_count_list
  field_marker <- hash$field_marker
  n1 <- hash$n1
  n2 <- hash$n2

  # Priors
  alpha <- rep(1, length(field_marker))
  Beta <- rep(1, length(field_marker))
  alpha_pi <- 1
  beta_pi <- 1

  # Initialize
  a <- rep(1, length(field_marker))
  if(b_init == T){
    b <- hash$ohe %>%
      sweep(., 1, hash$total_counts, "*") %>%
      colSums() + Beta
  } else {
    b = rep(1, length(field_marker))
  }
  a_pi <- 1
  b_pi <- 1

  t <- 1
  ratio <- 1
  elbo_seq <- vector()
  adjustment <- n2 / B

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
    phi <- exp(digamma(a_pi) - digamma(n1) + weights)
    single <- exp(digamma(b_pi))

    # Phi_j
    batch <- sample(1:n2, B, replace = F)
    C <- sapply(batch, function(j){
      hash_count_list[[j]] %*% phi + single
    })

    # S(Phi, B)
    total_nonmatch <- adjustment * sum(single/ C)

    # N_p(B)
    total_counts <- lapply(batch, function(j){
      hash_count_list[[j]]
    }) %>%
      do.call(rbind, .) %>%
      colSums() * adjustment

    # N_p(Psi, B)
    K <- sapply(seq_along(batch), function(i){
      hash_count_list[[batch[i]]] / C[i]
    }) %>%
      rowSums() * adjustment

    epsilon <- (t + tau) ^ (-k)

    AZ <- ohe %>%
      sweep(., 1, phi * K, "*") %>%
      colSums()
    a <- (1 - epsilon) * a + epsilon * (alpha + AZ)

    BZ <- ohe %>%
      sweep(., 1, total_counts - (phi * K), "*") %>%
      colSums()

    b <- (1 - epsilon) * b + epsilon * (Beta + BZ)

    a_pi <- (1 - epsilon) * a_pi  +
      epsilon * (alpha_pi + n2 - total_nonmatch)
    b_pi <- (1 - epsilon) * b_pi +
      epsilon * (beta_pi + total_nonmatch)

    # ELBO
    if(t %% store_every == 0 | t == 1){

      C_full <- sapply(1:n2, function(j){
        hash_count_list[[j]] %*% phi + single
      })
      full_nonmatch <- sum(single/ C_full)

      elbo_pieces <- vector(length = 6)
      elbo_pieces[1] <- sapply(1:n2, function(j){
        sum(hash_count_list[[j]] *
              (phi *(weights - log(phi) + log(C_full[j]))/ C_full[j] + u_p))
      }) %>%
        sum(.)
      elbo_pieces[2] <-  single * sum(1/C_full *log(C_full)) + full_nonmatch * (log(n1) - log(single)) -log(n1)*n2
      elbo_pieces[3] <- lbeta(a_pi, b_pi) - lbeta(alpha_pi, beta_pi)
      elbo_pieces[4] <- sapply(list(a, b), function(y){
        split(y, field_marker) %>%
          sapply(., function(x){
            sum(lgamma(x)) - lgamma(sum(x))
          })%>%
          sum(.)
      }) %>%
        sum(.)
      elbo_pieces[5] <- -sapply(list(alpha, Beta), function(y){
        split(y, field_marker) %>%
          sapply(., function(x){
            sum(lgamma(x)) - lgamma(sum(x))
          })%>%
          sum(.)
      }) %>%
        sum(.)
      elbo_pieces[6] <- sum((alpha - a) * a_chunk + (Beta - b) * b_chunk)
      elbo <- sum(elbo_pieces)
      elbo_seq <- c(elbo_seq, elbo)
    }


    if(is.null(fixed_iterations)){
      if(t %% check_every == 0){
        t_elbo <- length(elbo_seq)
        ratio <- abs((elbo_seq[t_elbo] - elbo_seq[t_elbo - 1])/
                       elbo_seq[t_elbo - 1])
      }
      if(ratio < threshold){
        break
      }
    }

    t <- t + 1
    if(t > tmax){
      print("Max iterations have passed before convergence")
      break
    }

    if(!is.null(fixed_iterations)){
      if(t == fixed_iterations){
        break
      }
    }


  }

  C <- sapply(1:n2, function(j){
    hash_count_list[[j]] %*% phi + single
  })

  list(pattern_weights = phi,
       C = C,
       a = a,
       b = b,
       a_pi = a_pi,
       b_pi = b_pi,
       elbo_seq = elbo_seq,
       t = t)

}
