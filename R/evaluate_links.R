#' @export
#'
evaluate_links <- function(Z_hat, Z_true, n1, format = "Z"){
  # Z_hat = Bayes Estimate of linkage structure (from BRL)
  # Z_true = True linkage structure
  # n1 = size of larger file

  if(format == "Z"){

  n_links <- sum(Z_hat <= n1 & Z_hat > 0)
  n_matches <- sum(Z_true <= n1 & Z_true > 0)
  n_correct_links <- sum(Z_hat[Z_hat<=n1 & Z_hat > 0] ==
                           Z_true[Z_hat<=n1 & Z_hat > 0])
  recall <- n_correct_links/n_matches
  precision <- n_correct_links/n_links
  fmeasure <- 2 * (recall * precision) / (recall + precision)
  eval <- c(recall, precision, fmeasure)
  names(eval) <- c("Recall", "Precision", "Fmeasure")
  return(eval)
  }

  if(format == "pairs"){

    n_links <- dim(Z_hat)[1]
    n_matches <- dim(Z_true)[1]
    Z_hat_pair <- Z_hat %>%
      data.frame() %>%
      tidyr::unite("pair")

    Z_true_pair <- Z_true %>%
      data.frame() %>%
      tidyr::unite("pair")

    n_correct_links <- intersect(Z_hat_pair, Z_true_pair) %>%
      length()
    # n_correct_links <- apply(Z_hat, 1, function(x){
    #   any(Z_true[, 1] == x[1] & Z_true[, 2] == x[2])
    # }) %>%
    #   sum()

    recall <- n_correct_links/n_matches
    precision <- n_correct_links/n_links
    fmeasure <- 2 * (recall * precision) / (recall + precision)
    eval <- c(recall, precision, fmeasure)
    names(eval) <- c("Recall", "Precision", "Fmeasure")
    return(eval)
  }
}
