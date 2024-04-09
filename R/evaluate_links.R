#' @export
#'
evaluate_links <- function(Z_hat, Z_true, n1){
  # Z_hat = Bayes Estimate of linkage structure (from BRL)
  # Z_true = True linkage structure
  # n1 = size of larger file

  n_links <- sum(Z_hat <= n1 & Z_hat > 0)
  n_matches <- sum(Z_true <= n1 & Z_true > 0)
  n_correct_links <- sum(Z_hat[Z_hat<=n1 & Z_hat > 0] ==
                           Z_true[Z_hat<=n1 & Z_hat > 0])
  recall <- n_correct_links/n_matches
  precision <- n_correct_links/n_links
  fmeasure <- 2 * (recall * precision) / (recall + precision)
  eval <- c(recall, precision, fmeasure)
  names(eval) <- c("Recall", "Precision", "Fmeasure")
  eval
}
