#' @export
#'

make_Zhat_pairs <- function(Z_hat){
  data.frame(id_1 = Z_hat,
             id_2 = 1:(length(Z_hat))) %>%
    filter(id_1 > 0)
}
