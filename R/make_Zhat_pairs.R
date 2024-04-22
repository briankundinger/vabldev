#' @export
#'

make_Zhat_pairs <- function(Z_hat){
  data.frame(id_1 = Z_hat,
             id_2 = 1:n2) %>%
    filter(id_1 > 0)
}
