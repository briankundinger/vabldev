#' @export
#'
swap_base_file <- function(cd){
  old_n1 <- cd$n1
  old_n2 <- cd$n2
  ohe <- cd$comparisons

  ids <- expand.grid(1:old_n1, 1:old_n2)
  new_gamma <- data.frame(ohe, id_1 = ids[, 1], id_2 = ids[, 2]) %>%
    arrange(id_1, id_2) %>%
    select(-id_1, -id_2) %>%
    as.matrix()
  Z_true_long <- NULL
  if("Z_true_long" %in% names(cd)){
    Z_true_long <- cd$Z_true_long
  names(Z_true_long) <- c("base", "target")
  }

  list(comparisons = new_gamma,
       n1 = old_n2,
       n2 = old_n1,
       nDisagLevs = cd$nDisagLevs,
       Z_true_long, Z_true_long)

  #TODO Add functionality to convert Z_true
}
