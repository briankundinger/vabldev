#'
#' Internal Functions
#'
#'
#'
fs_to_ohe <- function(gamma, levels){
  unname(unlist(mapply(function(z, y){
    gamma_f <- rep(0, y)
    if(z == 0){
      return(gamma_f)
    }
    gamma_f[z] <- 1
    gamma_f
  }, z = gamma, y = levels)))
}
