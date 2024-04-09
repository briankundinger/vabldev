sei <- function(x, R){
  if(length(x) <= R){
    return(x)
  } else {
    return(sample(x, R, replace = FALSE))
  }
}
