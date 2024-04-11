#' @export
#'

jaro <- function(complete_df, threshold = .5){
  pairs <- data.frame(a = complete_df$id_1,
                      b = complete_df$id_2,
                      gamma = NA,
                      r = complete_df$weight,
                      r_star = complete_df$fs_prob) %>%
    filter(a != 0)

  # command1: application of a preliminary filter to the input data
  filtered=pairs[pairs[,5]>threshold,]
  fs_matches <- filtered

  # command2: input preprocessing
  # counting of unique identifiers of records
  nA= length(unique(filtered[,1]))
  nB= length(unique(filtered[,2]))
  A=cbind(a=unique(filtered[,1]),A=1:nA)
  B=cbind(b=unique(filtered[,2]),B=1:nB)
  filtered =merge(B, filtered)
  filtered =merge(A, filtered)
  dat=t(filtered)

  # command3: preparing constraint matrix
  constr=array(rep(0,(nA+nB)*ncol(dat)), dim=c(nA+nB,ncol(dat)))
  p=rbind(matrix(rep(dat[2,],nA),c(nA,ncol(constr)),byrow=TRUE),
          matrix(rep(as.numeric(dat[4,])+nA,nB),c(nB,ncol(constr)),byrow=TRUE))
  constr [as.numeric(p)==row(constr)]=1

  # command4: preparing other LP parameters
  diseq=rep('<=',nA+nB)
  ones=rep(1,nA+nB)

  # target function
  coeff=dat[6,]

  # command5: LP execution
  ret= lpSolve::lp("max", coeff, constr, diseq, ones)
  # preparing the reduced set of pairs
  reduc <- t(dat[,ret$solution>0.9])

  Z_hat <- rep(0, n2)
  Z_hat[reduc[, 3]] <- reduc[, 1]

  return(list(Z_hat = Z_hat,
              fs_no_jaro = fs_matches))

  }
