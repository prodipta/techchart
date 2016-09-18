
#' Find the top level trends in a time series
#' @param x xts, dataframe, matrix or vector, representing a time series
#' @param penalty starting penalty value, should be aBIG number
#' @return change points data (object of class cpttrends)
#' @seealso \code{\link[techchart]{cpt.trend}}
#' @export
find.major.trends <- function(x,penalty=20){
  n <- penalty + 1
  cpts <- list(); class(cpts) <- "cpttrend"
  for(i in 1:n){
    cpts <- cpt.trend(x,50,20,n-i)
    if(NROW(cpts$cpts > 0)) break
  }
  return(cpts)
}

#' Scale a dataframe to a unit square
#' @param x input dataframe with columns x and y
#' @export
range01 <- function(x){
  zeroshift <- min(x)
  x <- x - min(x)
  x <- x/max(x)
  return(x)
}
