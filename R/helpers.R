
#' Find the top level trends in a time series
#' @param x xts, dataframe, matrix or vector, representing a time series
#' @param penalty starting penalty value, should be a BIG number (default 20)
#' @return change points data (object of class cpttrends)
#' @seealso \code{\link[techchart]{cpt.trend}}
#' @export
find.major.trends <- function(x,penalty=20){
  n <- penalty + 1
  cpts <- list(); class(cpts) <- "cpttrend"
  for(i in 1:n){
    tryCatch({
      cpts <- cpt.trend(x,50,20,n-i)
      if(NROW(cpts$cpts > 0)) break
    }, error=function(cond){
      message(cond)
      break
    }, warning=function(cod){
      # do nothing
    })

  }
  return(cpts)
}

