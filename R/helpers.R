annualize <- function(x){
  if(!is.xts(x) | NROW(x)<2){
    stop("input is not in xts format or has less than 2 points")
  }
  period <- as.double(index(x)[NROW(x)] - index(x)[1])
  years <- period/365.25
  r <- as.numeric(x[NROW(x),1])/as.numeric(x[1,1])
  r <- (r)^(1/years) - 1
  return(r)
}


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


#' Find the past, current and forming technical patter
#' @param x xts, dataframe, matrix or vector, representing a time series
#' @return pattern match object (object of class tpattern)
#' @seealso \code{\link[techchart]{find.tpattern}}
#' @export
find.pattern <- function(x, pattern=pattern.db("HS")[[1]], tolerance=c(0.25,0.5),
                         pip.tolerance=c(1.5,1.25)){

  m <- NROW(pip.tolerance); n <- NROW(tolerance)

  for(i in 1:m){
    for(j in 1:n){
      tryCatch({
        # check if we have a pattern at the end
        tpattern <- find.tpattern(x,pattern,tolerance[j],pip.tolerance[i], find.all = FALSE)
        if(NROW(tpattern$matches)>0) return(tpattern)

        #if not, check if one is forming
        last.extrema <- as.numeric(tpattern$imppts$results$sign[NROW(tpattern$imppts$results)])
        move <- tpattern$vol*pip.tolerance[i]
        x1 <- as.numeric(Cl(x)[NROW(x)])*(1+last.extrema*move)
        x1 <- data.frame(t(rep(x1,NCOL(x)))); colnames(x1) <- colnames(x)
        x1 <- as.xts(x1,index(x)[NROW(x)]+1)
        y <- rbind(x,x1);
        tryCatch({
          tpattern <- find.tpattern(y,pattern,tolerance[j],pip.tolerance[i], find.all = FALSE)
          if(NROW(tpattern$matches)>0){
            tpattern$virtual <- TRUE
            return(tpattern)
          }
        }, error=function(e){}
        )
      }, error=function(e){}
      )
    }
  }

  # no pattern found
  return(NULL)
}
