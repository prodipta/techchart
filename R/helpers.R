annualize <- function(x){
  if(!xts::is.xts(x) | NROW(x)<2){
    stop("input is not in xts format or has less than 2 points")
  }
  period <- as.double(zoo::index(x)[NROW(x)] - zoo::index(x)[1])
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

#' Find the last best-fit enveloping lines of a given time series
#' @param x xts object representing a time series
#' @param tscore envelope score threshold
#' @return returns the trend channel object (of type class tchannel)
#' @seealso \code{\link[techchart]{find.tchannel}}
#' @export
find.trend.channel <- function(x, tscore=(0.25)^2){
  tolerance <- seq(1.1,1.5,0.05)
  error_threshold <- tscore
  score <- 999

  # set up the return object
  tchannel <- list()
  tchannel$xlines <- NA
  tchannel$name <- NA
  tchannel$dir <- NA
  tchannel$upperlimit <- NA
  tchannel$lowerlimit <- NA
  tchannel$duration <- NA
  tchannel$midlinemove <- NA
  tchannel$maxlinemove <-NA
  tchannel$minlinemove <- NA
  tchannel$aspectratio <- NA
  tchannel$score <- NA
  tchannel$strength <- NA
  tchannel$fit <- NA
  class(tchannel) <- "tchannel"
  na_tchannel <- tchannel

  for(tol in tolerance){
    tryCatch({
      new_tchannel <- find.tchannel(x,tol)
      if(!is.na(new_tchannel$name)){
        if(new_tchannel$score >= score){
          if(score < error_threshold)return(tchannel)
          next
        }
        tchannel <- new_tchannel
        score <- tchannel$score
      }
    }, error=function(cond){
      # do nothing
    }, warning=function(w){
      # do nothing
    })
  }

  if(!is.na(tchannel$score)){
    if(tchannel$score < (1.5*error_threshold))return(tchannel)
  }
  return(na_tchannel)
}



#' Find current and/ or forming technical patter
#' @param x xts, a time series in xts format
#' @param pattern a list returned from a call to pattern.db()
#' @param tolerance a list tolerance used for pattern matching definition (vol multiple)
#' @param pip.tolerance a list tolerance used for searching important points (vol multiple)
#' @return a list of pattern matches (object of type tpattern)
#' @seealso \code{\link[techchart]{pattern.db}}
#' @seealso \code{\link[techchart]{find.tpattern}}
#' @export
find.pattern <- function(x, pattern=pattern.db("HS")[[1]], tolerance=c(0.25,1),
                         pip.tolerance=c(1.1,1.5)){

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
        x1 <- as.numeric(quantmod::Cl(x)[NROW(x)])*(1+last.extrema*move)
        x1 <- data.frame(t(rep(x1,NCOL(x)))); colnames(x1) <- colnames(x)
        x1 <- xts::as.xts(x1,zoo::index(x)[NROW(x)]+1)
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
