# cpt trend functions
#'@importFrom Rcpp evalCpp
#'@useDynLib techchart

#' @export
print.cpttrend <- function(x,...){
  print(x$cpts, ...)
}

#' @export
summary.cpttrend <- function(object, ...){
  x <- object
  n <- NROW(x$cpts)+1
  seg.len <- rep(0,n)
  seg.return <- rep(0,n)
  seg.offset <- rep(0,n)
  for(i in 1:n){
    seg.len[i] <- NROW(x$segments[[i]])
    seg.return[i] <- annualize(x$segments[[i]])
    if(i>1){
      start <- as.numeric(x$segments[[i]][1])
      end <- as.numeric(x$segments[[i-1]][NROW(x$segments[[i-1]])])
      seg.offset[i] <- start/end -1
    }
  }
  y <- list(cpts=x$cpts, seg.len=seg.len, seg.return=seg.return,
            seg.offset=seg.offset)
  class(y) <- "summary.cpttrend"
  return(y)
}

#' @export
print.summary.cpttrend <- function(x,...){
  cat("change points:\n")
  print(x$cpts)
  cat("segments length summary:\n")
  print(summary(x$seg.len))
  cat("segments returns summary:\n")
  cat(x$seg.return); cat("\n")
  cat("segments offset summary:\n")
  cat(x$seg.offset)
}

#' @export
plot.cpttrend <- function(x,...){
  n <- NROW(x$cpts)+1
  if(xts::is.xts(x$data)){
    if(quantmod::is.OHLC(x$data)){
      plot(quantmod::Cl(x$data), main = "trend change points", xlab="x", ylab="y")
    }else{
      plot(x$data[,1], main = "trend change points", xlab="x", ylab="y")
    }
    for(i in 1:n){
      lines(xts::as.xts(x$segments[[i]]), ...)
    }
  } else{
    plot(x$data, type = "n", main="trend change points", xlab="x", ylab="y")
    lines(x$data)
    for(i in 1:n){
      lines(x$segments[[i]], ...)
    }
  }
  Sys.sleep(0)
}

#' change point analysis using binary segmentation
#' @param x xts, dataframe, matrix or vector, representing a time series
#' @param Q Maximum number of change points required
#' @param minseglen Minimum length of a trend segment
#' @param penalty Penalty value, increasing it reduces number of segments
#' @return change points data (object of class cpttrends)
#' @examples
#' x <- quantmod::getSymbols("^GSPC", auto.assign = FALSE)
#' x <- x["2015/"]
#' cpts <- cpt.trend(x,50,minseglen = 20, penalty = 5)
#' summary(cpts)
#' quantmod::chart_Series(x)
#' quantmod::add_TA(cpts$segments[[NROW(cpts$segments)]],on=1,lty=3)
#' @seealso \code{\link[changepoint]{cpt.mean}}
#' @export
cpt.trend <- function(x, Q=10, minseglen=10, penalty=1){

  if.xts <- F
  if(xts::is.xts(x)){
    if(quantmod::is.OHLC(x)){
      y <- data.frame(x=seq(1:NROW(x)),y=as.numeric(zoo::coredata(quantmod::Cl(x))))
    } else{
      y <- data.frame(x=seq(1:NROW(x)),y=as.numeric(zoo::coredata(x[,1])))
    }

    if.xts <- T
  } else if(NCOL(x)==1){
    y <- data.frame(x=seq(1:NROW(x)), y=as.numeric(x))
  } else{
    y <- data.frame(x=as.numeric(x[,1]), y=as.numeric(x[,2]))
  }

  z <- cpt_trend(y$x, y$y, Q, minseglen, penalty)
  z <- z[z!=0]

  if(NROW(z)<1){
    warning("No change point found, try relaxing parameters")
    cpts <- list()
    cpts$data <- x
    cpts$cpts <- NULL
    cpts$segments <- NULL
    class(cpts) <- "cpttrend"
    return(cpts)
  }
  if(NROW(z)==Q){
    warning("number of change points found same as requested,
            consider increasing the value of Q or penalty")
  }

  z <- sort(z)
  n <- NROW(z)+1
  segments <- list()

  for(i in 1:n){
    a <- ifelse(i==1,1,z[i-1]+1)
    b <- ifelse(i==n,NROW(y),z[i])
    m <- stats::lm(y$y[a:b]~y$x[a:b])
    fit <- predict(m)
    if(if.xts){
      segments[[i]] <- xts::as.xts(fit,zoo::index(x)[a:b])
    } else{
      segments[[i]] <- data.frame(x=y$x[a:b],y=fit)
    }
  }

  cpts <- list()
  cpts$data <- x
  cpts$cpts <- z
  cpts$segments <- segments
  class(cpts) <- "cpttrend"
  return(cpts)
}
