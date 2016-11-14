# pattern search and find
#'@importFrom Rcpp evalCpp
#'@importFrom grDevices adjustcolor
#'@useDynLib techchart

#' @title Database of pattern definitions
#' @description definition of standard technical patterns
#' @param type string, any of "HS", "IHS", "BTOP", "BBOT" or "all"
#' @return a list if pattern definition
#' @export
pattern.db <- function(type="all"){
  patterns = list()

  #Head and Shoulder
  name = "Head and shoulder"
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!tchart.trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     (R[1] > R[2]) & (R[5] > R[4]) & (R[3] > R[1]) & (R[3] > R[5]) &
                     abs(R[1] - (R[1]+R[5])/2) < tolerance * (R[1]+R[5])/2 &
                     abs(R[5] - (R[1]+R[5])/2) < tolerance * (R[1]+R[5])/2 &
                     abs(R[2] - (R[2]+R[4])/2) < tolerance * (R[2]+R[4])/2 &
                     abs(R[4] - (R[2]+R[4])/2) < tolerance * (R[2]+R[4])/2
                     ")
  length <- 5
  start <- 1
  threshold <- quote((R[2]+R[4])/2)
  aspect.1 <- quote((R[5]-R[4])/(R[1]-R[2]))
  aspect.2 <- quote((R[5]-R[4])/(R[3]-R[2]))
  plot <- list(color=adjustcolor("red",alpha.f = 0.4),lwd=5)
  patterns$HS <- list(name=name, call=call.expr, length=length, start=start,
                      threshold=threshold, plot=plot, aspect.1=aspect.1, aspect.2=aspect.2)

  #Inverse Head and Shoulder
  name = "Inverse head and shoulder"
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!tchart.trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     (R[1] < R[2]) & (R[5] < R[4]) & (R[3] < R[1]) & (R[3] < R[5]) &
                     abs(R[1] - (R[1]+R[5])/2) < tolerance * (R[1]+R[5])/2 &
                     abs(R[5] - (R[1]+R[5])/2) < tolerance * (R[1]+R[5])/2 &
                     abs(R[2] - (R[2]+R[4])/2) < tolerance * (R[2]+R[4])/2 &
                     abs(R[4] - (R[2]+R[4])/2) < tolerance * (R[2]+R[4])/2
                     ")
  length <- 5
  start <- -1
  threshold <- quote((R[2]+R[4])/2)
  aspect.1 <- quote((R[5]-R[4])/(R[1]-R[2]))
  aspect.2 <- quote((R[5]-R[4])/(R[3]-R[2]))
  plot <- list(color=adjustcolor("blue",alpha.f = 0.4),lwd=5)
  patterns$IHS <- list(name=name, call=call.expr, length=length, start=start,
                       threshold=threshold, plot=plot, aspect.1=aspect.1, aspect.2=aspect.2)

  #Broadening Top
  name = "Broadening top"
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!tchart.trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     R[1] < R[3] & R[3] < R[5] & R[2] > R[4]")
  length <- 5
  start <- 1
  threshold <- quote((R[2]+R[4])/2)
  aspect.1 <- quote((R[5]-R[4])/(R[1]-R[2]))
  aspect.2 <- quote((R[5]-R[4])/(R[3]-R[2]))
  plot <- list(color=adjustcolor("orange",alpha.f = 0.4),lwd=5)
  patterns$BTOP <- list(name=name, call=call.expr, length=length, start=start,
                        threshold=threshold, plot=plot, aspect.1=aspect.1, aspect.2=aspect.2)

  #Broadening Bottom
  name = "Broadening bottom"
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!tchart.trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     R[1] > R[3] & R[3] > R[5] & R[2] < R[4]")
  length <- 5
  start <- -1
  threshold <- quote((R[2]+R[4])/2)
  aspect.1 <- quote((R[5]-R[4])/(R[1]-R[2]))
  aspect.2 <- quote((R[5]-R[4])/(R[3]-R[2]))
  plot <- list(color=adjustcolor("brown",alpha.f = 0.4),lwd=5)
  patterns$BBOT <- list(name=name, call=call.expr, length=length, start=start,
                        threshold=threshold, plot=plot, aspect.1=aspect.1, aspect.2=aspect.2)

  #done, now check and return required patterns
  ret.patterns <- list()
  index <-1
  if(type[1] != "all"){
    for(i in 1:NROW(patterns)){
      if(names(patterns)[i] %in% type){
        ret.patterns[[index]] <- patterns[[i]]
        index <- index + 1
      }
    }
  } else {
    ret.patterns <- patterns
  }

  return(ret.patterns)
}

#' Find matching patterns in a time series
#' @param x xts, a time series in xts format
#' @param pattern a list returned from a call to pattern.db()
#' @param tolerance tolerance used for pattern matching definition (vol multiple)
#' @param pip.tolerance tolerance used for searching important points (vol multiple)
#' @param find.all boolean, whether to return all matches or the latest one
#' @param trend.adjusted boolean, whether to adjust for trends during matching
#' @return a list of pattern matches (object of type tpattern)
#' @seealso \code{\link[techchart]{pattern.db}}
#' @seealso \code{\link[techchart]{find.pattern}}
#' @export
find.tpattern <- function(x, pattern=pattern.db("HS")[[1]], tolerance=0.25,
                          pip.tolerance=2, find.all= TRUE, trend.adjusted=FALSE){
  if(!xts::is.xts(x)){
    stop("require an xts time series object")
  }
  if(NROW(pattern)>1 && NROW(lengths(pattern))==1){
    warning("more than one patterns specified, first will be used")
    pattern <- pattern[[1]]
  }

  if(quantmod::is.OHLC(x)){
    ret <- na.omit(TTR::ROC(quantmod::Cl(x)))
  } else{
    ret <- na.omit(TTR::ROC(x[,1]))
  }
  vol <- sd(ret)
  tolerance <- vol*tolerance

  idx <- zoo::index(x)
  imppts <- find.imppoints(x,pip.tolerance)
  if(NROW(imppts$results)<pattern$length){
    stop("not enough points, time series too short")
  }
  z <- imppts$results

  k <- min(which(z$sign == pattern$start))
  if(!find.all){
    pattern.end <- ((-1)^(pattern$length-1))*pattern$start
    k <- max(which(z$sign==pattern.end)) - (pattern$length - 1)
  }
  l <- NROW(z)-pattern$length+1
  kl <- seq(k,l,2)

  pattern.env <- new.env(parent=globalenv())

  tpattern <- list(); tpattern$imppts <- imppts;
  tpattern$vol <- vol; class(tpattern) <- "tpattern"
  matches <- list(); n <- 0

  # loop through the important points
  for(k in kl){
    #assign optimas to pattern environment
    for(j in 1:pattern$length){
      assign(paste("E",j,sep=""),as.numeric(z$value[k+j-1]),envir=pattern.env)
      assign(paste("T",j,sep=""),as.numeric(z$pos[k+j-1]),envir=pattern.env)
    }
    # assign control variables
    curren.point <- ifelse(quantmod::is.OHLC(x),quantmod::Cl(x)[NROW(x)],x[NROW(x),1])
    assign("tolerance",tolerance,envir=pattern.env)
    assign("tchart.trend.adjusted",trend.adjusted,envir=pattern.env)
    assign("current.point", as.numeric(curren.point),envir = pattern.env)

    tryCatch({
      if(eval(pattern$call,pattern.env)){
        x.in.idx <- idx[z$pos[k:(k+pattern$length-1)]]
        x.in <- z$pos[k:(k+pattern$length-1)]
        y.in <- z$value[k:(k+pattern$length-1)]
        x.out <- z$pos[k]:z$pos[k+pattern$length-1]
        x.out.idx <- idx[z$pos[k]:z$pos[k+pattern$length-1]]
        y.out <- approx(x.in,y.in,x.out, n=NROW(x.out))
        patterns.out <- xts::as.xts(y.out$y, x.out.idx)
        colnames(patterns.out) <- c("pattern")
        points.out <- xts::as.xts(y.in, x.in.idx)
        colnames(points.out) <- c("points")
        xdate <- zoo::index(patterns.out)[NROW(patterns.out)]
        threshold <- eval(pattern$threshold,pattern.env)
        aspect.1 <- eval(pattern$aspect.1,pattern.env)
        aspect.2 <- eval(pattern$aspect.2,pattern.env)
        name <- pattern$name
        move <- annualize(patterns.out)
        duration <- as.double(zoo::index(patterns.out)[NROW(patterns.out)]
                              - zoo::index(patterns.out)[1])
        type= "complete"
        n <- n+1; matches[[n]] <- list(data=patterns.out,points=points.out,
                                       name=name, move=move, duration=duration,
                                       threshold=threshold, aspect.1=aspect.1,
                                       aspect.2=aspect.2, date=xdate, type=type)
      }
    }, error=function(e){
      #cat(paste("loop:",k,"\n"))
      #print(e)
      #print(zoo::index(x)[NROW(x)])
      #print(tail(z))
      #print(c(pattern.env$E1,pattern.env$E2,pattern.env$E3,pattern.env$E4))
    }
    )
    #k <- min(which(z$sign[(k+1):NROW(z)]==pattern$start)) + k
  }

  tpattern$matches <- matches
  return(tpattern)
}

#'@export
print.tpattern <- function(x,...){
  n <- NROW(x$matches)
  for(i in 1:n){
    y <- x$matches[[i]]
    cat(paste("------pattern matched on:",y$date,"--------")); cat("\n")
    cat(paste("name:",y$name)); cat("\n")
    cat(paste("type:",y$type)); cat("\n")
    cat(paste("move:",round(100*y$move,2)),"(percentage annualized)"); cat("\n")
    cat(paste("threshold:",round(y$threshold,2))); cat("\n")
    cat(paste("duration:",round(y$duration,2)),"(days)"); cat("\n")
  }
}

#'@export
summary.tpattern <- function(object, ...){
  x <- object
  return(x)
}

#'@export
print.summary.tpattern <- function(x,...){
  n <- NROW(x$matches)
  for(i in 1:n){
    y <- x$matches[[i]]
    cat(paste("------pattern matched on:",y$date,"--------")); cat("\n")
    cat(paste("name:",y$name)); cat("\n")
    cat(paste("type:",y$type)); cat("\n")
    cat(paste("move:",round(100*y$move,2)),"(percentage annualized)"); cat("\n")
    cat(paste("threshold:",round(y$threshold,2))); cat("\n")
    cat(paste("duration:",round(y$duration,2)),"(days)"); cat("\n")
  }
}
