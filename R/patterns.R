# technical pattern definition
#'@importFrom Rcpp evalCpp
#'@useDynLib techchart

#'@export
pattern.db <- function(type="all"){
  patterns = list()

  ############################# Head and Shoulder #########################
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
  plot <- list(color=adjustcolor("red",alpha.f = 0.4),lwd=5)
  patterns$HS <- list(call=call.expr, length=length, start=start,
                      threshold=threshold, plot=plot)

  ######################### Inverse Head and Shoulder ####################
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
  plot <- list(color=adjustcolor("blue",alpha.f = 0.4),lwd=5)
  patterns$IHS <- list(call=call.expr, length=length, start=start,
                       threshold=threshold, plot=plot)

  ######################### Broadening Top ####################
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!tchart.trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     R[1] < R[3] & R[3] < R[5] & R[2] > R[4]")
  length <- 5
  start <- 1
  threshold <- quote((R[2]+R[4])/2)
  plot <- list(color=adjustcolor("green",alpha.f = 0.4),lwd=5)
  patterns$BTOP <- list(call=call.expr, length=length, start=start,
                        threshold=threshold, plot=plot)

  ######################### Broadening Bottom ####################
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!tchart.trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     R[1] > R[3] & R[3] > R[5] & R[2] < R[4]")
  length <- 5
  start <- -1
  threshold <- quote((R[2]+R[4])/2)
  plot <- list(color=adjustcolor("brown",alpha.f = 0.4),lwd=5)
  patterns$BBOT <- list(call=call.expr, length=length, start=start,
                        threshold=threshold, plot=plot)

  ################## done, now check and return required patterns ###########
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

#'@export
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
    stop("not enough important points, time series too short")
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
  tpattern$vol <- vol
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
        n <- n+1; matches[[n]] <- list(data=patterns.out,points=points.out)
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
