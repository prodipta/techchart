# cpt trend functions
#'@importFrom Rcpp evalCpp
#'@useDynLib techchart

range01 <- function(x){
  zeroshift <- min(x)
  x <- x - min(x)
  x <- x/max(x)
  return(x)
}
cluster_grouping <-function(d,cluster,tolerance, FUN="max", is.rounded=TRUE){
  n <- NROW(d);
  all_idx <- (1:n)
  groups <- rep(0,n)
  for(i in 1:n){
    if(NROW(which(!is.na(all_idx)))<1)break
    cut <- cutree(cluster,i)
    for(j in 1:i){
      idx <- which(cut==j)
      idx <- idx[which(idx %in% all_idx)]
      x <- d[idx]
      if(NROW(x)<1)next
      if(NROW(x)==1){
        groups[idx] <- ifelse(is.rounded,round(x),x)
        all_idx[idx] <- NA
        next
      }
      if(max(x)-min(x)<tolerance){
        y <- do.call(FUN,list(x))
        groups[idx] <- ifelse(is.rounded,round(y),y)
        all_idx[idx] <- NA
      }
    }
  }
  return(groups)
}
merge_lines <- function(xlines, tolerance=2){

  #return if only one member
  if(NROW(xlines)<2)return(xlines)
  # cluster along start and end lines
  d <- unique(xlines$end);
  # if less than two, no clustering
  if(NROW(d)<2)return(xlines)

  # run clutering
  cluster <- hclust(dist(d))
  d1 <- cluster_grouping(d,cluster,2*tolerance, FUN="max",FALSE)
  endtag <- d1[match(xlines$end,d)]
  d <- unique(xlines$start); cluster <- hclust(dist(d))
  d1 <- cluster_grouping(d,cluster,2*tolerance, FUN="min",FALSE)
  starttag <- d1[match(xlines$start,d)]

  #merge same start and end points based on score
  xlines$key <- paste(starttag,"-",endtag)
  keys <- unique(xlines$key)
  lapply(1:NROW(keys), FUN=function(i){
    cut <- xlines[which(xlines$key == keys[i]),]
    if(NROW(cut)>1){
      idx <- which(cut$score==min(cut$score))
      if(NROW(idx)>1)idx <- which(cut$fit==max(cut$fit))
      x <- cut[idx,]
      xlines$r[which(xlines$key==keys[i])] <<- 0
      xlines <<- rbind(xlines,x)
    }
  })
  xlines$key <- NULL
  xlines <- xlines[xlines$r != 0,]

  return(xlines)
}
filter_lines <- function(xlines, imppts, flag=1, ptheta=0.25, pscore=0.3,
                         pfit=0.95, nsize=10, nlines=10){

  # nothing to do if no lines found
  if(NROW(xlines)<1)return(xlines)

  # compute the slope...
  #if(flag==1){
    #xlines$xtheta <- (180/pi)*atan((as.numeric(imppts$maxima$value[xlines$end]) - as.numeric(imppts$maxima$value[xlines$start]))
    #                      /(as.numeric(imppts$maxima$pos[xlines$end]) - as.numeric(imppts$maxima$pos[xlines$start])))
  #} else{
    #xlines$xtheta <- (180/pi)*atan((as.numeric(imppts$minima$value[xlines$end]) - as.numeric(imppts$minima$value[xlines$start]))
    #                      /(as.numeric(imppts$minima$pos[xlines$end]) - as.numeric(imppts$minima$pos[xlines$start])))
  #}
  xlines$xtheta <- xlines$theta - 90
  # ... and eliminate near-vertical lines
  if(!is.na(ptheta)){
    xlines <- xlines[abs(cos(xlines$xtheta*pi/180))>ptheta,]
  }

  # return if too few lines left
  if(NROW(xlines)<1)return(xlines)
  if(NROW(xlines)<nlines){
    xlines <- xlines[with(xlines, order(-end,score,-fit)),]
    return(xlines)
  }

  # check if envelop score is acceptable
  if(!is.na(pscore)){
    cutoffscore <- as.numeric(quantile(unique(xlines$score), pscore))
    tmplines <- xlines[xlines$score < pscore,]
    if(NROW(tmplines)<1)return(xlines)
    xlines <- tmplines
  }
  # return if too few lines left
  if(NROW(xlines)<1)return(xlines)
  if(NROW(xlines)<nlines){
    xlines <- xlines[with(xlines, order(-end,score,-fit)),]
    return(xlines)
  }

  # check if fit score is acceptable
  if(!is.na(pfit)){
    tmplines <- xlines[xlines$fit > pfit,]
    if(NROW(tmplines)<1)return(xlines)
    xlines <- tmplines
  }
  # return if too few lines left
  if(NROW(xlines)<1)return(xlines)
  if(NROW(xlines)<nlines){
    xlines <- xlines[with(xlines, order(-end,score,-fit)),]
    return(xlines)
  }

  # finally filter based on size
  if(!is.na(nsize)){
    s1 <- as.numeric(imppts$maxima$pos)[xlines$end]
    s2 <- as.numeric(imppts$maxima$pos)[xlines$start]
    if(flag==-1){
      s1 <- as.numeric(imppts$minima$pos)[xlines$end]
      s2 <- as.numeric(imppts$minima$pos)[xlines$start]
    }
    tmplines <- xlines
    tmplines$size <- s1 - s2
    tmplines <- tmplines[tmplines$size > nsize,]
    if(NROW(tmplines)<1)return(xlines)
    xlines <- tmplines
    xlines$size <- NULL
  }
  # return if too few lines left
  if(NROW(xlines)<1)return(xlines)
  if(NROW(xlines)<nlines){
    xlines <- xlines[with(xlines, order(-end,score,-fit)),]
    return(xlines)
  }

  # enough lines left, sort, cut and return
  xlines <- xlines[with(xlines, order(-end,score,-fit)),]
  xlines <- xlines[1:nlines,]

  return(xlines)
}
hough_lines <- function(x, flag, r.tol=0.02, theta.tol=1, s=2){

  rbucket <- seq(-1.42,1.42,r.tol); abucket <- seq(1,180,theta.tol)
  xlines <- houghtransform(range01(x$x),range01(x$y),flag,rbucket, abucket,s)
  xlines$r <- rbucket[xlines$r+1]; xlines$theta <- abucket[xlines$theta+1]
  xlines <- xlines[xlines$r != 0,]
  xlines$start <- xlines$start + 1
  xlines$end <- xlines$end + 1
  xlines <- xlines[which(xlines$end != xlines$start),]
  return(xlines)
}


#' Find enveloping lines of a given time series
#' @param x xts object representing a time series
#' @param tolerance tolerance specification for important points
#' @param nlines max number of lines to return
#' @param pscore envelope score, default value is NA (no filtering done)
#' @param pfit fit socre, default value is NA (no filtering done)
#' @param r.tol r (polar) coordinate tolerance in unit square
#' @param theta.tol theta coordinate tolerance in degrees
#' @param s minimum number less than the # of colinear points on the line
#' @param merge.tol tolerance for merging lines with similar start or end points
#' @seealso \code{\link[techchart]{find.tchannel}}
#' @export
find_lines <- function(x, tolerance=2, nlines=10, pfit=NA, pscore=NA,
                       r.tol=0.02, theta.tol=2, s=2, merge.tol=2){
  if.xts <- F
  if(xts::is.xts(x)){
    if.xts <- T
  }

  # find the important points and switch to 1 unit square map
  imppts <- find.imppoints(x,tolerance)
  xmin <- imppts$minima; xmax <- imppts$maxima
  xmin <- data.frame(x=as.numeric(zoo::coredata(xmin$pos)),y=as.numeric(zoo::coredata(xmin$value)))
  xmax <- data.frame(x=as.numeric(zoo::coredata(xmax$pos)),y=as.numeric(zoo::coredata(xmax$value)))

  # start with min points, find lines and merge and filter
  minlines <- hough_lines(xmin,-1,r.tol,theta.tol,s)
  if(NROW(minlines)>1){
    minlines <- merge_lines(minlines,merge.tol)
  }
  if(NROW(minlines)>1){
    minlines <- merge_lines(minlines,merge.tol)
  }
  if(NROW(minlines)>1){
    minlines <- filter_lines(minlines,imppts, -1,nlines=nlines, pfit = pfit, pscore = pscore)
  }

  # now same with max points
  maxlines <- hough_lines(xmax,1,r.tol,theta.tol,s)
  if(NROW(maxlines)>1){
    maxlines <- merge_lines(maxlines,merge.tol)
  }
  if(NROW(maxlines)>1){
    maxlines <- merge_lines(maxlines,merge.tol)
  }
  if(NROW(maxlines)>1){
    maxlines <- filter_lines(maxlines,imppts, 1,nlines=nlines, pfit = pfit, pscore = pscore)
  }

  # build the return object
  xlines <- list()
  xlines$imppts <- imppts
  xlines$maxlist <- maxlines
  xlines$minlist <- minlines

  # compute the lines for plot for minima points...
  minlist <- list()

  if(NROW(minlines)>0){
    for(i in 1:NROW(minlines)){
      r <- minlines$r[i]
      theta <- minlines$theta[i]
      start <- minlines$start[i]; end <- minlines$end[i]
      s2 <- as.numeric(imppts$minima$pos)[minlines$end[i]]
      s1 <- as.numeric(imppts$minima$pos)[minlines$start[i]]
      y2 <- as.numeric(quantmod::Lo(x)[s2]);y1 <- as.numeric(quantmod::Lo(x)[s1])
      idx1 <- zoo::index(x)[s1:s2]
      idx <- zoo::index(x)[s1:NROW(x)]
      z <- Hmisc::approxExtrap(c(1,NROW(idx1)),c(y1,y2),seq(1:NROW(idx)),n=NROW(idx))
      lastval <- as.numeric(z$y[NROW(idx)]); lastx <- as.numeric(quantmod::Cl(x)[NROW(x)])
      maxx <- max(as.numeric(quantmod::Cl(x))); minx <- min(as.numeric(quantmod::Cl(x)))
      b_out_range <- lastval > maxx | lastval < minx
      b_no_extrapolate <- lastval < 0.9*lastx | lastval > 1.1*lastx
      if(b_out_range & b_no_extrapolate){
        z <- approx(c(1,NROW(idx1)),c(y1,y2),seq(1:NROW(idx)),n=NROW(idx))
      }
      minlist[[i]] <- na.omit(xts::as.xts(z$y,idx))
    }
  }

  # and for maxima points...
  maxlist <- list()
  if(NROW(maxlines)>0){
    for(i in 1:NROW(maxlines)){
      r <- maxlines$r[i]
      theta <- maxlines$theta[i]
      start <- maxlines$start[i]; end <- maxlines$end[i]
      s2 <- as.numeric(imppts$maxima$pos)[maxlines$end[i]]
      s1 <- as.numeric(imppts$maxima$pos)[maxlines$start[i]]
      y2 <- as.numeric(quantmod::Hi(x)[s2]);y1 <- as.numeric(quantmod::Hi(x)[s1]);
      idx1 <- zoo::index(x)[s1:s2]
      idx <- zoo::index(x)[s1:NROW(x)]
      z <- Hmisc::approxExtrap(c(1,NROW(idx1)),c(y1,y2),seq(1:NROW(idx)),n=NROW(idx))
      lastval <- as.numeric(z$y[NROW(idx)]); lastx <- as.numeric(quantmod::Cl(x)[NROW(x)])
      maxx <- max(as.numeric(quantmod::Cl(x))); minx <- min(as.numeric(quantmod::Cl(x)))
      b_out_range <- lastval > maxx | lastval < minx
      b_no_extrapolate <- lastval < 0.9*lastx | lastval > 1.1*lastx
      if(b_out_range & b_no_extrapolate){
        z <- approx(c(1,NROW(idx1)),c(y1,y2),seq(1:NROW(idx)),n=NROW(idx))
      }
      maxlist[[i]] <- na.omit(xts::as.xts(z$y,idx))
    }
  }

  xlines$maxlines <- maxlist
  xlines$minlines <- minlist

  # we are done
  return(xlines)
}

find.lines <- function(x, tolerance=1.5, n=3, pscore=(0.05)^2, pfit=0.85,
                            force.one=!return.all, return.all=FALSE){

  last.idx <- zoo::index(x)[NROW(x)]

  # get the lines
  xlines <- find_lines(x,tolerance,n)
  if(return.all)return(xlines)

  # find the best fit max line
  if(NROW(xlines$maxlist)>0){
    nn <- NROW(xlines$maxlist)
    last.pos <- max(xlines$maxlist$end)
    scorecard <- data.frame(pos=rep(0,nn), score=rep(0,nn), valid=rep(0,nn),
                            mid=rep(0,nn))
    for(i in 1:nn){
      idx <- NROW(xlines$maxlines[[i]])
      b_ending <- zoo::index(xlines$maxlines[[i]])[idx] == last.idx
      b_pos <- b_ending & abs(xlines$maxlist$end[i] - max(xlines$maxlist$end)) < 3
      b_score <- xlines$maxlist$score[i] == min(xlines$maxlist$score) |
        xlines$maxlist$score[i] < pscore
      b_fit <- xlines$maxlist$fit[i] == min(xlines$maxlist$fit) |
        xlines$maxlist$fit[i] > pfit
      scorecard$pos[i] <- b_pos
      scorecard$score[i] <- b_score
      scorecard$valid[i] <- b_pos & b_score & b_fit
      scorecard$mid[i] <- !b_ending & b_score
    }
    idx <- which(scorecard$valid==TRUE)
    if(NROW(idx)< 1 & force.one){
      idx <- which(scorecard$pos==TRUE & scorecard$score==TRUE)
      if(NROW(idx)< 1){
        idx <- which(scorecard$pos==TRUE)
      }
    }
    if(NROW(idx)>1){
      idx <- which(xlines$maxlist$score==min(xlines$maxlist$score))
      if(NROW(idx)>1){
        idx <- which(xlines$maxlist[idx,]$strength
                     ==max(xlines$maxlist[idx,]$strength))
      }

      if(NROW(idx)>1)idx <- idx[1]
    }
    idx2 <- which(scorecard$mid==TRUE)
    idx <- c(idx,idx2)
    if(NROW(idx)>0){
      xlines$maxlist <- xlines$maxlist[idx,]
      for(i in 1:nn){
        if(!(i%in%idx))xlines$maxlines[[i]] <- NULL
      }
    } else{
      xlines$maxlist <- xlines$maxlist[0,]
      xlines$maxlines <- list()
    }
  }

  # now repeat for minima
  if(NROW(xlines$minlist)>0){
    nn <- NROW(xlines$minlist)
    last.pos <- max(xlines$minlist$end)
    scorecard <- data.frame(pos=rep(0,nn), score=rep(0,nn), valid=rep(0,nn),
                            mid=rep(0,nn))
    for(i in 1:nn){
      idx <- NROW(xlines$minlines[[i]])
      b_ending <- zoo::index(xlines$minlines[[i]])[idx] == last.idx
      b_pos <- b_ending & abs(xlines$minlist$end[i] - max(xlines$minlist$end)) < 2
      b_score <- xlines$minlist$score[i] == min(xlines$minlist$score) |
        xlines$minlist$score[i] < pscore
      b_fit <- xlines$minlist$fit[i] == min(xlines$minlist$fit) |
        xlines$minlist$fit[i] > pfit
      scorecard$pos[i] <- b_pos
      scorecard$score[i] <- b_score
      scorecard$valid[i] <- b_pos & b_score & b_fit
      scorecard$mid[i] <- !b_ending & b_score
    }
    idx <- which(scorecard$valid==TRUE)
    if(NROW(idx)<1 & force.one){
      idx <- which(scorecard$pos==TRUE & scorecard$score==TRUE)
      if(NROW(idx)< 1){
        idx <- which(scorecard$pos==TRUE)
      }
    }
    if(NROW(idx)>1){
      idx <- which(xlines$minlist$score==min(xlines$minlist$score))
      if(NROW(idx)>1){
        idx <- which(xlines$minlist[idx,]$strength==
                       max(xlines$minlist[idx,]$strength))
      }

      if(NROW(idx)>1)idx <- idx[1]
    }
    idx2 <- which(scorecard$mid==TRUE)
    idx <- c(idx,idx2)
    if(NROW(idx)>0){
      xlines$minlist <- xlines$minlist[idx,]
      for(i in 1:nn){
        if(!(i%in%idx))xlines$minlines[[i]] <- NULL
      }
    } else{
      xlines$minlist <- xlines$minlist[0,]
      xlines$minlines <- list()
    }
  }

  #done, return the lines
  return(xlines)

}

#' Find most current best-fit enveloping lines of a given time series
#' @param x xts object representing a time series
#' @param tolerance tolerance specification for important points, expressed as times the standard deviation
#' @param n number of lines to to choose the best ones from
#' @param pscore envelope score, defines how best the lines envelopes the time series
#' @param pfit fit socre, defines how much the lines deviates from the extreme points it connects
#' @return returns the trend channel object (of type class tchannel)
#' @examples
#' x <- quantmod::getSymbols("^GSPC", auto.assign = FALSE)
#' x <- x["2015/"]
#' quantmod::chart_Series(x)
#' tchannel <- find.tchannel(x)
#' tchannel
#' quantmod::add_TA(tchannel$xlines$maxlines[[1]],on=1)
#' quantmod::add_TA(tchannel$xlines$minlines[[1]],on=1)
#' @seealso \code{\link[techchart]{find_lines}}
#' @export
find.tchannel <- function(x, tolerance=1.5, n=3, pscore=(0.05)^2,
                               pfit=0.85){

  try(xlines <- find.lines(x, tolerance, pscore, pfit, force.one = TRUE))
  tchannel <- list()
  tchannel$xlines <- xlines
  tchannel$name <- NA
  tchannel$type <- NA
  tchannel$dir <- NA
  tchannel$threshold <- NA

  if(NROW(xlines$maxlist)<1 | NROW(xlines$minlist)<1){
    warning("no envelopes found, try changing the tolerance, increasing n or pscore or reducing pfit")
    return(tchannel)
  }

  maxx <- as.numeric(xlines$maxlines[[1]][NROW(xlines$maxlines[[1]])])
  minx <- as.numeric(xlines$minlines[[1]][NROW(xlines$minlines[[1]])])

  if(maxx < minx)return(tchannel)

  vol <- sd(na.omit(TTR::ROC(quantmod::Cl(x)))); tol <- 0.25*vol

  idx <- max(zoo::index(xlines$maxlines[[1]])[1],zoo::index(xlines$minlines[[1]])[1])
  max0 <- as.numeric(xlines$maxlines[[1]][idx])
  min0 <- as.numeric(xlines$minlines[[1]][idx])
  startdev <- 100*(max0/min0-1)
  enddev <- 100*(maxx/minx-1)

  if(startdev > enddev & abs(startdev-enddev)>(100*tol)){
    tchannel$name <- "triangle"
  } else if(startdev < enddev & abs(startdev-enddev)>(100*tol)){
    tchannel$name <- "megaphone"
  } else{
    tchannel$name <- "channel"
  }

  startmean <- mean(min0,max0)
  endmean <- mean(minx,maxx)

  if(startmean < endmean & abs(endmean/startmean-1)>tol ){
    if(tchannel$name == "triangle"){
      tchannel$type <- "continuation"
      tchannel$dir <- 1
    }else if(tchannel$name == "megaphone"){
      tchannel$type <- "breakout"
      tchannel$dir <- -1
    } else{
      tchannel$type <- "neutral"
      tchannel$dir <- 0
    }
  } else if(startmean > endmean & abs(endmean/startmean-1)>tol){
    if(tchannel$name == "triangle"){
      tchannel$type <- "continuation"
      tchannel$dir <- -1
    }else if(tchannel$name == "megaphone"){
      tchannel$type <- "breakout"
      tchannel$dir <- 1
    } else{
      tchannel$type <- "neutral"
      tchannel$dir <- 0
    }
  } else{
    tchannel$type <- "neutral"
    tchannel$dir <- 0
  }

  if(tchannel$dir==1)tchannel$threshold <- maxx
  if(tchannel$dir==-1)tchannel$threshold <- minx

  class(tchannel) <- "tchannel"

  return(tchannel)

}

#'@export
print.tchannel <- function(x,...){
  cat(paste("name:",x$name)); cat("\n")
  cat(paste("type:",x$type)); cat("\n")
  cat(paste("direction:",x$dir)); cat("\n")
  cat(paste("threshold:",round(x$threshold,3)))
}
