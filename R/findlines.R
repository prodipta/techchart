# cpt trend functions
#'@importFrom Rcpp evalCpp
#'@useDynLib techchart

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

#' Merge lines as computed from a find.lines call
#' @param xlines a findlines object
#' @param tolerance Maximum number nearby points to merge
#' @seealso \code{\link[techchart]{find.lines}}
#' @seealso \code{\link[techchart]{sanitize.lines}}
#' @export
mergelines <- function(xlines, tolerance=2){

  # cluster along start and end lines
  d <- unique(xlines$end); cluster <- hclust(dist(d))
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

#' Apply filters to lines as computed from a find.lines or merge.lines call
#' @param xlines a findlines object
#' @param ptheta cosine of the angle of near-vertical lines to filter out
#' @param pscore maximum percentile of the scores to keep
#' @param pfit minimum (unadjusted) R-square value to filter the fit
#' @param nlines maximum number of lines to keep
#' @seealso \code{\link[techchart]{find.lines}}
#' @seealso \code{\link[techchart]{mergelines}}
#' @export
sanitize.lines <- function(xlines, ptheta=0.9, pscore=0.3, pfit=0.95, nlines=10){

  xlines <- xlines[abs(cos(xlines$theta*pi/180))<ptheta,]
  if(NROW(xlines)<nlines)return(xlines)

  cutoffscore <- as.numeric(quantile(unique(xlines$score), pscore))
  xlines <- xlines[xlines$score < pscore,]
  if(NROW(xlines)<nlines)return(xlines)

  xlines <- xlines[xlines$fit > pfit,]
  if(NROW(xlines)<nlines)return(xlines)

  xlines <- xlines[with(xlines, order(score,-fit)),]
  xlines <- xlines[1:nlines,]

  return(xlines)
}

#' Find major lines from a set of given extremum points
#' @param x a dataframe, with column names x and y
#' @param flag set to 1 for finding handling maximum points and -1 for minimum points
#' @param r.tol minimum step for r coordinate quantization
#' @param theta.tol minimum step for theta coordinate quantization (degrees)
#' @seealso \code{\link[techchart]{sanitize.lines}}
#' @seealso \code{\link[techchart]{mergelines}}
#' @export
find.lines <- function(x, flag, r.tol=0.02, theta.tol=2){

  rbucket <- seq(0,1.42,r.tol); abucket <- seq(0,360,theta.tol)
  xlines <- houghtransform(range01(x$x),range01(x$y),flag,rbucket, abucket)
  xlines$r <- rbucket[xlines$r+1]; xlines$theta <- abucket[xlines$theta+1]
  xlines <- xlines[xlines$r != 0,]
  return(xlines)
}

