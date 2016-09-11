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
                     if(!trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     (R[1] > R[2]) & (R[5] > R[4]) & (R[3] > R[1]) & (R[3] > R[5]) &
                     abs(R[1] - (R[1]+R[5])/2) < tolerance * (R[1]+R[5])/2 &
                     abs(R[5] - (R[1]+R[5])/2) < tolerance * (R[1]+R[5])/2 &
                     abs(R[2] - (R[2]+R[4])/2) < tolerance * (R[2]+R[4])/2 &
                     abs(R[4] - (R[2]+R[4])/2) < tolerance * (R[2]+R[4])/2
                     ")
  length <- 5
  start <- 1
  trigger <- quote((R[2]+R[4])/2)
  patterns$HS <- list(call=call.expr, length=length, start=start,
                      trigger=trigger)

  ######################### Inverse Head and Shoulder ####################
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     (R[1] < R[2]) & (R[5] < R[4]) & (R[3] < R[1]) & (R[3] < R[5]) &
                     abs(R[1] - (R[1]+R[5])/2) < tolerance * (R[1]+R[5])/2 &
                     abs(R[5] - (R[1]+R[5])/2) < tolerance * (R[1]+R[5])/2 &
                     abs(R[2] - (R[2]+R[4])/2) < tolerance * (R[2]+R[4])/2 &
                     abs(R[4] - (R[2]+R[4])/2) < tolerance * (R[2]+R[4])/2
                     ")
  length <- 5
  start <- -1
  trigger <- quote((R[2]+R[4])/2)
  patterns$IHS <- list(call=call.expr, length=length, start=start,
                       trigger=trigger)

  ######################### Broadening Top ####################
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     R[1] < R[3] & R[3] < R[5] & R[2] > R[4]")
  length <- 5
  start <- 1
  trigger <- quote((R[2]+R[4])/2)
  patterns$BTOP <- list(call=call.expr, length=length, start=start,
                        trigger=trigger)

  ######################### Broadening Bottom ####################
  call.expr <- parse(text="
                     fit <- lm(c(E1,E2,E3,E4,E5)~c(T1,T2,T3,T4,T5));
                     R <- fit$residuals; R <- R+E1;
                     if(!trend.adjusted){R <- c(E1,E2,E3,E4,E5)};
                     R[1] > R[3] & R[3] > R[5] & R[2] < R[4]")
  length <- 5
  start <- -1
  trigger <- quote((R[2]+R[4])/2)
  patterns$BBOT <- list(call=call.expr, length=length, start=start,
                        trigger=trigger)

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
