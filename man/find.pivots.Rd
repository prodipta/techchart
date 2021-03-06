% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supports.R
\name{find.pivots}
\alias{find.pivots}
\title{Find supports and resitance for a time series}
\usage{
find.pivots(x, type = c("SR", "FIB"), tolerance = 0.02, strength = 3,
  maxline = 10, lookback = 20)
}
\arguments{
\item{x}{xts object, or vector, representing a time series}

\item{type}{either FIB (Fibonacci) or SR. SR is based on best fit lines of multiple peaks and troughs}

\item{tolerance}{threshold for percentage change or vol multiple for extreme points}

\item{strength}{minimum number of extreme points defining a support}

\item{maxline}{maximum number of support/ resistance lines to return}

\item{lookback}{Used for volatility dependent adaptive threshold}
}
\value{
support/ resistance object (object of class supports)
}
\description{
Find supports and resitance for a time series
}
\examples{
x <- quantmod::getSymbols("^GSPC", auto.assign = FALSE)
x <- x["2015/"]
sups <- find.pivots(x, type = "FIB")
summary(sups)
sups <- find.pivots(x, type = "SR", strength = 5)
summary(sups)
quantmod::chart_Series(x)
quantmod::add_TA(sups$lines[[1]],on=1, lty=2)
}

