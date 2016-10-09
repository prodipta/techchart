## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  fig.align = "center"
)

## ---- fig.show='hold', fig.width = 5-------------------------------------
spx <- quantmod::getSymbols("^GSPC", auto.assign = FALSE)
spx <- spx["2014::2015"]
imppts <- techchart::find.imppoints(spx,2)
head(imppts)
quantmod::chart_Series(spx)
points(as.numeric(imppts$maxima$pos),as.numeric(imppts$maxima$value),bg="green",pch=24,cex=1.25)
points(as.numeric(imppts$minima$pos),as.numeric(imppts$minima$value),bg="red",pch=25,cex=1.25)

## ------------------------------------------------------------------------
spx <- quantmod::getSymbols("^GSPC", auto.assign = FALSE)
spx <- spx["2014::2015"]
cpts <- techchart::find.major.trends(spx)
summary(cpts)
quantmod::chart_Series(spx)
quantmod::add_TA(cpts$segments[[1]],on=1,lty=3, col="red")
quantmod::add_TA(cpts$segments[[2]],on=1,lty=3, col="red")

## ---- fig.show='hold'----------------------------------------------------
spx <- quantmod::getSymbols("^GSPC", auto.assign = FALSE)
spx <- spx["2014::2015"]
sups <- techchart::find.pivots(spx, type = "FIB")
summary(sups)
sups <- techchart::find.pivots(spx, type = "SR", strength = 5)
summary(sups)

## ------------------------------------------------------------------------
spx <- quantmod::getSymbols("^GSPC", auto.assign = FALSE)
spx <- spx["2016-01-01::2016-09-30"]
tchannel <- techchart::find.tchannel(spx,1.25)
tchannel
quantmod::chart_Series(spx)
quantmod::add_TA(tchannel$xlines$maxlines[[1]],on=1, lty=3, col="brown")
quantmod::add_TA(tchannel$xlines$minlines[[1]],on=1, lty=3, col="brown")

