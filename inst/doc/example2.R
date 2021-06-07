## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SEAGLE)

## -----------------------------------------------------------------------------
dat <- makeSimData(H=cosihap, n=5000, L=100, gammaG=1, gammaGE=0, causal=40, seed=1)


## -----------------------------------------------------------------------------
objSEAGLE <- prep.SEAGLE(y=dat$y, X=dat$X, intercept=1, E=dat$E, G=dat$G)

## -----------------------------------------------------------------------------
res <- SEAGLE(objSEAGLE, init.tau=0.5, init.sigma=0.5)
res$T
res$pv

