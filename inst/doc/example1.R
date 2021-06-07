## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SEAGLE)

## -----------------------------------------------------------------------------
y_loc <- system.file("extdata", "y.txt", package = "SEAGLE")
y <- as.numeric(unlist(read.csv(y_loc)))

X_loc <- system.file("extdata", "X.txt", package = "SEAGLE")
X <- as.matrix(read.csv(X_loc))

E_loc <- system.file("extdata", "E.txt", package = "SEAGLE")
E <- as.numeric(unlist(read.csv(E_loc)))

G_loc <- system.file("extdata", "G.txt", package = "SEAGLE")
G <- as.matrix(read.csv(G_loc))

## -----------------------------------------------------------------------------
objSEAGLE <- prep.SEAGLE(y=as.matrix(y), X=X, intercept=1, E=E, G=G)

## -----------------------------------------------------------------------------
res <- SEAGLE(objSEAGLE, init.tau=0.5, init.sigma=0.5)
res$T
res$pv

