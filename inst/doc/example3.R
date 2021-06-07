## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SEAGLE)

## -----------------------------------------------------------------------------
acr_loc <- system.file("extdata", "ACR.raw", package = "SEAGLE")
acr <- read.delim(acr_loc, sep=" ")

## -----------------------------------------------------------------------------
dim(acr)
head(acr)

## -----------------------------------------------------------------------------
G <- as.matrix(acr[,-c(1:6)])
dim(G)

## -----------------------------------------------------------------------------
# Determine number of individuals and loci
n <- dim(G)[1]
L <- dim(G)[2]

# Generate synthetic phenotype and covariate data
set.seed(1)
y <- 2 * rnorm(n)

set.seed(2)
X <- rnorm(n)

set.seed(3)
E <- rnorm(n)

## -----------------------------------------------------------------------------
objSEAGLE <- prep.SEAGLE(y=y, X=X, intercept=0, E=E, G=G)

## -----------------------------------------------------------------------------
res <- SEAGLE(objSEAGLE, init.tau=0.5, init.sigma=0.5)
res$T
res$pv

