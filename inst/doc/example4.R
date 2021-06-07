## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SEAGLE)

## -----------------------------------------------------------------------------
# Read in gene list
dir <- "../inst/extdata/"
genelist <- read.delim(paste0(dir, "glist-hg19_chr22_example"), sep=" ", header=FALSE)

# Identify number of genes in genelist
num_genes <- dim(genelist)[1]

# Generate synthetic phenotype and covariate data
n <- 1092 # number of study participants

set.seed(1)
y <- 2 * rnorm(n)

set.seed(2)
X <- as.matrix(rnorm(n))

set.seed(3)
E <- as.matrix(rnorm(n))

## -----------------------------------------------------------------------------
# Initialize output containers for T and p-value for each gene
T.list <- numeric(length=num_genes)
pv.list <- numeric(length=num_genes)

# Run SEAGLE on each gene in genelist
for (i in 1:num_genes) {
  
  # Identify current gene
  gene_name <- genelist[i,4]
  
  # Obtain G
  Gtemp <- read.delim(paste0(dir, gene_name, ".raw"), sep=" ")
  G <- as.matrix(Gtemp[,-c(1:6)])
  L <- dim(G)[2]                ## number of SNPs
  
  # Make weights
  avg_newsnp <- colMeans(G)
  fA  = avg_newsnp/2            ## freq of allele "A"
  fa  = 1-fA                    ## freq of allele "a"
  maf = fA; maf[fA>0.5]= fa[fA>0.5]## maf should be b/w 0 and 0.5
  G = G[ ,maf>0]                ## only keep those SNPs with MAF>0
  maf <- maf[maf > 0]           ## Keep only MAF > 0)
  wt   = sqrt(maf^(-3/4))       ## Note we take the square root here
  if (length(wt) > 1) {
    G_final    = G %*% diag(wt) ## Input this G
  } else {
    T.list[i] <- NA
    pv.list[i] <- NA
    next
  }

  # Run SEAGLE
  objSEAGLE <- prep.SEAGLE(y=y, X=X, intercept=0, 
                           E=E, G=G_final)
  res <- SEAGLE(objSEAGLE, init.tau=0.5, init.sigma=0.5)
  
  # Save T and p-value into output lists
  T.list[i] <- res$T
  pv.list[i] <- res$pv
}

## -----------------------------------------------------------------------------
resMat <- cbind(T.list, pv.list)
colnames(resMat) <- c("T", "p-value")
resMat

