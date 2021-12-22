#' Compute variance decomposition
#' @param x measurements by subjects by variables
#' @keywords internal
#' @return The variance decomposition
var_decomp <- function(x) {

  # Get data dimensions.
  # cat("\tChecking data dimensions and missing values presence.\n")
  d <- dim(x)
  if (length(d) == 2) { x <- array(x, dim=c(dim(x), 1)); d <- dim(x) }
  stopifnot(length(d) == 3)
	nM <- d[1]
  nS <- d[2]

  # Handle missing values: for each var, remove subjects without complete data
  na_mask <- apply(is.na(x), seq(2, 3), any)
  if (any(na_mask)) {
    # message("`NA`s detected. For each variable, removing subjects w/o complete data.\n")
    x[rep(na_mask, each=nM)] <- NA
  }

  # Variance decomposition
  # cat("\tCalculating means.\n")
  sub_mean <- colMeans(x, na.rm=TRUE)
  visit_mean <- apply(x, c(1,3), mean, na.rm=TRUE)
  grand_mean <- colMeans(sub_mean, na.rm=TRUE)
  grand_mean2 <- array(rep(grand_mean, each=nM*nS), dim=dim(x))

  ### Sum of squares
  # cat("\tCalculating variance decomposition.\n")
  SST <- apply((x - grand_mean2)^2, 3, sum, na.rm=TRUE)
  SSW <- apply((x - rep(sub_mean, each=nM))^2, 3, sum, na.rm=TRUE)
  SSB <- nM * colSums((sub_mean - rep(grand_mean, each=nS))^2, na.rm=TRUE)
  SSM <- nS * colSums((visit_mean - rep(grand_mean, each=nM))^2, na.rm=TRUE)
  colnames(SSM) <- colnames(SST)
  SSR <- SSW - SSM # RESIDUAL/ERROR

  list(
    nS = nS,
    nM = nM,
    grand_mean = grand_mean,
    SST = SST,
    SSB = SSB, # delete?
    SSW = SSW, # delete?
    SSM = SSM,
    SSR = SSR
  )
}

#' Compute mean squares from variance decomposition
#' @param vd The variance decomposition
#' @keywords internal
#' @return The mean squares
mean_squares <- function(vd){
  n <- vd$nS; v <- vd$nM
  SST <- vd$SST
  SSW <- vd$SSM + vd$SSR
  SSB <- vd$SST - SSW
  SSM <- vd$SSM
  SSR <- vd$SSR

  ### Mean squares
  MST <- SST / (v*n - 1)
  MSW <- SSW / ((v-1)*(n))
  MSB <- SSB / (n-1)
  MSM <- SSM / (v-1)
  MSR <- SSR / ((v-1)*(n-1))
  # covXY <- (MSB - MSR) / 2 # TRUE for v != 2?

  list(
    MST = MST,
    MSW = MSW,
    MSB = MSB,
    MSM = MSM,
    MSR = MSR
  )
}