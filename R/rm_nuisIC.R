#' Check \code{Q2_max}
#'
#' Check \code{Q2_max} and set it if \code{NULL}.
#'
#' @param Q2_max,nQ,nT The args
#' @return \code{Q2_max}, clamped to acceptable range of values.
#' @keywords internal
Q2_max_check <- function(Q2_max, nQ, nT){
  if (!is.null(Q2_max)) {
    if (round(Q2_max) != Q2_max || Q2_max <= 0) {
      stop('`Q2_max` must be `NULL` or a non-negative integer.')
    }
  } else {
    Q2_max <- pmax(round(nT*.50 - nQ), 1)
  }

  # This is to avoid the area of the pesel objective function that spikes close
  #   to rank(X), which often leads to nPC close to rank(X)
  if (Q2_max > round(nT*.75 - nQ)) {
    warning('`Q2_max` too high, setting to 75% of T.')
    Q2_max <- round(nT*.75 - nQ)
  }

  Q2_max
}

#' Remove nuisance ICs from data
#'
#' Subtract estimated nuisance ICs from data matrix. If the number of 
#'  nuisance ICs is not provided, will estimate using PESEL (the nuisance ICs
#'  are a rotation of the nuisance PCs).
#' 
#' @param BOLD the row-centered \eqn{V} by \eqn{T} data
#' @param DR,template_mean We need an initial estimate of the group ICs, so 
#'  that we can avoid removing them during removal of the estimated noise ICs.
#'  Provide either \code{DR} if dual regression has already been calculated, or
#'  \code{template_mean} (in this context, equivalent to a GICA result) if it 
#'  hasn't. Exactly one must be provided.
#' @param Q2 The number of nuisance ICs. If \code{NULL} (default) will estimate 
#'  using PESEL.
#' @param Q2_max If \code{Q2} is \code{NULL}, PESEL's estimate will be less than
#'  or equal to \code{Q2_max}. If \code{Q2_max} is \code{NULL} (default), do not
#'  limit PESEL's estimate. 
#' @param checkRowCenter Check row means, and raise an error if they are 
#'  nonzero? Default: \code{TRUE}.
#' @param verbose If \code{TRUE}, display progress updates.
#' @param return_Q2 Return (estimated) \code{Q2} too? Default: \code{FALSE}.
#' 
#' @return The \eqn{V} by \eqn{T} data with the estimated nuisance ICs 
#'  subtracted from it. If \code{return_Q2}, a list of length two: the second 
#'  entry will be \code{Q2}.
#' 
#' @importFrom pesel pesel
#' @importFrom fMRItools colCenter dual_reg
#' @keywords internal 
rm_nuisIC <- function(BOLD, DR=NULL, template_mean=NULL, Q2=NULL, Q2_max=NULL, 
  checkRowCenter=TRUE, verbose=FALSE, return_Q2=FALSE){

  stopifnot(is.matrix(BOLD))
  if (checkRowCenter) { stopifnot(all(rowMeans(BOLD) < 1e-8)) }

  # Get `nQ` and check that `Q2` makes sense.
  if (is.null(DR)) {
    if (is.null(template_mean)) { stop("Need either `DR` or `template_mean`.") }
    nQ <- ncol(template_mean)
  } else {
    nQ <- ncol(DR$A)
  }
  Q2_max <- Q2_max_check(Q2_max, nQ=nQ, nT=ncol(BOLD))

  if ( (!is.null(Q2) && Q2==0) || (!is.null(Q2_max) && Q2_max==0) ) {
    # if (verbose) { cat("`Q2` and/or `Q2_max` specifies no nuisance ICs. Skipping denoising.\n") }
    if (return_Q2) {
      return(list(BOLD=BOLD, Q2=0))
    } else {
      return(BOLD)
    }
  }

  # i. PERFORM DUAL REGRESSION TO GET INITIAL ESTIMATE OF TEMPLATE ICS
  if (is.null(DR)) { DR <- dual_reg(BOLD, template_mean) }

  # ii. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD2
  BOLD2 <- BOLD - t(DR$A %*% DR$S) #data without template ICs

  # Remove the global signal, for the purpose of estimating the nuisance ICs.
  #   It will be "added back" later when we subtract the nuisance ICs estimate
  #   from the original `BOLD`.
  BOLD2 <- fMRItools::colCenter(BOLD2)

  # iii. ESTIMATE THE NUMBER OF REMAINING ICS
  #   pesel function expects nxp data and will determine asymptotic framework
  #   here, we consider n=V (vertices) and p=T (timepoints). (it will use n-asymptotic framework)
  if (is.null(Q2)) {
    if(verbose) cat(paste0('Estimating number of nuisance components... '))
    Q2 <- suppressWarnings(pesel::pesel(BOLD2, npc.max=Q2_max, method='homogenous')$nPCs) #estimated number of nuisance ICs
    if(verbose) cat(paste0(Q2,'\n'))
  }

  # Not sure if this actually happens?
  if (Q2 == 0) {
    if (return_Q2) {
      return(list(BOLD=BOLD, Q2=Q2))
    } else {
      return(BOLD)
    }
  }

  # iv. ESTIMATE THE NUISANCE ICS USING GIFT/INFOMAX
  #   if(verbose) cat(paste0('ESTIMATING AND REMOVING ',Q2,' NUISANCE COMPONENTS\n'))
  #   ICA_BOLD2 <- icaimax(BOLD2, nc=Q2, center=TRUE)
  #   fit <- ICA_BOLD2$M %*% t(ICA_BOLD2$S)

  # iv. INSTEAD OF ESTIMATING ICS, JUST ESTIMATE PCS!
  # v. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA
  #   THE RESIDUAL IS THE EXACT SAME BECAUSE THE ICS ARE JUST A ROTATION OF THE PCS
  #   IF THE NUISANCE ICS ARE NOT OF INTEREST, CAN TAKE THIS APPROACH
  #   Y2 * Y2' = U * D^2 * U'
  #   V' = (1/D) * U' * Y2
  #   UDV' = U * U' * Y2
  if (verbose) { cat('Estimating & subtracting nuisance components.') }
  BOLD2 <- BOLD - (BOLD2 %*% tcrossprod(svd(crossprod(BOLD2), nu=Q2, nv=0)$u))

  if (return_Q2) {
    return(list(BOLD=BOLD2, Q2=Q2))
  } else {
    return(BOLD2)
  }
}