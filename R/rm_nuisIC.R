#' Remove nuisance ICs from data
#'
#' Subtract estimated nuisance ICs from data matrix. If the number of nuisance ICs is not provided, 
#'  will estimate using PESEL (the nuisance ICs are a rotation of the nuisance PCs).
#' 
#' @param BOLD the \eqn{V} by \eqn{T} data
#' @param DR,template_mean We need an initial estimate of the group ICs, so that we can avoid removing them
#'  during removal of the estimated noise ICs. Provide either \code{DR} if dual regression has already
#'  been calculated, or \code{template_mean} (in this context, equivalent to GICA) if it hasn't. One is required. 
#' @param Q2 the number of nuisance ICs. If \code{NULL} (default) will estimate using PESEL.
#' @param Q2_max If \code{Q2} is \code{NULL}, PESEL's estimate will be less than or equal to \code{Q2_max}.
#'  If \code{Q2_max} is \code{NULL}, do not limit PESEL's estimate. 
#' @param verbose If \code{TRUE}, display progress updates
#' 
#' @return The \eqn{V} by \eqn{T} data with the estimated nuisance ICs subtracted from it
#' 
#' @keywords internal 
rm_nuisIC <- function(BOLD, DR=NULL, template_mean=NULL, Q2=NULL, Q2_max=NULL, verbose=FALSE){

  # Get `nQ` and check that `Q2` makes sense.
  if (is.null(DR)) {
    if (is.null(template_mean)) { stop("Need either `DR` or `template_mean`.") }
    nQ <- ncol(template_mean)
  } else {
    nQ <- ncol(DR$A)
  }
  Q2_max <- Q2_max_check(Q2_max, nQ=nQ, nT=ncol(BOLD))

  if ( (!is.null(Q2) && Q2==0) || (!is.null(Q2_max) && Q2_max==0) ) {
    if (verbose) { cat("`Q2` and/or `Q2_max` specifies no nuisance ICs. Returning original `BOLD`.") }
    return(BOLD)
  }

  #i. PERFORM DUAL REGRESSION TO GET INITIAL ESTIMATE OF TEMPLATE ICS
  if (is.null(DR)) { DR <- dual_reg(BOLD, template_mean) }

  #ii. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD2
  BOLD2 <- BOLD - t(DR$A %*% DR$S) #data without template ICs

  #iii. ESTIMATE THE NUMBER OF REMAINING ICS
  #pesel function expects nxp data and will determine asymptotic framework
  #here, we consider n=T (volumes) and p=V (vertices), and will use p-asymptotic framework
  if (is.null(Q2)) {
    if(verbose) cat(paste0('DETERMINING NUMBER OF NUISANCE COMPONENTS.... '))
    Q2 <- suppressWarnings(pesel(BOLD2, npc.max=Q2_max, method='homogenous')$nPCs) #estimated number of nuisance ICs
    if(verbose) cat(paste0(Q2,'\n'))
  }

  if (Q2 == 0) { return(BOLD) } # Not sure if this actually happens?

  #iv. ESTIMATE THE NUISANCE ICS USING GIFT/INFOMAX
  # if(verbose) cat(paste0('ESTIMATING AND REMOVING ',Q2,' NUISANCE COMPONENTS\n'))
  # ICA_BOLD2 <- icaimax(BOLD2, nc=Q2, center=TRUE)
  # fit <- ICA_BOLD2$M %*% t(ICA_BOLD2$S)

  #iv. INSTEAD OF ESTIMATING ICS, JUST ESTIMATE PCS!
  #v. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD3
  #THE RESIDUAL (BOLD3) IS THE EXACT SAME BECAUSE THE ICS ARE JUST A ROTATION OF THE PCS
  #IF THE NUISANCE ICS ARE NOT OF INTEREST, CAN TAKE THIS APPROACH
  BOLD - (BOLD2 %*% tcrossprod(svd(crossprod(BOLD2), nu=Q2, nv=0)$u))
}