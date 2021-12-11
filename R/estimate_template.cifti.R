#' Estimate CIFTI template from DR estimates
#'
#' Estimate template for Template or Diagnostic ICA based on CIFTI-format data
#'
#' @param cifti_fnames Vector of file paths of CIFTI-format fMRI timeseries
#'  (*.dtseries.nii) for template estimation
#' @param cifti_fnames2 Vector of file paths of "retest" CIFTI-format fMRI
#'  timeseries (*.dtseries.nii) for template estimation.  Must be from the same
#'  subjects and in the same order as cifti_fnames.  If none specified, will
#'  create pseudo test-retest data from single session.
#' @param GICA_fname File path of CIFTI-format group ICA maps (ending in .d*.nii)
#' @param inds Numeric indices of the group ICs to include in template. If NULL,
#'  use all group ICs.
#' 
#'  The dual regression estimate (as well as the cleaned dual regression estimate)
#'  will be calculated using all ICs. ICs not in \code{inds} are removed afterward.
#'  This is because removing the ICs prior to dual regression would leave unmodeled
#'  signals, which could bias the dual regression results (and thus the templates).
#' @param scale Logical indicating whether BOLD data should be scaled by the
#'  spatial standard deviation before template estimation. Will also scale each
#'  component in the ICA maps. Default: \code{TRUE}.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param keep_DR Keep the DR estimates? If \code{TRUE}, the DR estimates are
#'  returned too. If a single file path, save the DR estimates as an .rds file at
#'  that location. If \code{FALSE} (default), do not save the DR estimates and
#'  only return the templates.
#   If a list of two vectors of file paths, each of which the
#   length of \code{cifti_fnames(2)}, save the DR estimates as individual .rds
#   files at these locations.
#' @param Q2 The number of nuisance ICs to identify. If NULL, will be estimated. Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T). Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param var_method \code{"unbiased"} (default) or \code{"non-negative"}
#' @param verbose If \code{TRUE}, display progress updates
#' @param out_fname (Required if templates are to be resampled to a lower spatial
#' resolution, usually necessary for spatial template ICA.) The path and base name
#' prefix of the CIFTI files to write. Will be appended with "_mean.dscalar.nii" for
#' template mean maps and "_var.dscalar.nii" for template variance maps.
#' @param FC If TRUE, include template for functional connectivity
#'
#' @importFrom ciftiTools read_cifti write_cifti newdata_xifti
#' @importFrom stats cov quantile
#'
#' @return List of two elements: template mean of class xifti and template
#'  variance of class xifti
#'
#' @export
#'
estimate_template.cifti <- function(
  cifti_fnames,
  cifti_fnames2=NULL,
  GICA_fname,
  inds=NULL,
  scale=TRUE,
  brainstructures=c("left","right"),
  var_method=c("unbiased", "non-negative"),
  keep_DR=FALSE,
  Q2=NULL, maxQ=NULL,
  verbose=TRUE,
  out_fname=NULL,
  FC = FALSE){

  #TO DOs:
  # Create function to print and check template, template_cifti and template_nifti objects


  if (is.character(cifti_fnames)) {
    notthere <- sum(!file.exists(cifti_fnames))
    if(notthere == length(cifti_fnames)) stop('The files in cifti_fnames do not exist.')
    if(notthere > 0) warning(paste0('There are ', notthere, ' files in cifti_fnames that do not exist. These scans will be excluded from template estimation.'))
  }
  if (retest && is.character(cifti_fnames)) {
    notthere2 <- sum(!file.exists(cifti_fnames2))
    if(notthere2 == length(cifti_fnames2)) stop('The files in cifti_fnames2 do not exist.')
    if(notthere2 > 0) warning(paste0('There are ', notthere2, ' files in cifti_fnames2 that do not exist. These scans will be excluded from template estimation.'))
  }

  # Check arguments.
  if (!is.logical(scale) || length(scale) != 1) { stop('scale must be a logical value') }
  brainstructures <- match_input(
    brainstructures, c("left","right","subcortical","all"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }

  # Read GICA result
  if(verbose) { cat('\n Reading in GICA result') }
  GICA <- read_cifti(GICA_fname, brainstructures=brainstructures)
  GICA <- newdata_xifti(GICA, scale(as.matrix(GICA), scale=FALSE))
  V <- nrow(GICA); Q <- ncol(GICA)
  # Center each IC map.
  if(verbose) {
    cat(paste0('\n Number of data locations: ',V))
    cat(paste0('\n Number of original group ICs: ',Q))
  }

  L <- Q
  if(!is.null(inds)){
    if(any(!(inds %in% 1:Q))) stop('Invalid entries in inds argument.')
    L <- length(inds)
  } else {
    inds <- 1:Q
  }

  N <- length(cifti_fnames)

  if(verbose){
    cat(paste0('\n Number of template ICs: ',L))
    cat(paste0('\n Number of training subjects: ',N))
  }

  # PERFORM DUAL REGRESSION ON (PSEUDO) TEST-RETEST DATA
  DR1 <- DR2 <- array(NA, dim=c(N, L, V))
  if(FC) FC1 <- FC2 <- array(NA, dim=c(N, L, L)) #for functional connectivity template

  missing_data <- NULL

  for(ii in 1:N){
    if(verbose) cat(paste0('\n Reading and analyzing data for subject ',ii,' of ',N, ".\n"))
    if (retest) { BOLD2 <- cifti_fnames2[ii] } else { BOLD2 <- NULL }

    DR_ii <- dual_reg2(
      cifti_fnames[ii], BOLD2=BOLD2, GICA=as.matrix(GICA),
      scale=scale, format="CIFTI", dim_expect=V,
      brainstructures=brainstructures, 
      Q2=Q2, maxQ=maxQ, verbose=verbose
    )

    if (!isFALSE(DR_ii$missing)) {
      missing_data <- c(missing_data, DR_ii$missing)
      next
    }

    DR1[ii,,] <- DR_ii$test[inds,]
    DR2[ii,,] <- DR_ii$retest[inds,]
    if(FC) FC1[ii,,] <- cov(DR_ii$test[,inds])
    if(FC) FC2[ii,,] <- cov(DR_ii$retest[,inds])
  }

  # Estimate template
  if (verbose) { cat("Estimating template.\n") }
  var_method <- match.arg(var_method, c("unbiased", "non-negative"))
  template <- estimate_template_from_DR(DR1, DR2, var_method=var_method)

  if(FC){

    mean_FC <- (apply(FC1, c(2,3), mean, na.rm=TRUE) + apply(FC2, c(2,3), mean, na.rm=TRUE))/2
    var_FC_tot  <- (apply(FC1, c(2,3), var, na.rm=TRUE) + apply(FC2, c(2,3), var, na.rm=TRUE))/2
    var_FC_within  <- 1/2*(apply(FC1-FC2, c(2,3), var, na.rm=TRUE))
    var_FC_between <- var_FC_tot - var_FC_within
    var_FC_between[var_FC_between < 0] <- NA

    #function to minimize w.r.t. k
    fun <- function(nu, p, var_ij, xbar_ij, xbar_ii, xbar_jj){
      LHS <- var_ij
      phi_ij <- xbar_ij*(nu-p-1)
      phi_ii <- xbar_ii*(nu-p-1)
      phi_jj <- xbar_jj*(nu-p-1)
      RHS_numer <- (nu-p+1)*phi_ij^2 + (nu-p-1)*phi_ii*phi_jj
      RHS_denom <- (nu-p)*((nu-p-1)^2)*(nu-p-3)

      sq_diff <- (LHS - RHS_numer/RHS_denom)^2
      return(sq_diff)
    }

    nu_est <- matrix(NA, L, L)
    for(q1 in 1:L){
      for(q2 in q1:L){

        #estimate k = nu - p - 1
        nu_opt <- optimize(f=fun, interval=c(L+1,L*10), p=L, var_ij=var_FC_between[q1,q2], xbar_ij=mean_FC[q1,q2], xbar_ii=mean_FC[q1,q1], xbar_jj=mean_FC[q2,q2])
        nu_est[q1,q2] <- nu_opt$minimum
      }
    }
    nu_est[lower.tri(nu_est, diag=FALSE)] <- NA
    nu_est1 <- quantile(nu_est[upper.tri(nu_est, diag=TRUE)], 0.1, na.rm = TRUE)

    template_FC <- list(nu = nu_est1,
                        psi = mean_FC*(nu_est1 - L - 1))
  } else {
    template_FC <- NULL
  }

  # Keep DR
  if (!isFALSE(keep_DR)) {
    if (is.character(keep_DR)) {
      if (length(keep_DR) > 1) {
        warning("Using first entry of `keep_DR`.")
        keep_DR <- keep_DR[1]
      }
      if (!endsWith(keep_DR, ".rds")) { keep_DR <- paste0(keep_DR, ".rds") }
      saveRDS(list(DR1=DR1, DR2=DR2), keep_DR)
      keep_DR <- FALSE # no longer need it.
    } else if (!isTRUE(keep_DR)) {
      warning("`keep_DR` should be `TRUE`, `FALSE`, or a file path. Using `FALSE`.")
      keep_DR <- FALSE
    }
  }
  if (!keep_DR) { rm(DR1, DR2) }

  # Format template as "xifti"s
  GICA <- select_xifti(GICA, inds)
  GICA$meta$cifti$names <- paste0("IC ", inds)
  xifti_mean <- newdata_xifti(GICA, template$mean)
  xifti_mean$meta$cifti$misc <- list(template="mean")
  xifti_var <- newdata_xifti(GICA, template$var)
  xifti_var$meta$cifti$misc <- list(template="var")

  if(!is.null(out_fname)){
    write_cifti(xifti_mean, paste0(out_fname, '_mean.dscalar.nii'), verbose=verbose)
    write_cifti(xifti_var, paste0(out_fname, '_var.dscalar.nii'), verbose=verbose)
  }

  result <- list(
    template_mean=xifti_mean, template_var=xifti_var, template_FC=template_FC,
    scale=scale, inds=inds, var_method=var_method
  )
  if (keep_DR) { result <- c(result, list(DR=list(DR1=DR1, DR2=DR2))) }

  class(result) <- 'template_cifti'
  result
}

#' Summarize a \code{"template_cifti"} object
#'
#' Summary method for class \code{"template_cifti"}
#'
#' @param object Object of class \code{"template_cifti"}. 
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary template_cifti
summary.template_cifti <- function(object, ...) {

  x <- c(
    summary(object$template_mean),
    list(has_DR="DR" %in% names(object)),
    object[c("scale", "inds", "var_method")]
  )

  class(x) <- "summary.template_cifti"
  return(x)
}

#' @rdname summary.template_cifti
#' @export
#' 
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.template_cifti
print.summary.template_cifti <- function(x, ...) {

  cat("====TEMPLATE INFO====================\n")
  cat("Spatially scaled:", x$scale, "\n")
  cat("Variance method: ", x$var_method, "\n")
  cat("\n")

  class(x) <- "summary.xifti"
  print(x) 
}

#' @rdname summary.template_cifti
#' @export
#' 
#' @method print template_cifti
print.template_cifti <- function(x, ...) {
  print.summary.template_cifti(summary(x))
}

#' Plot template
#' 
#' @param x The template from \code{estimate_template.cifti}
#' @param stat \code{"mean"} (default) or \code{"var"}
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @importFrom ciftiTools view_xifti
#' @method plot template_cifti
plot.template_cifti <- function(x, stat=c("mean", "var"), ...) {
  stopifnot(inherits(x, "template_cifti"))
  stat <- match.arg(stat, c("mean", "var"))
  view_xifti(x[[paste0("template_", stat)]], ...)
}