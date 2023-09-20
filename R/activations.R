#' Activations of (spatial) template ICA
#'
#' Identify areas of activation in each independent component map from the
#'  result of (spatial) template ICA.
#'
#' @param tICA Fitted (spatial) template ICA object from \code{\link{templateICA}}.
#' @param u,z Set a threshold value for activation? A threshold value can be
#'  specified directly with \code{u}, or a z-score-like threshold in terms of
#'  standard deviations (the SD of values in the mean template) can be specified
#'  with \code{z}. Only one type of threshold can be used. Default: \code{NULL}
#'  (do not use a threshold). Either argument can also be a vector of the same
#'  length as \code{which.ICs}, to use different thresholds for each IC.
#' @param alpha Significance level for joint PPM. Default: \code{0.01}.
#' @param type Type of region: \code{">"}, \code{"<"}, \code{"!="}, or 
#'  \code{"abs >"}. The last case tests for magnitude by taking the absolute
#'  value and then testing if they are greater than... Default: \code{"abs >"}.
#' @param method_p If the input is a \code{"tICA.[format]"} model object, the type of
#'  multiple comparisons correction to use for p-values, or \code{NULL} for no
#'  correction. See \code{help(p.adjust)}. Default: \code{"BH"} (Benjamini &
#'  Hochberg, i.e. the false discovery rate).
#' @param verbose If \code{TRUE}, display progress of algorithm. Default:
#'  \code{FALSE}.
#' @param which.ICs Indices of ICs for which to identify activations.  If
#'  \code{NULL} (default), use all ICs.
#' @param deviation If \code{TRUE} identify significant deviations from the
#'  template mean, rather than significant areas of engagement. Default:
#'  \code{FALSE}.
#'
#' @return A list containing activation maps for each IC, the joint and
#'  marginal PPMs for each IC, and the parameters used for computing activation.
#'  If the input represented CIFTI- or NIFTI-format data, then the activations
#'  maps will be formatted accordingly.
#'
#'  Use \code{summary} to obtain information about the activations results.
#'  For CIFTI-format activations, use \code{plot} to visualize the activation
#'  maps.
#'
#' @export
#'
#' @importFrom fMRItools is_1 is_integer is_posNum
#' @importFrom excursions excursions
#' @importFrom stats pnorm p.adjust
#'
#' @examples
#' \dontrun{
#'  activations(tICA_result, alpha=.05, deviation=TRUE)
#' }
activations <- function(
  tICA, u=NULL, z=NULL, alpha=0.01,
  type=c("abs >", ">", "<", "!="),
  method_p='BH',
  verbose=FALSE, which.ICs=NULL, deviation=FALSE){

  # Setup ----------------------------------------------------------------------
  is_tICA <- inherits(tICA, "tICA.matrix") || inherits(tICA, "tICA.cifti") || inherits(tICA, "tICA.nifti")
  is_stICA <- inherits(tICA, "stICA.matrix") || inherits(tICA, "stICA.cifti") || inherits(tICA, "stICA.nifti")
  if (!(xor(is_tICA, is_stICA))) { stop("tICA must be of class stICA or tICA") }

  FORMAT <- class(tICA)[grepl("tICA", class(tICA))]
  if (length(FORMAT) != 1) { stop("Not a tICA.") }
  FORMAT <- switch(FORMAT,
    tICA.cifti = "CIFTI",
    tICA.gifti = "GIFTI",
    tICA.nifti = "NIFTI",
    tICA.matrix = "DATA",
    stICA.cifti = "CIFTI",
    stICA.gifti = "GIFTI",
    stICA.nifti = "NIFTI",
    stICA.matrix = "DATA"
  )

  if (FORMAT == "CIFTI") {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to read NIFTI data. Please install it.", call. = FALSE)
    }
  }

  # Simple argument checks
  if ((!is.null(u)) && (!is.null(z))) { stop("Set only one of `u` or `z`.") }
  stopifnot(is.null(u) || is.numeric(u))
  stopifnot(is.null(z) || is.numeric(z))
  stopifnot(fMRItools::is_posNum(alpha))
  stopifnot(alpha <= 1)
  type <- match.arg(type, c("abs >", ">", "<", "!="))
  if(alpha <= 0 | alpha >= 1) stop('alpha must be between 0 and 1')
  stopifnot(fMRItools::is_1(verbose, "logical"))
  stopifnot(fMRItools::is_1(deviation, "logical"))

  nL <- ncol(tICA$A)
  if(is.null(which.ICs)) which.ICs <- seq(nL)
  stopifnot(is.numeric(which.ICs))
  stopifnot(which.ICs == round(which.ICs))
  if(min((which.ICs) %in% (1:nL))==0) stop('Invalid entries in which.ICs')

  if (deviation) {
    if (!is.null(z)) { stop("`z` not compatible with `deviation=TRUE`.") }
    if (u != 0) { warning("`u != 0` not advised for `deviation=TRUE`. Proceeding anyway.") }
  }

  # Get needed metadata from `tICA`.
  Q <- tICA$omega
  mask <- tICA[["mask"]] # avoid grabbing mask_nii

  # Vectorize data.
  if (FORMAT == "CIFTI") {
    xii1 <- ciftiTools::newdata_xifti(ciftiTools::select_xifti(tICA$subjICmean,1), 0)
    tICA$subjICmean <- as.matrix(tICA$subjICmean)
    tICA$subjICse <- as.matrix(tICA$subjICse)
  } else if (FORMAT == "NIFTI") {
    mask_nii <- tICA$mask_nii
    tICA$subjICmean <- matrix(tICA$subjICmean[rep(mask_nii, nL)], ncol=nL)
    tICA$subjICse <- matrix(tICA$subjICse[rep(mask_nii, nL)], ncol=nL)
  }

  tICA <- tICA[c("template_mean", "template_var", "subjICmean", "subjICse")]
  names(tICA) <- c("t_mean", "t_var", "s_mean", "s_se")

  # Apply data mask.
  use_mask <- (!is.null(mask)) && (!all(mask))
  if (use_mask) {
    tICA$s_mean <- tICA$s_mean[mask,]
    tICA$s_se <- tICA$s_se[mask,]
  }

  # Compute activations. -------------------------------------------------------
  u_og <- u
  if (!is.null(z)) {
    stopifnot(length(z) %in% c(1, length(which.ICs)))
    u <- z * sqrt(colVars(tICA$t_mean[,which.ICs,drop=FALSE]))
  } else if (!is.null(u)) {
    stopifnot(length(u) %in% c(1, length(which.ICs)))
  } else {
    u <- 0
  }
  if(length(u)==1) u <- rep(u, length(which.ICs))
  nvox <- nrow(tICA$s_mean)
  u_mat <- matrix(u, nrow=nvox, ncol=nL, byrow = TRUE)

  if (type == "abs >") {
    tICA$s_mean <- abs(tICA$s_mean)
    type <- ">"
  }

  if (is_stICA) {
    if(verbose) cat('Determining areas of activations based on joint posterior distribution of latent fields\n')

    #identify areas of activation in each IC
    active <- jointPPM <- marginalPPM <- vars <- matrix(NA, nrow=nvox, ncol=nL)

     for(q in which.ICs){
      if(verbose) cat(paste0('.. IC ',q,' (',which(which.ICs==q),' of ',length(which.ICs),') \n'))
      inds_q <- (1:nvox) + (q-1)*nvox
      if(deviation){
        Dinv_mu_s <-  (as.vector(tICA$s_mean) - as.vector(tICA$t_mean) - u_mat)/as.vector(sqrt(tICA$t_var))
      } else {
        Dinv_mu_s <- (as.vector(tICA$s_mean) - u_mat)/as.vector(sqrt(tICA$t_var))
      }

      if(q==which.ICs[1]) {
        #we scale mu by D^(-1) to use Omega for precision (actual precision of s|y is D^(-1) * Omega * D^(-1) )
        #we subtract u first since rescaling by D^(-1) would affect u too
        #save rho from first time running excursions, pass into excursions for other ICs
        tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = Q, type = type, u = u[1], ind = inds_q)) #I think u does not matter, should check because may differ across fields
        if(verbose) print(tmp)
        rho <- res_q$rho
      } else {
        tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = Q, type = type, u = u[1], ind = inds_q, rho=rho))
        if(verbose) print(tmp)
      }
      active[,q] <- res_q$E[inds_q]
      jointPPM[,q] <- res_q$F[inds_q]
      marginalPPM[,q] <- res_q$rho[inds_q]
      vars[,q] <- res_q$vars[inds_q]
    }

    if (length(unique(u))==1) { u <- u[1] }
    if (length(unique(z))==1) { z <- z[1] }
    result <- list(
      active=active, jointPPM=jointPPM, marginalPPM=marginalPPM,
      vars=vars, u = u, z = z, alpha = alpha, type = type, deviation=deviation
    )
  }

  if (is_tICA) {

    if(verbose) cat('Determining areas of activations based on hypothesis testing at each location\n')

    nL <- ncol(tICA$s_mean)
    if(deviation){
      t_stat <- (as.matrix(tICA$s_mean) - tICA$t_mean - u_mat) / as.matrix(tICA$s_se)
    } else {
      t_stat <- (as.matrix(tICA$s_mean) - u_mat) / as.matrix(tICA$s_se)
    }

    if(type=='>') pvals <- 1-pnorm(t_stat)
    if(type=='<') pvals <- pnorm(t_stat)
    if(type=='!=') pvals <- 2*(1-pnorm(abs(t_stat)))

    if(verbose) cat(paste0('Correcting for multiple comparisons with method ', method_p))

    pvals_adj <- active <- matrix(NA, nrow=nvox, ncol=nL)
    if(is.null(method_p)) method_p <- 'none'
    for(q in which.ICs){
      pvals_adj[,q] <- p.adjust(pvals[,q], method=method_p)
      active[,q] <- (pvals_adj[,q] < alpha)
    }

    if (length(unique(u))==1) { u <- u[1] }
    if (length(unique(z))==1) { z <- z[1] }
    result <- list(
      active = active,
      pvals = pvals, pvals_adj = pvals_adj,
      se = tICA$s_se, tstats = t_stat,
      alpha = alpha, method_p = method_p,
      type=type, u = u, z = z, deviation=deviation
    )
  }

  # Format result. -------------------------------------------------------------
  # Unmask data.
  if (use_mask) { result$active <- fMRItools::unmask_mat(result$active, mask) }

  # Un-vectorize data.
  if (FORMAT == "CIFTI") {
    result$active <- ciftiTools::newdata_xifti(xii1, as.numeric(result$active))
    result$active <- ciftiTools::move_from_mwall(result$active, -1)
    result$active <- ciftiTools::convert_xifti(
      result$active, "dlabel",
      values=c(-1, 0, 1),
      colors=c("grey", "white", "red"),
      add_white=FALSE
    )

    for (ii in seq(ncol(result$active))) {
      rownames(result$active$meta$cifti$labels[[ii]]) <- c("Medial Wall", "Inactive", "Active") # add "NA" to first
    }
    class(result) <- "tICA_act.cifti"
  } else if (FORMAT == "NIFTI") {
    result$active <- fMRItools::unvec_vol(result$active, mask_nii, fill=NA)
    class(result) <- "tICA_act.nifti"
  } else {
    class(result) <- "tICA_act.matrix"
  }

  result
}
