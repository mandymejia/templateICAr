#' Activations of (s)tICA
#'
#' Identify areas of activation in each independent component map
#'
#' @param tICA Fitted (spatial) template ICA object from \code{\link{templateICA}}.
#' @param u Activation threshold. Default: \code{0}.
#' @param alpha Significance level for joint PPM. Default: \code{0.01}.
#' @param type Type of region.  Default: \code{">"} (positive excursion region).
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
#' @return A list containing activation maps for each IC and the joint and
#'  marginal PPMs for each IC. If the input represented CIFTI- or NIFTI-format
#'  data, then the activations maps will be formatted accordingly.
#'
#'  Use \code{summary} to obtain information about the activations results.
#'  For CIFTI-format activations, use \code{plot} to visualize the activation
#'  maps.
#'
#' @export
#'
#' @importFrom excursions excursions
#' @importFrom stats pnorm p.adjust
#'
#' @examples
#' \dontrun{
#'  activations(tICA_result, alpha=.05, deviation=TRUE)
#' }
activations <- function(
  tICA, u=0, alpha=0.01, type=">", method_p='BH',
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

  if(!(type %in% c('>','<','!='))) stop("type must be one of: '>', '<', '!='")
  if(alpha <= 0 | alpha >= 1) stop('alpha must be between 0 and 1')

  L <- ncol(tICA$theta_MLE$A)
  if(is.null(which.ICs)) which.ICs <- 1:L
  if(min((which.ICs) %in% (1:L))==0) stop('Invalid entries in which.ICs')

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
    tICA$subjICmean <- matrix(tICA$subjICmean[rep(mask_nii, L)], ncol=L)
    tICA$subjICse <- matrix(tICA$subjICse[rep(mask_nii, L)], ncol=L)
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
  if (is_stICA) {

    if(verbose) cat('Determining areas of activations based on joint posterior distribution of latent fields\n')

    nvox <- nrow(tICA$s_mean)

    #identify areas of activation in each IC
    active <- jointPPM <- marginalPPM <- vars <- matrix(NA, nrow=nvox, ncol=L)

     for(q in which.ICs){
      if(verbose) cat(paste0('.. IC ',q,' (',which(which.ICs==q),' of ',length(which.ICs),') \n'))
      inds_q <- (1:nvox) + (q-1)*nvox
      if(deviation){
        Dinv_mu_s <-  (as.vector(tICA$s_mean) - as.vector(tICA$t_mean) - u)/as.vector(sqrt(tICA$t_var))
      } else {
        Dinv_mu_s <- (as.vector(tICA$s_mean) - u)/as.vector(sqrt(tICA$t_var))
      }

      if(q==which.ICs[1]) {
        #we scale mu by D^(-1) to use Omega for precision (actual precision of s|y is D^(-1) * Omega * D^(-1) )
        #we subtract u first since rescaling by D^(-1) would affect u too
        #save rho from first time running excursions, pass into excursions for other ICs
        tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = Q, type = type, u = 0, ind = inds_q))
        if(verbose) print(tmp)
        rho <- res_q$rho
      } else {
        tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = Q, type = type, u = 0, ind = inds_q, rho=rho))
        if(verbose) print(tmp)
      }
      active[,q] <- res_q$E[inds_q]
      jointPPM[,q] <- res_q$F[inds_q]
      marginalPPM[,q] <- res_q$rho[inds_q]
      vars[,q] <- res_q$vars[inds_q]
    }

    result <- list(
      active=active, jointPPM=jointPPM, marginalPPM=marginalPPM,
      vars=vars, u = u, alpha = alpha, type = type, deviation=deviation
    )
  }

  if (is_tICA) {

    if(verbose) cat('Determining areas of activations based on hypothesis testing at each location\n')

    nvox <- nrow(tICA$s_mean)
    L <- ncol(tICA$s_mean)
    if(deviation){
      t_stat <- (as.matrix(tICA$s_mean) - tICA$t_mean - u) / as.matrix(tICA$s_se)
    } else {
      t_stat <- (as.matrix(tICA$s_mean) - u) / as.matrix(tICA$s_se)
    }

    if(type=='>') pvals <- 1-pnorm(t_stat)
    if(type=='<') pvals <- pnorm(t_stat)
    if(type=='!=') pvals <- 2*(1-pnorm(abs(t_stat)))

    if(verbose) cat(paste0('Correcting for multiple comparisons with method ', method_p))

    pvals_adj <- active <- matrix(NA, nrow=nvox, ncol=L)
    if(is.null(method_p)) method_p <- 'none'
    for(q in which.ICs){
      pvals_adj[,q] <- p.adjust(pvals[,q], method=method_p)
      active[,q] <- (pvals_adj[,q] < alpha)
    }

    result <- list(
      active = active,
      pvals = pvals, pvals_adj = pvals_adj,
      se = tICA$s_se, tstats = t_stat,
      alpha = alpha, method_p = method_p,
      type=type, u = u, deviation=deviation
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
