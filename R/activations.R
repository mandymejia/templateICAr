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
#'  (do not use a threshold). Either argument can also be a vector to test
#'  multiple thresholds at once, as long as \code{type} is not \code{"!="}
#'  (to ensure the activation regions are successive subsets).
#' @param alpha Significance level for hypothesis testing. Default: \code{0.01}.
#' @param type Type of region: \code{">"} (default), \code{"abs >"}, \code{"<"},
#'  or \code{"!="}. \code{"abs >"} tests for magnitude by taking the absolute
#'  value and then testing if they are greater than... .
#' @param method_p If the input is a \code{"tICA.[format]"} model object, the type of
#'  multiple comparisons correction to use for p-values, or \code{NULL} for no
#'  correction. See \code{help(p.adjust)}. Default: \code{"BH"} (Benjamini &
#'  Hochberg, i.e. the false discovery rate). Note that multiple comparisons
#'  will account for data locations, but not ICs.
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
#' @importFrom grDevices colorRamp rgb
#'
#' @examples
#' \dontrun{
#'  activations(tICA_result, alpha=.05, deviation=TRUE)
#' }
activations <- function(
  tICA, u=NULL, z=NULL, alpha=0.01,
  type=c(">", "abs >", "<", "!="),
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
  type <- match.arg(type, c(">", "abs >", "<", "!="))
  if(alpha <= 0 | alpha >= 1) stop('alpha must be between 0 and 1')
  stopifnot(fMRItools::is_1(verbose, "logical"))
  stopifnot(fMRItools::is_1(deviation, "logical"))

  nL <- ncol(tICA$A)
  if(is.null(which.ICs)) which.ICs <- seq(nL)
  stopifnot(is.numeric(which.ICs))
  stopifnot(which.ICs == round(which.ICs))
  if(min((which.ICs) %in% (1:nL))==0) stop('Invalid entries in which.ICs')
  nI <- length(which.ICs)

  if (deviation) {
    if (!is.null(z)) { stop("`z` not compatible with `deviation==TRUE`.") }
    if (u != 0) { warning("`u != 0` not advised for `deviation==TRUE`. Proceeding anyway.") }
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

  nV <- nrow(tICA$s_mean)

  # Convert `z` to `u`.
  # Make `u` a nU x nL matrix (cutoffs by ICs).
  u_og <- u
  if (!is.null(z)) {
    stopifnot(!is.matrix(z))
    z <- sort(z, decreasing=type=="<")
    u_mat <- outer(z, sqrt(colVars(tICA$t_mean[,which.ICs,drop=FALSE])))
  } else if (!is.null(u)) {
    stopifnot(!is.matrix(u))
    u <- sort(u, decreasing=type=="<")
    u_mat <- outer(u, rep(1, nL))
  } else {
    u_mat <- matrix(0, nrow=1, ncol=nL)
  }
  nU <- nrow(u_mat)

  if (type == "!=" && nU > 1) {
    stop("Multiple u/z not compatible with '!=' test.")
  }

  if (type == "abs >") {
    tICA$s_mean <- abs(tICA$s_mean)
    type <- ">"
  }

  act_name <- format_activation_name(
    u=u, z=z, type=type, deviation=deviation, collapse=FALSE
  )

  # Loop over `u` to compute activations. --------------------------------------
  out <- vector("list", nU)
  names(out) <- act_name
  for (uu in seq(nU)) {
    if (verbose) { cat(act_name[uu], ".\n") }
    uu_mat <- matrix(u_mat[uu,], nrow=nV, ncol=nL, byrow=TRUE)

    # Spatial template ICA activations -----------------------------------------
    if (is_stICA) {
      if(verbose) cat('Determining areas of activations based on joint posterior distribution of latent fields\n')

      #identify areas of activation in each IC
      active <- jointPPM <- marginalPPM <- vars <- matrix(NA, nrow=nV, ncol=nL)

      for(q in which.ICs){
        if(verbose) cat(paste0('.. IC ',q,' (',which(which.ICs==q),' of ',length(which.ICs),') \n'))
        inds_q <- (1:nV) + (q-1)*nV
        if(deviation){
          Dinv_mu_s <-  (as.vector(tICA$s_mean) - as.vector(tICA$t_mean) - uu_mat)/as.vector(sqrt(tICA$t_var))
        } else {
          Dinv_mu_s <- (as.vector(tICA$s_mean) - uu_mat)/as.vector(sqrt(tICA$t_var))
        }

        if(q==which.ICs[1]) {
          #we scale mu by D^(-1) to use Omega for precision (actual precision of s|y is D^(-1) * Omega * D^(-1) )
          #we subtract u first since rescaling by D^(-1) would affect u too
          #save rho from first time running excursions, pass into excursions for other ICs
          tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = Q, type = type, u = u_mat[1], ind = inds_q)) #I think u does not matter, should check because may differ across fields
          if(verbose) print(tmp)
          rho <- res_q$rho
        } else {
          tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = Q, type = type, u = u_mat[1], ind = inds_q, rho=rho))
          if(verbose) print(tmp)
        }
        active[,q] <- res_q$E[inds_q]
        jointPPM[,q] <- res_q$F[inds_q]
        marginalPPM[,q] <- res_q$rho[inds_q]
        vars[,q] <- res_q$vars[inds_q]
      }

      if (length(unique(u))==1) { u <- u[1] }
      if (length(unique(z))==1) { z <- z[1] }
      out[[uu]] <- list(
        active=active, jointPPM=jointPPM, marginalPPM=marginalPPM, vars=vars
      )
    }

    # Template ICA activations -------------------------------------------------
    if (is_tICA) {
      if(verbose) cat('\tDetermining areas of activations based on hypothesis testing at each location\n')

      nL <- ncol(tICA$s_mean)
      if(deviation){
        t_stat <- (as.matrix(tICA$s_mean) - tICA$t_mean - uu_mat) / as.matrix(tICA$s_se)
      } else {
        t_stat <- (as.matrix(tICA$s_mean) - uu_mat) / as.matrix(tICA$s_se)
      }

      if(type=='>') pvals <- 1-pnorm(t_stat)
      if(type=='<') pvals <- pnorm(t_stat)
      if(type=='!=') pvals <- 2*(1-pnorm(abs(t_stat)))

      if(verbose) cat(paste0('\tCorrecting for multiple comparisons with method ', method_p, "\n"))

      pvals_adj <- active <- matrix(NA, nrow=nV, ncol=nL)
      if(is.null(method_p)) method_p <- 'none'
      for(q in which.ICs){
        pvals_adj[,q] <- p.adjust(pvals[,q], method=method_p)
        active[,q] <- (pvals_adj[,q] < alpha)
      }

      out[[uu]] <- list(
        active = active,
        pvals = pvals, pvals_adj = pvals_adj,
        se = tICA$s_se, tstats = t_stat
      )
    }
  }

  # Format result. -------------------------------------------------------------
  active <- rowSums(abind(lapply(out, "[[", "active"), along=3), dims=2)
  active[] <- as.numeric(active)
  dimnames(active) <- NULL

  # Unmask data.
  if (use_mask) { active <- fMRItools::unmask_mat(active, mask) }

  # Get colors.
  if (FORMAT == "CIFTI") {
    active_colors <- rev(ciftiTools::make_color_pal("plasma")$color)[round(seq(nU)/nU*256)]
    # ac_rgb <- grDevices::colorRamp(c("white", "red"), space="Lab")(seq(0, nU)/(nU))
    # active_colors <- vector("character", nU+1)
    # for (uu in seq(nU+1)) {
    #   active_colors[uu] <- rgb(ac_rgb[uu,1], ac_rgb[uu,2], ac_rgb[uu,3], 255, maxColorValue=255)
    # }
  }

  # Un-vectorize data.
  if (FORMAT == "CIFTI") {
    active <- ciftiTools::newdata_xifti(xii1, active)
    active <- ciftiTools::move_from_mwall(active, -1)
    active <- ciftiTools::convert_xifti(
      active, "dlabel",
      levels_old=c(-1, 0, seq(nU)),
      levels=c(-1, 0, seq(nU)),
      labels=c("Medial Wall", "Inactive", paste("Active:", act_name)),
      colors=c("#888888", "white", active_colors),
      add_white=FALSE
    )
  }

  result <- c(
    list(active=active),
    out,
    list(params=list(
      alpha=alpha, method_p=method_p, type=type, u=u, z=z, deviation=deviation
    ))
  )

  if (FORMAT == "CIFTI") {
    class(result) <- "tICA_act.cifti"
  } else if (FORMAT == "NIFTI") {
    active <- fMRItools::unvec_vol(active, mask_nii, fill=NA)
    class(result) <- "tICA_act.nifti"
  } else {
    class(result) <- "tICA_act.matrix"
  }

  result
}

