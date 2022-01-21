#' Template ICA
#'
#' Perform template independent component analysis (ICA) using expectation-maximization (EM)
#'
#' @param BOLD Vector of subject-level fMRI data in one of the following formats: 
#'  CIFTI file paths, \code{"xifti"} objects, NIFTI file paths, \code{"nifti"} objects, or
#'  \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data locations and
#'  \eqn{T} is the number of timepoints. 
#' 
#'  If multiple BOLD data are provided, they will be independently centered, scaled, and detrended (if applicable),
#'  and then they will be concatenated together prior to denoising (if applicable) and computing the initial dual regression estimate.
#' @param template_mean,template_var,template_FC Template estimates in a format
#'  compatible with \code{BOLD}. 
#' 
#'  \code{template_mean} gives the mean of the empirical population prior for each of the \eqn{L} independent components. 
#'  Its vectorized dimensions should be \eqn{V \times L}.
#'  
#'  \code{template_var} gives the between-subject variance of the empirical population prior for each of the \eqn{L} independent components.
#'  Its vectorized dimensions should be \eqn{V \times L}. If not provided and \code{template_mean}
#'  is a file path, will be inferred by substituting \code{"_mean"} in the file name with \code{"_var"}.
#' 
#'  \code{template_FC} is not yet supported. 
#' @param center_Bcols Center BOLD across columns (each image)? Default: \code{FALSE} (recommended).
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation).
#' @param detrend_DCT Detrend the data? This is the number of DCT bases to use for detrending. If \code{0} (default), do not detrend.
#' @param normA Scale each IC timeseries (column of \eqn{A}) in the dual regression 
#'  estimates? Default: \code{FALSE}. (The opposite scaling will be applied to \eqn{S}
#'  such that the product \eqn{A \times S} remains the same).
#' @param Q2,Q2_max Denoise the BOLD data? Denoising is based on modeling and 
#'  removing nuisance ICs. It may result in a cleaner estimate for smaller 
#'  datasets, but it may be unnecessary (and time-consuming) for larger datasets. 
#'  
#'  Set \code{Q2} to control denoising: use a positive integer to specify the
#'  number of nuisance ICs, \code{NULL} to have the number of nuisance ICs
#'  estimated by PESEL (default), or zero to skip denoising. 
#' 
#'  If \code{is.null(Q2)}, use \code{Q2_max} to specify the maximum number of
#'  nuisance ICs that should be estimated by PESEL. \code{Q2_max} must be less
#'  than \eqn{T * .75 - Q} where \eqn{T} is the number of timepoints in BOLD
#'  and \eqn{Q} is the number of group ICs. If \code{NULL} (default),
#'  \code{Q2_max} will be set to \eqn{T * .50 - Q}, rounded.
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI file paths. 
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI file paths or
#'  \code{"nifti"} objects. This is a brain map formatted as a binary array of the same 
#'  size as the fMRI data, with \code{TRUE} corresponding to in-mask voxels.
#' @param time_inds Subset of fMRI BOLD volumes to include in analysis. 
#'  If \code{NULL} (default), use all volumes. Subsetting is performed before
#'  any centering, scaling, detrending, and denoising.
#' @param spatial_model Should spatial modeling be performed? If \code{NULL}, assume
#'  spatial independence. Otherwise, provide meshes specifying the spatial priors assumed on
#'  each independent component. Each should represent a brain structure, between which
#'  spatial independence can be assumed.
#' 
#'  If \code{BOLD} represents CIFTI-format data, \code{spatial_model} should give the left and
#'  right cortex surface geometries (whichever one(s) are being used) as \code{"surf"}
#'  objects or GIFTI surface geometry file paths. Spatial modeling is not yet available for
#'  the subcortex. This argument can also be \code{TRUE}, in which case spatial modeling
#'  will be performed with the surfaces included in the first entry of \code{BOLD} if it is a 
#'  \code{"xifti"} object, or if those are not present available, the default inflated
#'  surfaces from \code{ciftiTools}.
#' 
#'  If \code{BOLD} represents NIFTI-format data, spatial modeling is not yet available.
#' 
#'  If \code{BOLD} is a numeric matrix, \code{spatial_model} should be a list of meshes
#'  (see \code{\link{make_mesh}}).
#' @param resamp_res Only applies if \code{BOLD} represents CIFTI-format data, and only recommended
#'  for spatial modeling. The target resolution for resampling (number of
#'  cortical surface vertices per hemisphere). A value less than 10000 is recommended for computational
#'  feasibility. If \code{NULL} (default), do not perform resampling.
#' @param rm_mwall Only applies if \code{BOLD} represents CIFTI-format data. Should medial wall (missing data) locations be removed from the mesh?  
#'  If \code{TRUE}, faster computation but less accurate estimates at the boundary of wall.
#' @param maxiter Maximum number of EM iterations. Default: \code{100}.
#' @param epsilon Smallest proportion change between iterations. Default: \code{.001}.
#' @param kappa_init Starting value for kappa.  Default: \code{0.2}.
#' @param usePar Parallelize the computation over data locations? Default: \code{FALSE}. Can be the number of cores
#'  to use or \code{TRUE}, which will use the number on the PC minus two.
#' @param verbose If \code{TRUE}. display progress of algorithm
# @param common_smoothness If \code{TRUE}. use the common smoothness version
#  of the spatial template ICA model, which assumes that all IC's have the same
#  smoothness parameter, \eqn{\kappa}
#'
#' @return A list containing the estimated independent components S
#'  (a VxL matrix), their mixing matrix A (a TxL matrix), and the number of
#'  nuisance ICs estimated (Q_nuis)
#'
#' @export
#'
# @importFrom INLA inla inla.spde.result inla.pardiso.check inla.setOption
#' @import ciftiTools
#' @importFrom stats optim
#' @importFrom matrixStats rowVars
#'
templateICA <- function(
  BOLD,
  template_mean, template_var=NULL, template_FC=NULL,
  scale=TRUE, detrend_DCT=0, 
  center_Bcols=FALSE, normA=FALSE,
  Q2=NULL, Q2_max=NULL,
  brainstructures=c("left","right"), mask=NULL, time_inds=NULL,
  spatial_model=NULL, resamp_res=NULL, rm_mwall=TRUE,
  maxiter=100,
  epsilon=0.001,
  kappa_init=0.2,
  #common_smoothness=TRUE,
  usePar=FALSE,
  verbose=TRUE){

  # Check arguments ------------------------------------------------------------
  
  # Simple argument checks.
  stopifnot(is.logical(scale) && length(scale)==1)
  if (isFALSE(detrend_DCT)) { detrend_DCT <- 0 }
  stopifnot(is.numeric(detrend_DCT) && length(detrend_DCT)==1)
  stopifnot(detrend_DCT >=0 && detrend_DCT==round(detrend_DCT))
  stopifnot(is.logical(center_Bcols) && length(center_Bcols)==1)
  stopifnot(is.logical(normA) && length(normA)==1)
  if (!is.null(Q2)) { stopifnot(Q2 >= 0) } # Q2_max checked later.
  if (!is.null(resamp_res)) {
    stopifnot(is.numeric(resamp_res) && length(resamp_res)==1)
    stopifnot(round(resamp_res) == resamp_res && resamp_res >= 0)
  }
  stopifnot(is.logical(rm_mwall) && length(rm_mwall)==1)
  stopifnot(is.numeric(maxiter) && length(maxiter)==1)
  stopifnot(round(maxiter) == maxiter && maxiter >= 0)
  stopifnot(is.numeric(epsilon) && length(epsilon)==1)
  if (!is.null(kappa_init)) {
    stopifnot(is.numeric(kappa_init) && length(kappa_init)==1)
    if(kappa_init <= 0) stop('kappa_init must be a positive scalar or `NULL`.')
  }
  stopifnot(is.logical(usePar) && length(usePar)==1)
  stopifnot(is.logical(verbose) && length(verbose)==1)
  xii1 <- NULL

  # `usePar`
  if (!isFALSE(usePar)) {
    parPkgs <- c("parallel", "doParallel")
    parPkgs_missing <- vapply(parPkgs, function(x){requireNamespace(x, quietly=TRUE)}, FALSE)
    if (any(parPkgs_missing)) {
      if (all(parPkgs_missing)) {
        stop(paste0(
          "Packages `parallel` and `doParallel` needed ",
          "for `usePar` to loop over voxels. Please install it.", call.=FALSE
        ))
      } else {
        stop(paste0(
          "Package `", parPkgs[parPkgs_missing], "` needed ",
          "for `usePar` to loop over voxels. Please install it.", call.=FALSE
        ))
      }
    }

    cores <- parallel::detectCores()
    if (isTRUE(usePar)) { nCores <- cores[1] - 2 } else { nCores <- usePar; usePar <- TRUE }
    if (nCores < 2) { 
      usePar <- FALSE
    } else {
      cluster <- parallel::makeCluster(nCores)
      doParallel::registerDoParallel(cluster)
    }
  }

  # `out_fname`?

  # `BOLD` ---------------------------------------------------------------------
  # Determine the format of `BOLD`. 
  format <- infer_BOLD_format(BOLD)
  FORMAT <- switch(format,
    CIFTI = "CIFTI",
    xifti = "CIFTI",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    data = "DATA"
  )
  FORMAT_extn <- switch(FORMAT, CIFTI=".dscalar.nii", NIFTI=".nii", DATA=".rds")

  # [TO DO] if FORMAT=="NIFTI", take intersection of template mask and provided mask.

  # If BOLD (and BOLD2) is a CIFTI or NIFTI file, check that the file paths exist.
  if (format %in% c("CIFTI", "NIFTI")) {
    missing_BOLD <- !file.exists(BOLD)
    if (all(missing_BOLD)) stop('The files in `BOLD` do not exist.')
    if (any(missing_BOLD)) {
      warning('There are ', missing_BOLD, ' scans in `BOLD` that do not exist. These scans will be excluded from template estimation.')
      BOLD <- BOLD[!missing_BOLD]
    }
  }

  # Make `BOLD` a list.
  if (is.character(BOLD)) { 
    BOLD <- as.list(BOLD)
  } else if (!is.list(BOLD)) {
    BOLD <- list(BOLD)
  }
  nN <- length(BOLD)

  # `brainstructures`
  if (FORMAT == "CIFTI") {
    brainstructures <- match.arg(brainstructures, c("left", "right", "subcortical", "all"), several.ok=TRUE)
    if ("all" %in% brainstructures) { brainstructures <- c("left", "right", "subcortical") }
    do_left <- "left" %in% brainstructures
    do_right <- "right" %in% brainstructures
    do_sub <- "subcortical" %in% brainstructures
  }

  # templates ------------------------------------------------------------------
  # Conver templates to numeric data matrices or arrays.
  # Check that the mean and variance template dimensions match.
  if (FORMAT == "CIFTI") {
    if (is.character(template_mean)) { 
      if (is.null(template_var)) { 
        template_var <- gsub(
          paste0("_mean", FORMAT_extn), paste0("_var", FORMAT_extn), 
          template_mean
        )
        if (!file.exists(template_var)) { 
          stop("Could not infer `template_var` file path; please provide it.")
        }
      }
      template_mean <- ciftiTools::read_xifti(template_mean, brainstructures=brainstructures) 
      template_var <- ciftiTools::read_xifti(template_var, brainstructures=brainstructures) 
    }
    if (is.xifti(template_mean)) { 
      xii1 <- select_xifti(template_mean, 1) # for formatting output
      template_mean <- as.matrix(template_mean)
    }
    if (is.xifti(template_var)) {
      if (is.null(xii1)) { xii1 <- select_xifti(template_var, 1) }
      template_var <- as.matrix(template_var)
    }
    stopifnot(is.matrix(template_mean))
    stopifnot(is.matrix(template_var))
  } else if (FORMAT == "NIFTI") {
    if (is.character(template_mean)) {
      if (is.null(template_var)) { 
        template_var <- gsub(
          paste0("_mean", FORMAT_extn), paste0("_var", FORMAT_extn), 
          template_mean
        )
        if (!file.exists(template_var)) { 
          stop("Could not infer `template_var` file path; please provide it.")
        }
      }      
      template_mean <- oro.nifti::readNIfTI(template_mean, reorient=FALSE)
      template_var <- oro.nifti::readNIfTI(template_var, reorient=FALSE)
    }
    stopifnot(length(dim(template_mean)) > 1)
    stopifnot(length(dim(template_var)) > 1)
  } else {
    stopifnot(is.matrix(template_mean))
    stopifnot(is.matrix(template_var))
  }
  stopifnot(length(dim(template_mean)) == length(dim(template_var)))
  stopifnot(all(dim(template_mean) == dim(template_var)))
  nL <- dim(template_mean)[length(dim(template_mean))]

  # `mask` ---------------------------------------------------------------------
  # Get `mask` as a logical array.
  # Check templates and `mask` dimensions match.
  # Vectorize templates.
  if (FORMAT == "NIFTI") {
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- oro.nifti::readNIfTI(mask, reorient=FALSE) }
    if (dim(mask)[length(dim(mask))] == 1) { mask <- array(mask, dim=dim(mask)[length(dim(mask))-1]) }
    if (is.numeric(mask)) { 
      cat("Coercing `mask` to a logical array with `as.logical`.\n")
      mask[] <- as.logical(mask)
    }
    nI <- dim(mask); nV <- sum(mask)
    stopifnot(length(dim(template_mean)) %in% c(2, length(nI)+1))
    if (length(dim(template_mean)) == length(nI)+1) {
      if (length(dim(template_mean)) != 2) {
        stopifnot(all(dim(template_mean)[length(dim(template_mean))-1] == nI))
      }
      if (all(dim(template_mean)[length(dim(template_mean))-1] == nI)) {
        template_mean <- matrix(template_mean[rep(mask, nL)], ncol=nL)
        template_var <- matrix(template_var[rep(mask, nL)], ncol=nL)
        stopifnot(nrow(template_mean) == nV)  
      }
    }
  } else {
    nI <- nV <- nrow(template_mean)
  }

  # `spatial_model` ------------------------------------------------------------
  do_spatial <- !is.null(spatial_model)
  if (isTRUE(spatial_model)) { spatial_model <- NULL }
  meshes <- spatial_model
  if (do_spatial) {
    # Check that `BOLD` format is compatible with the spatial model
    if (FORMAT == "NIFTI") { stop("`spatial_model` not available for NIFTI BOLD.") }
    if (FORMAT == "CIFTI") {
      if ("subcortical" %in% brainstructures) { 
        stop("Subcortical `brainstructures` not compatible with `spatial_model.`")
      }
    }

    # INLA
    INLA_check()
    flag <- INLA::inla.pardiso.check()
    if (any(grepl('FAILURE',flag))) {
      stop(
        'PARDISO IS NOT INSTALLED OR NOT WORKING. ',
        'PARDISO for R-INLA is required for computational efficiency. ',
        'If you already have a PARDISO / R-INLA License, run inla.setOption(pardiso.license = "/path/to/license") and try again. ',
        'If not, run inla.pardiso() to obtain a license.'
      )
    }
    INLA::inla.setOption(smtp='pardiso')

    if (verbose) {
      mesh_name <- ifelse(is.null(meshes), "the default inflated surface", "the provided mesh")
      cat(paste0(
        "Fitting a spatial model based on ", mesh_name, ". ",
        "Note that computation time and memory demands may be high.\n"
      ))
    }

    if (FORMAT == "CIFTI") {
      if (is.null(meshes)) {
        if (is.null(resamp_res)) {
          res <- ciftiTools:::infer_resolution(BOLD)
        } else {
          res <- resamp_res
        }
        if (do_left) {
          surf <- BOLD[[1]]$surf$cortex_left
          if (is.null(surf)) { surf <- make_surf(ciftiTools.files()$surf["left"], resamp_res=res) }
          if (!is.null(BOLD[[1]]$meta$cortex$medial_wall_mask$left)) {
            wall_mask <- which(BOLD[[1]]$meta$cortex$medial_wall_mask$left)
          } else {
            wall_mask <- NULL
          }
          if (rm_mwall) { 
            mesh <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
          } else {
            mesh <- make_mesh(surf = surf, inds_data = wall_mask) #retain wall in mesh for more accuracy along boundary with wall
          }
          meshes <- c(meshes, list(left=mesh))
        }
        if (do_right) {
          surf <- BOLD[[1]]$surf$cortex_right
          if (is.null(surf)) { surf <- make_surf(ciftiTools.files()$surf["right"], resamp_res=res) }
          if (!is.null(BOLD[[1]]$meta$cortex$medial_wall_mask$right)) {
            wall_mask <- which(BOLD[[1]]$meta$cortex$medial_wall_mask$right)
          } else {
            wall_mask <- NULL
          }
          if (rm_mwall) { 
            mesh <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
          } else {
            mesh <- make_mesh(surf = surf, inds_data = wall_mask) #retain wall in mesh for more accuracy along boundary with wall
          }
          meshes <- c(meshes, list(right=mesh))
        }
      }
    } else {
      if (is.null(meshes)) { 
        stop("`meshes` must be provided if the input format is not CIFTI.") 
      }
    }
    if (!is.list(meshes)) stop('`meshes` must be a list.')
    if (!all(vapply(meshes, inherits, what="templateICA_mesh", FALSE))) {
      stop('Each element of `meshes` should be of class `"templateICA_mesh"`. See `help(make_mesh)`.')
    }
    ndat_mesh <- sum(vapply(meshes, function(x){colSums(x$A)}, 0))
    if (ndat_mesh != nV) { 
      stop("Number of data locations in `meshes` does not match that of the BOLD data.") 
    }
    # [TO-DO]: Check that numbers of data locations on meshes (column sums of A) add up to match the number of data locations.
    #if(class(common_smoothness) != 'logical' | length(common_smoothness) != 1) stop('common_smoothness must be a logical value')
    #if(!do_spatial & !is.null(kappa_init)) stop('kappa_init should only be provided if mesh also provided for spatial modeling')
  }

  #TO DO: Define do_FC, pass in FC template
  do_FC <- FALSE

  # Process the scan -----------------------------------------------------------
  if (verbose) {
    cat('Number of data locations:      ', nV, "\n")
    cat('Number of template ICs:        ', nL, "\n")
    cat('Number of BOLD scans:          ', nN, "\n")
  }

  # Get each entry of `BOLD` as a data matrix or array. 
  if (format == "CIFTI") {
    for (bb in seq(nN)) {
      if (is.character(BOLD[[bb]])) { BOLD[[bb]] <- ciftiTools::read_xifti(BOLD[[bb]], brainstructures=brainstructures) }
      if (is.xifti(BOLD[[bb]])) { BOLD[[bb]] <- as.matrix(BOLD[[bb]]) }
      stopifnot(is.matrix(BOLD[[bb]]))
    }
  } else if (format == "NIFTI") {
    for (bb in seq(nN)) {
      if (is.character(BOLD[[bb]])) { BOLD[[bb]] <- oro.nifti::readNIfTI(BOLD[[bb]], reorient=FALSE) }
      stopifnot(length(dim(BOLD[[bb]])) > 1)
    }
  } else {
    for (bb in seq(nN)) {
      stopifnot(is.matrix(BOLD[[bb]]))
    }
  }

  # If `BOLD` is a list, ensure that all dimensions are the same.
  dBOLDs <- lapply(BOLD, dim)
  if (length(unique(vapply(dBOLDs, length, 0))) > 1) {
    stop("`BOLD` do not have the same dimensions.")
  }
  dBOLD <- dBOLDs[[1]]
  ldB <- length(dBOLD)
  if (nN > 1) {
    for (bb in seq(2, nN)) {
      stopifnot(all(dBOLD[seq(ldB-1)] == dBOLDs[[bb]][seq(ldB-1)]))
    }
  }
  nT <- vapply(dBOLDs, function(x){x[ldB]}, 0)

  # `time_inds`
  if (!is.null(time_inds)) {
    for (bb in seq(nN)) {
      if (!all(time_inds %in% seq(nT))) {
        warning('Not all `time_inds` available.') # [TO DO]: improve
        time_inds_bb <- intersect(time_inds, seq(nT))
      } else {
        time_inds_bb <- time_inds
      }
      BOLD[[bb]] <- BOLD[[bb]][,time_inds_bb,drop=FALSE]
    }
  }
  nT <- vapply(dBOLDs, function(x){x[ldB]}, 0)
  nTmin <- min(nT)

  # Check `BOLD` dimensions correspond with template and `mask`.
  stopifnot(ldB-1 == length(nI))
  stopifnot(all(dBOLD[seq(ldB-1)] == nI))

  # Vectorize `BOLD`.
  if (FORMAT=="NIFTI") {
    for (bb in seq(nN)) {
      BOLD[[bb]] <- matrix(BOLD[[bb]][rep(mask, dBOLD[ldB])], ncol=nT)
      stopifnot(nrow(BOLD[[bb]]) != nV)
    }
  }

  # Check that numbers of data locations (nV), time points (nT) and ICs (nL) makes sense, relatively
  if (sum(nT) > nV) warning('More time points than voxels. Are you sure?')
  if (nL > nV) stop('The arguments you supplied suggest that you want to estimate more ICs than you have data locations.  Please check the orientation and size of template_mean, template_var and BOLD.')
  if (nL > sum(nT)) stop('The arguments you supplied suggest that you want to estimate more ICs than you have time points.  Please check the orientation and size of template_mean, template_var and BOLD.')

  ### IDENTIFY AND REMOVE ANY BAD VOXELS/VERTICES
  keep <- rep(TRUE, nV)
  keep[rowSums(is.nan(template_mean)) > 0] <- FALSE
  keep[rowSums(is.na(template_mean)) > 0] <- FALSE
  keep[rowSums(is.nan(template_var)) > 0] <- FALSE
  keep[rowSums(is.na(template_var)) > 0] <- FALSE
  for (bb in seq(nN)) {
    keep[rowSums(is.nan(BOLD[[bb]])) > 0] <- FALSE
    keep[rowSums(is.na(BOLD[[bb]])) > 0] <- FALSE
    keep[rowVars(BOLD[[bb]]) == 0] <- FALSE
  }
  if(sum(!keep) > 0){
    stop('flat or NA voxels detected in data or templates')
    # For this part, would need to also update "A" matrix (projection from mesh to data locations)
    # template_mean_orig <- template_mean
    # template_var_orig <- template_var
    # nV <- sum(keep)
    # if(verbose) cat(paste0('Excluding ',sum(!keep),' bad (NA, NaN or flat) voxels/vertices from analysis.\n'))
    # template_mean <- template_mean[keep,]
    # template_var <- template_var[keep,]
    # BOLD <- BOLD[keep,]
  }

  # Normalize BOLD -------------------------------------------------------------
  # (Center, scale, and detrend)
  
  BOLD <- lapply(BOLD, norm_BOLD, 
    center_rows=TRUE, center_cols=center_Bcols, 
    scale=scale, detrend_DCT=detrend_DCT
  )

  # Concatenate the data.
  BOLD <- do.call(cbind, BOLD)
  nT <- sum(nT)

  # Estimate and deal with nuisance ICs ----------------------------------------
  if (is.null(Q2) || Q2>0) {
    BOLD <- rm_nuisIC(BOLD, template_mean=template_mean, Q2=Q2, Q2_max=Q2_max, verbose=verbose)
  }

  # Center and scale `BOLD` again, but do not detrend again. -------------------
  BOLD <- norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, detrend_DCT=FALSE
  )

  # Initialize with the dual regression-based estimate -------------------------
  BOLD_DR <- dual_reg(
    BOLD, template_mean, center_Bcols=FALSE, 
    scale=FALSE, detrend_DCT=0, normA=normA
  )

  # ----------------------------------------------------------------------------
  # EM -------------------------------------------------------------------------
  #Three algorithms to choose from:
  #1) FC Template ICA (new)
  #2) Template ICA
  #3) Spatial Template ICA (initialize with standard Template ICA)
  # ----------------------------------------------------------------------------

  #1) FC Template ICA ----------------------------------------------------------
  if(do_FC) {

    ## TO DO: FILL IN HERE
    prior_params <- c(0.001, 0.001) #alpha, beta (uninformative) -- note that beta is scale parameter in IG but rate parameter in the Gamma
    template_FC <- NULL
    EM_FCtemplateICA <- function(...){NULL}
    resultEM <- EM_FCtemplateICA(
      template_mean, template_var, template_FC,
      prior_params, #for prior on tau^2
      BOLD=BOLD,
      AS_init = BOLD_DR, #initial values for A and S
      maxiter=maxiter, epsilon=epsilon,
      verbose=verbose
    )
  } else {

    # Reduce data dimensions
    BOLD <- dim_reduce(BOLD, nL)
    err_var <- BOLD$sigma_sq
    BOLD2 <- BOLD$data_reduced
    H <- BOLD$H
    Hinv <- BOLD$H_inv
    # In original template ICA model nu^2 is separate
    #   for spatial template ICA it is part of C
    C_diag <- BOLD$C_diag
    if (do_spatial) { C_diag <- C_diag * (BOLD$sigma_sq) } #(nu^2)HH' in paper
    rm(BOLD)
    # Apply dimension reduction
    HA <- H %*% BOLD_DR$A 
    theta0 <- list(A = HA)
    # #initialize residual variance --- no longer do this, because we use dimension reduction-based estimate
    # theta0$nu0_sq = dat_list$sigma_sq
    # if(verbose) paste0('nu0_sq = ',round(theta0$nu0_sq,1)))

    #2) Template ICA -----------------------------------------------------------
    if (do_spatial && verbose) { cat('INITIATING WITH STANDARD TEMPLATE ICA\n') }
    theta00 <- theta0
    theta00$nu0_sq <- err_var
    resultEM <- EM_templateICA.independent(template_mean,
                                           template_var,
                                           BOLD=BOLD2,
                                           theta0=theta00,
                                           C_diag=C_diag,
                                           maxiter=maxiter,
                                           epsilon=epsilon,
                                           verbose=verbose)
    resultEM$A <- Hinv %*% resultEM$theta_MLE$A
    class(resultEM) <- 'tICA'

    #3) Spatial Template ICA ---------------------------------------------------
    if(do_spatial){
      resultEM_tICA <- resultEM
      theta0$kappa <- rep(kappa_init, nL)
      if(verbose) cat('ESTIMATING SPATIAL MODEL\n')
      t000 <- Sys.time()
      resultEM <- EM_templateICA.spatial(template_mean,
                                         template_var,
                                         meshes,
                                         BOLD=BOLD2,
                                         theta0,
                                         C_diag,
                                         maxiter=maxiter,
                                         epsilon=epsilon,
                                         verbose=verbose)
      #common_smoothness=common_smoothness)
      print(Sys.time() - t000)

      #organize estimates and variances in matrix form
      resultEM$subjICmean <- matrix(resultEM$subjICmean, ncol=nL)
      resultEM$subjICvar <- matrix(diag(resultEM$subjICcov), ncol=nL)
      resultEM$A <- Hinv %*% resultEM$theta_MLE$A

      resultEM$result_tICA <- resultEM_tICA
      class(resultEM) <- 'stICA'    
    }
  }

  # Return DR estimates.
  resultEM$result_DR <- BOLD_DR

  #This part problematic for spatial template ICA, but can bring back
  #for template ICA and FC template ICA.  When we check for bad locations,
  #can return an error only for spatial template ICA.

  #resultEM$keep <- keep
  # #map estimates & templates back to original locations
  # if(sum(!keep)>0){
  #   #estimates
  #   subjICmean <- subjICse <- matrix(nrow=length(keep), ncol=L)
  #   subjICmean[keep,] <- resultEM$subjICmean
  #   subjICse[keep,] <- resultEM$subjICse
  #   resultEM$subjICmean <- subjICmean
  #   resultEM$subjICse <- subjICse
  #   #templates
  #   resultEM$template_mean <- template_mean_orig
  #   resultEM$template_var <- template_var_orig
  # }

  # Format output.
  if (FORMAT=="CIFTI" && !is.null(xii1)) {
    resultEM$subjICmean <- newdata_xifti(xii1, resultEM$subjICmean)
    resultEM$subjICse <- newdata_xifti(xii1, resultEM$subjICse)
    if (do_spatial) {
      resultEM$result_tICA$subjICmean <- newdata_xifti(xii1, resultEM$result_tICA$subjICmean)
      resultEM$result_tICA$subjICse <- newdata_xifti(xii1, resultEM$result_tICA$subjICse)
    }
    class(resultEM) <- 'templateICA.cifti'

  } else if (FORMAT == "NIFTI") {
    subjICmean_nifti <- subjICvar_nifti <- template_mean
    newdata_nii <- function(dat, mask, nii_temp) {
      nL <- ncol(dat)
      for (qq in seq(nL)) {
        z <- mask
        z[mask] <- dat[,qq]
        nii_temp@.Data[,,,qq] <- z
      }
      nii_temp
    }
    resultEM$subjICmean <- newdata_nii(resultEM$subjICmean, mask, template_mean)
    resultEM$subjICvar <- newdata_nii(resultEM$subjICvar, mask, template_mean)
    if (spatial_model) {
      resultEM$subjICmean <- newdata_nii(resultEM$subjICmean, mask, template_mean)
      resultEM$subjICvar <- newdata_nii(resultEM$subjICvar, mask, template_mean)
    }
    resultEM$mask <- mask
    class(resultEM) <- 'templateICA.nifti'
  }

  resultEM
}
