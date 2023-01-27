#' Template ICA
#'
#' Perform template independent component analysis (ICA) using
#'  expectation-maximization (EM).
#'
#' @param BOLD Vector of subject-level fMRI data in one of the following
#'  formats: CIFTI file paths, \code{"xifti"} objects, NIFTI file paths,
#'  \code{"nifti"} objects, or \eqn{V \times T} numeric matrices, where \eqn{V}
#'  is the number of data locations and \eqn{T} is the number of timepoints.
#'
#'  If multiple BOLD data are provided, they will be independently centered,
#'  scaled, and detrended (if applicable), and then they will be concatenated
#'  together followed by denoising (if applicable) and computing the initial
#'  dual regression estimate.
#' @param template Template estimates in a format compatible with \code{BOLD},
#'  from \code{\link{estimate_template}}.
#' @param tvar_method Which calculation of the template variance to use:
#'  \code{"non-negative"} (default) or \code{"unbiased"}. The unbiased template
#'  variance is based on the assumed mixed effects/ANOVA model, whereas the
#'  non-negative template variance adds to it to account for greater potential
#'  between-subjects variation. (The template mean is the same for either choice
#'  of \code{tvar_method}.)
#' @inheritParams center_Bcols_Param
#' @inheritParams scale_Param
#' @param scale_sm_surfL,scale_sm_surfR,scale_sm_FWHM Only applies if
#'  \code{scale=="local"} and \code{BOLD} represents CIFTI-format data. To
#'  smooth the standard deviation estimates used for local scaling, provide the
#'  surface geometries along which to smooth as GIFTI geometry files or
#'  \code{"surf"} objects, as well as the smoothing FWHM (default: \code{2}).
#'
#'  If \code{scale_sm_FWHM==0}, no smoothing of the local standard deviation
#'  estimates will be performed.
#'
#'  If \code{scale_sm_FWHM>0} but \code{scale_sm_surfL} and
#'  \code{scale_sm_surfR} are not provided, the default inflated surfaces from
#'  the HCP will be used.
#'
#'  To create a \code{"surf"} object from data, see
#'  \code{\link[ciftiTools]{make_surf}}. The surfaces must be in the same
#'  resolution as the \code{BOLD} data.
#' @inheritParams detrend_DCT_Param
#' @inheritParams normA_Param
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
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI
#'  file paths. This is a character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI
#'  file paths or \code{"nifti"} objects. This is a brain map formatted as a
#'  binary array of the same spatial dimensions as the fMRI data, with
#'  \code{TRUE} corresponding to in-mask voxels.
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
#' @inheritParams varTol_Param
#' @param resamp_res Only applies if \code{BOLD} represents CIFTI-format data.
#'  The target resolution for resampling (number of cortical surface vertices
#'  per hemisphere). For spatial modelling, a value less than 10000 is
#'  recommended for computational feasibility. If \code{NULL} (default), do not
#'  perform resampling.
#' @param rm_mwall Only applies if \code{BOLD} represents CIFTI-format data.
#'  Should medial wall (missing data) locations be removed from the mesh?
#'  If \code{TRUE}, faster computation but less accurate estimates at the
#'  boundary of wall.
#' @param reduce_dim Reduce the temporal dimension of the data using PCA?
#'  Default: \code{TRUE}. Skipping dimension reduction will slow the model
#'  estimation, but may result in more accurate results.
#' @param maxiter Maximum number of EM iterations. Default: \code{100}.
#' @param epsilon Smallest proportion change between iterations. Default: \code{.01}.
#' @param kappa_init Starting value for kappa. Default: \code{0.2}.
#' @param usePar Parallelize the computation over data locations? Default:
#'  \code{FALSE}. Can be the number of cores to use or \code{TRUE}, which will
#'  use the number on the PC minus two.
#' @param verbose If \code{TRUE}, display progress of algorithm
# @param common_smoothness If \code{TRUE}. use the common smoothness version
#  of the spatial template ICA model, which assumes that all IC's have the same
#  smoothness parameter, \eqn{\kappa}
#'
#' @return A (spatial) template ICA object, which is a list containing:
#'  \code{subjICmean}, the \eqn{V \times L} estimated independent components
#'  \strong{S}; \code{subjICse}, the standard errors of \strong{S}; the
#'  \code{mask} of locations without template values due to too many low
#'  variance or missing values; and the function \code{params} such as
#'  the type of scaling and detrending performed.
#'
#'  If \code{BOLD} represented CIFTI or NIFTI data, \code{subjICmean} and
#'  \code{subjICse} will be formatted as \code{"xifti"} or \code{"nifti"}
#'  objects, respectively.
#'
#' @export
#'
# @importFrom INLA inla inla.spde.result inla.pardiso.check inla.setOption
#' @importFrom fMRItools infer_format_ifti unmask_mat unvec_vol dim_reduce
#' @importFrom stats optim
#' @importFrom matrixStats rowVars
#'
#' @examples
#' \dontrun{
#'  tm <- estimate_template(cii1_fnames, cii2_fnames, gICA_fname)
#'  templateICA(newcii_fname, tm, spatial_model=TRUE, resamp_res=2000)
#' }
templateICA <- function(
  BOLD, template, tvar_method=c("non-negative", "unbiased"),
  scale=c("global", "local", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  detrend_DCT=0,
  center_Bcols=FALSE,
  normA=FALSE,
  Q2=NULL,
  Q2_max=NULL,
  brainstructures=c("left","right"), mask=NULL, time_inds=NULL,
  varTol=1e-6,
  spatial_model=NULL, resamp_res=NULL, rm_mwall=TRUE,
  reduce_dim=TRUE,
  maxiter=100,
  epsilon=0.01,
  kappa_init=0.2,
  #common_smoothness=TRUE,
  usePar=FALSE,
  verbose=TRUE){

  # Check arguments ------------------------------------------------------------

  # Simple argument checks.
  tvar_method <- match.arg(tvar_method, c("non-negative", "unbiased"))
  tvar_name <- switch(tvar_method, `non-negative`="varNN", unbiased="varUB")
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("global", "local", "none"))
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  if (isFALSE(detrend_DCT)) { detrend_DCT <- 0 }
  stopifnot(is.numeric(detrend_DCT) && length(detrend_DCT)==1)
  stopifnot(detrend_DCT >=0 && detrend_DCT==round(detrend_DCT))
  stopifnot(is.logical(center_Bcols) && length(center_Bcols)==1)
  stopifnot(is.logical(normA) && length(normA)==1)
  if (!is.null(Q2)) { stopifnot(Q2 >= 0) } # Q2_max checked later.
  stopifnot(is.numeric(varTol) && length(varTol)==1)
  if (varTol < 0) { cat("Setting `varTol=0`."); varTol <- 0 }
  if (isFALSE(spatial_model)) { spatial_model <- NULL }
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
  stopifnot((is.logical(usePar)||is.numeric(usePar)) && length(usePar)==1)
  stopifnot(is.logical(verbose) && length(verbose)==1)

  # `usePar`
  if (!isFALSE(usePar)) {
    parPkgs <- c("parallel", "doParallel")
    parPkgs_missing <- !vapply(parPkgs, function(x){requireNamespace(x, quietly=TRUE)}, FALSE)
    if (any(parPkgs_missing)) {
      if (all(parPkgs_missing)) {
        stop(paste0(
          "Packages `parallel` and `doParallel` needed ",
          "for `usePar` to loop over voxels. Please install them."), call.=FALSE
        )
      } else {
        stop(paste0(
          "Package `", parPkgs[parPkgs_missing], "` needed ",
          "for `usePar` to loop over voxels. Please install it."), call.=FALSE
        )
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
  format <- fMRItools::infer_format_ifti(BOLD)[1]
  FORMAT <- get_FORMAT(format)
  FORMAT_extn <- switch(FORMAT, CIFTI=".dscalar.nii", GIFTI=".func.gii", NIFTI=".nii", MATRIX=".rds")

  check_req_ifti_pkg(FORMAT)

  # If BOLD (and BOLD2) is a CIFTI, GIFTI, NIFTI, or RDS file, check that the file paths exist.
  if (format %in% c("CIFTI", "GIFTI", "NIFTI", "RDS")) {
    missing_BOLD <- !file.exists(BOLD)
    if (all(missing_BOLD)) stop('The files in `BOLD` do not exist.')
    if (any(missing_BOLD)) {
      warning(
        'There are ', missing_BOLD,
        ' scans in `BOLD` that do not exist. ',
        'These scans will be excluded from template estimation.'
      )
      BOLD <- BOLD[!missing_BOLD]
    }
  }

  # Make `BOLD` a list.
  if (is.character(BOLD)) {
    BOLD <- as.list(BOLD)
  } else if (!is.list(BOLD)) {
    BOLD <- list(BOLD)
  } else if (format == "xifti" && inherits(BOLD, "xifti")) {
    BOLD <- list(BOLD)
  } else if (format == "gifti" && inherits(BOLD, "gifti")) {
    BOLD <- list(BOLD)
  } else if (format == "nifti" && inherits(BOLD, "nifti")) {
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

  # Read in CIFTI, GIFTI, or NIFTI files.
  # (Need to do now rather than later, so that CIFTI resolution info can be used.)
  if (format == "CIFTI") {
    for (bb in seq(nN)) {
      if (is.character(BOLD[[bb]])) {
        BOLD[[bb]] <- ciftiTools::read_cifti(
          BOLD[[bb]], resamp_res=resamp_res,
          brainstructures=brainstructures
        )
      }
      stopifnot(ciftiTools::is.xifti(BOLD[[bb]]))
    }
  } else if (format == "GIFTI") {
    for (bb in seq(nN)) {
      if (is.character(BOLD[[bb]])) {
        BOLD[[bb]] <- gifti::readgii(BOLD[[bb]])
      }
      stopifnot(gifti::is.gifti(BOLD[[bb]]))
      if (bb == 1) {
        ghemi <- BOLD[[bb]]$file_meta["AnatomicalStructurePrimary"]
        if (!(ghemi %in% c("CortexLeft", "CortexRight"))) {
          stop("AnatomicalStructurePrimary metadata missing or invalid for BOLD #", bb, ".")
        }
      } else {
        if (ghemi != BOLD[[bb]]$file_meta["AnatomicalStructurePrimary"]) {
          stop("AnatomicalStructurePrimary metadata missing or invalid for BOLD #", bb, ".")
        }
      }
    }
    ghemi <- switch(ghemi, CortexLeft="left", CortexRight="right")
  } else if (format == "NIFTI") {
    for (bb in seq(nN)) {
      if (is.character(BOLD[[bb]])) {
        BOLD[[bb]] <- RNifti::readNifti(BOLD[[bb]])
      }
      # [TO DO] check?
    }
  } else if (format == "RDS") {
    for (bb in seq(nN)) {
      if (is.character(BOLD[[bb]])) {
        BOLD[[bb]] <- readRDS(BOLD[[bb]])
      }
    }
  }

  # templates ------------------------------------------------------------------
  if (is.character(template)) { template <- readRDS(template) }

  # Check template format matches BOLD format.
  TFORMAT <- class(template)[grepl("template", class(template))]
  if (length(TFORMAT) != 1) { stop("`template` is not a template.") }
  TFORMAT <- toupper(gsub("template.", "", TFORMAT, fixed=TRUE))
  if (TFORMAT != FORMAT) {
    stop("The BOLD format is '", FORMAT, ",' but the template format is '", TFORMAT, ".'")
  }

  # Check that parameters match.
  pmatch <- c(
    scale=scale, #scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT,
    center_Bcols=center_Bcols, normA=normA,
    # Q2=Q2, Q2_max=Q2_max,
    varTol=varTol
  )
  for (pp in seq(length(pmatch))) {
    pname <- names(pmatch)[pp]
    if (pmatch[pname] != template$params[[pname]]) {
      warning(paste0(
        "The `", pname, "` parameter was ",
        template$params[[pname]], " for the template, but is ",
        pmatch[pname], " here. These should match. (Proceeding anyway.)\n"
      ))
    }
  }

  #check for FC template
  do_FC <- FALSE; template_FC <- NULL
  if('FC' %in% names(template$template)) {
    do_FC <- TRUE
    template_FC <- template$template$FC
  }

  # Get `nI`, `nL`, and `nV`.
  # Check brainstructures and resamp_res if CIFTI, where applicable.
  # Convert to CIFTI if GIFTI.
  # Check mask if NIFTI.
  xii1 <- NULL
  if (FORMAT == "CIFTI") {
    xii1 <- template$dat_struct
    # Check brainstructures.
    tbs <- names(xii1$data)[!vapply(xii1$data, is.null, FALSE)]
    bs2 <- c(left="cortex_left", right="cortex_right", subcortical="subcort")[brainstructures]
    if (!all(bs2 %in% tbs)) {
      bs_missing <- bs2[!(bs2 %in% tbs)]
      stop(paste0(
        ifelse(length(bs_missing) > 1, "These brain structures are", "This brain structure is"),
        " not included in the template: ",
        paste(bs_missing, collapse=", "), ". Adjust the `brainstructures` argument accordingly."
      ))
    } else if (!all(tbs %in% bs2)) {
      bs_missing <- tbs[!(tbs %in% bs2)]
      if (verbose) {
        cat(paste0(
          "Ignoring ", ifelse(length(bs_missing) > 1, "these brain structures", "This brain structure"),
          " in the template:", paste(bs_missing, collapse=", "), "."
        ))
      }
      template <- removebs_template(template, bs_missing)
    }
    # Check `resamp_res`.
    if (!is.null(resamp_res)) {
      template <- resample_template(template, resamp_res=resamp_res)
    }
    nI <- nrow(template$template$mean)

  } else if (FORMAT == "GIFTI") {
    if (ghemi == "left") {
      xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexL=do.call(cbind, BOLD[[1]]$data)), 1) * 0
      for (bb in seq(length(BOLD))) {
        BOLD[[bb]] <- ciftiTools::as.xifti(cortexL=do.call(cbind, BOLD[[bb]]$data))
      }
    } else if (ghemi == "right") {
      xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexR=do.call(cbind, BOLD[[1]]$data)), 1) * 0
      for (bb in seq(length(BOLD))) {
        BOLD[[bb]] <- ciftiTools::as.xifti(cortexR=do.call(cbind, BOLD[[bb]]$data))
      }
    } else { stop() }
    # xii1 <- move_to_mwall(xii1)
    nI <- nrow(template$template$mean)

  } else if (FORMAT == "NIFTI") {
    # Check `mask`
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- RNifti::readNifti(mask); mask <- array(as.logical(mask), dim=dim(mask)) }
    if (dim(mask)[length(dim(mask))] == 1) { mask <- array(mask, dim=dim(mask)[length(dim(mask))-1]) }
    if (is.numeric(mask)) {
      cat("Coercing `mask` to a logical array.\n")
      mask <- array(as.logical(mask), dim=dim(mask))
    }
    nI <- dim(mask)
    # Check its compatibility with the template
    tds_dim <- dim(template$dat_struct) # has an empty last dim
    if (length(tds_dim) == length(nI) + 1 && tds_dim[length(tds_dim)] == 1) {
      tds_dim <- tds_dim[seq(length(tds_dim)-1)]
    }
    stopifnot(length(tds_dim) == length(nI))
    stopifnot(all(tds_dim == nI))
  } else {
    nI <- nrow(template$template$mean)
  }

  # Get IC inds.
  IC_inds <- template$params$inds

  # Get the data mask based on missing values, & low variance locations
  mask2 <- template$mask
  use_mask2 <- (!is.null(mask2)) && (!all(mask2))
  # Get the selected variance.
  template <- list(
    mean = template$template$mean,
    var = template$template[[tvar_name]]
  )
  # Make variance values non-negative (for unbiased template.)
  template$var[] <- pmax(0, template$var)
  # Get the template dimensions.
  nV <- nrow(template$mean)
  nL <- ncol(template$mean)

  # `spatial_model` meshes -----------------------------------------------------
  do_spatial <- !is.null(spatial_model)
  if (isTRUE(spatial_model)) { spatial_model <- NULL }
  meshes <- spatial_model
  if (do_spatial) {
    # Check that `BOLD` format is compatible with the spatial model
    if (FORMAT == "NIFTI") { stop("`spatial_model` not available for NIFTI BOLD.") }
    if (FORMAT == "CIFTI") {
      if ("subcortical" %in% brainstructures) {
        stop("Subcortex not compatible with `spatial_model.`")
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

    if (FORMAT %in% c("GIFTI", "CIFTI")) {
      if (is.null(meshes)) {
        if (is.null(resamp_res)) {
          res <- ciftiTools::infer_resolution(BOLD[[1]])
        } else {
          res <- resamp_res
        }
        if (do_left) {
          surf <- BOLD[[1]]$surf$cortex_left
          if (is.null(surf)) { surf <- ciftiTools::read_surf(ciftiTools::ciftiTools.files()$surf["left"], resamp_res=res) }
          if (!is.null(BOLD[[1]]$meta$cortex$medial_wall_mask$left)) {
            wall_mask <- BOLD[[1]]$meta$cortex$medial_wall_mask$left
            if (length(wall_mask) != nrow(surf$vertices)) { stop("Could not make surface of compatible resolution with data.") }
            wall_mask <- which(wall_mask)
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
          if (is.null(surf)) { surf <- ciftiTools::read_surf(ciftiTools::ciftiTools.files()$surf["right"], resamp_res=res) }
          if (!is.null(BOLD[[1]]$meta$cortex$medial_wall_mask$right)) {
            wall_mask <- BOLD[[1]]$meta$cortex$medial_wall_mask$right
            if (length(wall_mask) != nrow(surf$vertices)) { stop("Could not make surface of compatible resolution with data.") }
            wall_mask <- which(wall_mask)
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
    ndat_mesh <- sum(vapply(meshes, function(x){sum(x$A)}, 0))
    if (ndat_mesh != nV) {
      stop(
        "Total number of data locations in `meshes` (", ndat_mesh,
        ") does not match that of the templates (", nV, ")."
      )
    }
    # [TO-DO]: Check that numbers of data locations on meshes (column sums of A) add up to match the number of data locations.
    #if(class(common_smoothness) != 'logical' | length(common_smoothness) != 1) stop('common_smoothness must be a logical value')
    #if(!do_spatial & !is.null(kappa_init)) stop('kappa_init should only be provided if mesh also provided for spatial modeling')
  }


  # Process the scan -----------------------------------------------------------
  # Get each entry of `BOLD` as a data matrix or array.
  if (FORMAT %in% c("CIFTI", "GIFTI")) {
    for (bb in seq(nN)) {
      if (ciftiTools::is.xifti(BOLD[[bb]])) { BOLD[[bb]] <- as.matrix(BOLD[[bb]]) }
      stopifnot(is.matrix(BOLD[[bb]]))
    }
  } else if (FORMAT == "NIFTI") {
    for (bb in seq(nN)) {
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
  dBOLDs <- lapply(BOLD, dim)
  nT <- vapply(dBOLDs, function(x){x[ldB]}, 0)
  nTmin <- min(nT)

  if (verbose) {
    cat("Data input format:             ", format, "\n")
    cat('Number of data locations:      ', nV, "\n")
    if (FORMAT == "NIFTI") {
      cat("Unmasked dimensions:           ", paste(nI, collapse=" x "), "\n")
    }
    cat('Number of template ICs:        ', nL, "\n")
    cat('Number of BOLD scans:          ', nN, "\n")
    cat('Total number of timepoints:    ', sum(nT), "\n")
    cat('\n')
  }

  # Check `BOLD` dimensions correspond with template and `mask`.
  stopifnot(ldB-1 == length(nI))
  stopifnot(all(dBOLD[seq(ldB-1)] == nI))

  # Vectorize `BOLD` and apply `mask2`.
  for (bb in seq(nN)) {
    if (FORMAT == "NIFTI") {
      BOLD[[bb]] <- matrix(BOLD[[bb]][rep(mask, dBOLD[ldB])], ncol=nT)
      stopifnot(nrow(BOLD[[bb]]) == nV)
    }
    if (use_mask2) { BOLD[[bb]] <- BOLD[[bb]][mask2,,drop=FALSE] }
  }
  if (use_mask2) {
    template$mean <- template$mean[mask2,,drop=FALSE]
    template$var <- template$var[mask2,,drop=FALSE]
  }

  # Check that numbers of data locations (nV), time points (nT) and ICs (nL) makes sense, relatively
  if (sum(nT) > nV) warning('More time points than voxels. Are you sure?')
  if (nL > nV) stop('The arguments you supplied suggest that you want to estimate more ICs than you have data locations.  Please check the orientation and size of `BOLD` and `template`.')
  if (nL > sum(nT)) stop('The arguments you supplied suggest that you want to estimate more ICs than you have time points.  Please check the orientation and size of `BOLD` and `template`.')

  # Check that no NA or NaN values remain in template after masking.
  if (any(is.na(template$mean))) { stop("`NA` values in template mean.") }
  if (any(is.nan(template$mean))) { stop("`NaN` values in template mean.") }
  if (any(is.na(template$var))) { stop("`NA` values in template var.") }
  if (any(is.nan(template$mean))) { stop("`NaN` values in template mean.") }

  # Mask out additional locations due to data mask.
  mask3 <- apply(do.call(rbind, lapply(BOLD, make_mask, varTol=varTol)), 2, all)

  if (any(!mask3)) {
    if (do_spatial) {
      stop('Not supported yet: flat or NA voxels in data, after applying template mask, with spatial model.')
    }

    # [NOTE] For same results, would have needed to also update "A" matrix
    #   (projection from mesh to data locations)

    if(verbose) {
      cat(paste0(
        'Excluding ', sum(!mask3),
        ' bad (flat or NA) voxels/vertices from analysis.\n'
      ))
    }
    template_orig <- template
    BOLD <- lapply(BOLD, function(x){x[mask3,]})
    dBOLDs <- lapply(BOLD, dim); dBOLD <- dBOLDs[[1]]
    template$mean <- template$mean[mask3,]
    template$var <- template$var[mask3,]
    mask2[mask2][!mask3] <- FALSE
  }

  # Normalize BOLD -------------------------------------------------------------
  # (Center, scale, and detrend)
  if (verbose) { cat("Pre-processing BOLD data.\n") }

  if (!is.null(xii1) && scale=="local" && scale_sm_FWHM > 0) {
    xii1 <- ciftiTools::add_surf(xii1, surfL=scale_sm_surfL, surfR=scale_sm_surfR)
  }

  BOLD <- lapply(BOLD, norm_BOLD,
    center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT
  )

  # Concatenate the data.
  BOLD <- do.call(cbind, BOLD)
  nT <- sum(nT)

  # Estimate and deal with nuisance ICs ----------------------------------------
  x <- rm_nuisIC(BOLD, template_mean=template$mean, Q2=Q2, Q2_max=Q2_max, verbose=verbose, return_Q2=TRUE)
  BOLD <- x$BOLD
  Q2_est <- x$Q2
  rm(x)

  # Center and scale `BOLD` again, but do not detrend again. -------------------
  BOLD <- norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=FALSE
  )

  # Initialize with the dual regression-based estimate -------------------------
  if (verbose) { cat("Computing DR.\n") }

  BOLD_DR <- dual_reg(
    BOLD, template$mean, center_Bcols=FALSE,
    scale=FALSE, detrend_DCT=0, normA=normA
  )

  # ----------------------------------------------------------------------------
  # EM -------------------------------------------------------------------------
  #Three algorithms to choose from:
  #1) FC Template ICA
  #2) Template ICA
  #3) Spatial Template ICA (initialize with standard Template ICA)
  # ----------------------------------------------------------------------------

  #1) FC Template ICA ----------------------------------------------------------
  if(do_FC) {
    if (verbose) { cat("Estimating FC Template ICA model.\n") }

    template_mean = template$mean
    template_var = template$var
    prior_params = c(0.001, 0.001) #alpha, beta (uninformative) -- note that beta is scale parameter in IG but rate parameter in the Gamma

    #EM_FCtemplateICA <- function(...){NULL}
    resultEM <- EM_FCtemplateICA(
      template_mean = template$mean,
      template_var = template$var,
      template_FC = template_FC,
      prior_params, #for prior on tau^2
      BOLD=BOLD,
      AS_0 = BOLD_DR, #initial values for A and S
      maxiter=maxiter,
      epsilon=epsilon,
      verbose=verbose
    )
  } else {

    if (reduce_dim) {
      if (verbose) { cat("Reducing data dimensions.\n") }
      # Reduce data dimensions
      BOLD <- dim_reduce(BOLD, nL)
      err_var <- BOLD$sigma_sq # spw: need to run dim red to get this quantity
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
    } else {
      # [TO DO]: what if just compute eigenvalues? faster, right?
      err_var <- dim_reduce(BOLD, nL)$sigma_sq
      BOLD2 <- BOLD; rm(BOLD)
      theta0 <- list(A = BOLD_DR$A)
      C_diag <- rep(1, nT)
    }

    #2) Template ICA -----------------------------------------------------------
    if (verbose) {
      if (do_spatial) { cat("Initializing with standard Template ICA.\n") }
      if (!do_spatial) { cat("Computing Template ICA.\n") }
    }
    theta00 <- theta0
    theta00$nu0_sq <- err_var
    resultEM <- EM_templateICA.independent(template$mean,
                                           template$var,
                                           BOLD=BOLD2,
                                           theta0=theta00,
                                           C_diag=C_diag,
                                           maxiter=maxiter,
                                           usePar=usePar,
                                           epsilon=epsilon,
                                           verbose=verbose)
    if (reduce_dim) { resultEM$A <- Hinv %*% resultEM$theta_MLE$A }
    class(resultEM) <- 'tICA'

    #3) Spatial Template ICA ---------------------------------------------------
    if(do_spatial){
      resultEM_tICA <- resultEM
      theta0$kappa <- rep(kappa_init, nL)
      if(verbose) cat('ESTIMATING SPATIAL MODEL\n')
      t000 <- Sys.time()
      resultEM <- EM_templateICA.spatial(template$mean,
                                         template$var,
                                         meshes,
                                         BOLD=BOLD2,
                                         theta0,
                                         C_diag,
                                         maxiter=maxiter,
                                         usePar=usePar,
                                         epsilon=epsilon,
                                         verbose=verbose)
      #common_smoothness=common_smoothness)
      print(Sys.time() - t000)

      #organize estimates and variances in matrix form
      resultEM$subjICmean <- matrix(resultEM$subjICmean, ncol=nL)
      resultEM$subjICse <- sqrt(matrix(diag(resultEM$subjICcov), ncol=nL))
      resultEM$A <- Hinv %*% resultEM$theta_MLE$A

      resultEM$result_tICA <- resultEM_tICA
      class(resultEM) <- 'stICA'
    }
    if (usePar) { doParallel::stopImplicitCluster() }
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

  # Params, formatted as length-one character vectors to put in "xifti" metadata
  tICA_params <- list(
    time_inds=time_inds, center_Bcols=center_Bcols,
    scale=scale, detrend_DCT=detrend_DCT, normA=normA,
    Q2=Q2, Q2_max=Q2_max, Q2_est = Q2_est,
    brainstructures=brainstructures,
    tvar_method=tvar_method,
    spatial_model=do_spatial,
    rm_mwall=rm_mwall,
    reduce_dim=reduce_dim,
    maxiter=maxiter,
    epsilon=epsilon,
    kappa_init=kappa_init
  )
  tICA_params <- lapply(
    tICA_params,
    function(x) {
      if (is.null(x)) { x <- "NULL"};
      paste0(as.character(x), collapse=" ")
    }
  )

  # Format output.
  if (use_mask2) {
    resultEM$subjICmean <- fMRItools::unmask_mat(resultEM$subjICmean, mask2)
    resultEM$subjICse <- fMRItools::unmask_mat(resultEM$subjICse, mask2)
  }

  if (FORMAT %in% c("CIFTI", "GIFTI") && !is.null(xii1)) {
    xiiL <- ciftiTools::select_xifti(xii1, rep(1, nL))
    if (grepl("all", IC_inds)) {
      IC_inds <- seq(nL)
    } else {
      IC_inds <- strsplit(IC_inds, " ")[[1]]
      if (length(IC_inds) != nL) { IC_inds <- rep("?", nL) } # TO-DO: improve
    }
    xiiL$meta$cifti$names <- paste("IC", IC_inds)

    # [TO DO]: fMRItools update: simplify this
    if (any(!mask3)) {
      resultEM$subjICmean <- ciftiTools::newdata_xifti(
        xiiL,
        fMRItools::unmask_mat(resultEM$subjICmean, mask3),
      )
      resultEM$subjICse <- ciftiTools::newdata_xifti(
        xiiL,
        fMRItools::unmask_mat(resultEM$subjICse, mask3),
      )
    } else {
      resultEM$subjICmean <- ciftiTools::newdata_xifti(xiiL, resultEM$subjICmean)
      resultEM$subjICse <- ciftiTools::newdata_xifti(xiiL, resultEM$subjICse)
    }
    if (FORMAT == "GIFTI") {
      # Apply `mask2`.
      resultEM$subjICmean <- ciftiTools::move_to_mwall(resultEM$subjICmean)
      resultEM$subjICse <- ciftiTools::move_to_mwall(resultEM$subjICse)
      mask2 <- NULL
    }
    if (do_spatial) {
      resultEM$result_tICA$subjICmean <- ciftiTools::newdata_xifti(xiiL, resultEM$result_tICA$subjICmean)
      resultEM$result_tICA$subjICse <- ciftiTools::newdata_xifti(xiiL, resultEM$result_tICA$subjICse)
    }
    class(resultEM) <- 'tICA.cifti'

  } else if (FORMAT == "NIFTI") {
    resultEM$subjICmean <- RNifti::asNifti(
      fMRItools::unvec_vol(resultEM$subjICmean, mask, fill=NA)
    )
    resultEM$subjICse <- RNifti::asNifti(
      fMRItools::unvec_vol(resultEM$subjICse, mask, fill=NA)
    )
    if (do_spatial) {
      resultEM$result_tICA$subjICmean <- RNifti::asNifti(
        fMRItools::unvec_vol(resultEM$result_tICA$subjICmean, mask, fill=NA)
      )
      resultEM$result_tICA$subjICse <- RNifti::asNifti(
        fMRItools::unvec_vol(resultEM$result_tICA$subjICse, mask, fill=NA)
      )
    }
    resultEM$mask_nii <- mask
    # resultEM$mask <- mask
    class(resultEM) <- 'tICA.nifti'
  } else {
    class(resultEM) <- 'tICA.matrix'
  }

  resultEM$mask <- mask2
  resultEM$params <- tICA_params
  resultEM
}
