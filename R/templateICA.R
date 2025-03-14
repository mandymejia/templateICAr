#' Template ICA
#'
#' Perform template independent component analysis (ICA) using variational Bayes
#'  (VB) or expectation-maximization (EM).
#'
#' @param BOLD Vector of subject-level fMRI data in one of the following
#'  formats: CIFTI file paths, \code{"xifti"} objects, NIFTI file paths,
#'  \code{"nifti"} objects, or \eqn{V \times T} numeric matrices, where \eqn{V}
#'  is the number of data locations and \eqn{T} is the number of timepoints.
#'
#'  If multiple BOLD data are provided, they will be independently centered,
#'  scaled, detrended (if applicable), and denoised (if applicable). Then they
#'  will be concatenated together followed by computing the initial dual
#'  regression estimate.
#' @param template Template estimates in a format compatible with \code{BOLD},
#'  from \code{\link{estimate_template}}.
#' @param tvar_method Which calculation of the template variance to use:
#'  \code{"non-negative"} (default) or \code{"unbiased"}. The unbiased template
#'  variance is based on the assumed mixed effects/ANOVA model, whereas the
#'  non-negative template variance adds to it to account for greater potential
#'  between-subjects variation. (The template mean is the same for either choice
#'  of \code{tvar_method}.)
#' @param GSR Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default:
#'  \code{"template"}, to use the same option used for estimation of the
#'  \code{template}.
#' @param scale \code{"global"}, \code{"local"}, or \code{"none"}.
#'  Global scaling will divide the entire data matrix by the mean image standard
#'  deviation (\code{mean(sqrt(rowVars(BOLD)))}). Local scaling will divide each
#'  data location's time series by its estimated standard deviation. Default:
#'  \code{"template"}, to use the same option used for estimation of the
#'  \code{template}.
#' @param scale_sm_surfL,scale_sm_surfR,scale_sm_FWHM Only applies if
#'  \code{scale=="local"} and \code{BOLD} represents CIFTI-format data. To
#'  smooth the standard deviation estimates used for local scaling, provide the
#'  surface geometries along which to smooth as GIFTI geometry files or
#'  \code{"surf"} objects, as well as the smoothing FWHM (default:
#'  \code{"template"} to use the same option used for estimation of the
#'  \code{template}).
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
#' @param nuisance (Optional) Signals to regress from the data, given as a
#'  numeric matrix with the same number of rows as there are volumes in the
#'  \code{BOLD} data. If multiple \code{BOLD} sessions are provided,
#'  this argument can be a list to use different nuisance regressors for
#'  different sessions. Nuisance regression is performed as a first step, before
#'  centering, scaling, and denoising. An intercept column will automatically be
#'  added to \code{nuisance}. If \code{NULL}, no extra nuisance signals will be
#'  regressed from the data, but a nuisance regression will still be used if
#'  warranted by \code{scrub} or \code{hpf}.
#'
#' @param scrub (Optional) A numeric vector of integers indicating the indices
#'  of volumes to scrub from the BOLD data. (List the volumes to remove, not the
#'  ones to keep.) If multiple \code{BOLD} sessions are provided, this
#'  argument can be a list to remove different volumes for different sessions.
#'  Scrubbing is performed within nuisance regression by adding a spike
#'  regressor to the nuisance design matrix for each volume to scrub. If
#'  \code{NULL}, do not scrub.
#'
#' @param TR,hpf These arguments control detrending. \code{TR} is the temporal
#'  resolution of the data, i.e. the time between volumes, in seconds;
#'  \code{hpf} is the frequency of the high-pass filter, in Hertz. Detrending
#'  is performed via nuisance regression of DCT bases. Default:
#'  \code{"template"} to use the values from the template. Be sure to set the
#'  correct \code{TR} if it's different for the new data compared to the data
#'  used in \code{estimate_template}.
#'
#'  Note that if multiple \code{BOLD} sessions are provided, their
#'  \code{TR} and \code{hpf} must be the same; both arguments accept only one
#'  value.
#' @param Q2,Q2_max Denoise the BOLD data? Denoising is based on modeling and
#'  removing nuisance ICs. It may result in a cleaner estimate for smaller
#'  datasets, but it may be unnecessary (and time-consuming) for larger datasets.
#'
#'  Set \code{Q2} to control denoising: use a positive integer to specify the
#'  number of nuisance ICs, \code{NULL} to have the number of nuisance ICs
#'  estimated by PESEL, or zero to skip denoising.
#'
#'  If \code{is.null(Q2)}, use \code{Q2_max} to specify the maximum number of
#'  nuisance ICs that should be estimated by PESEL. \code{Q2_max} must be less
#'  than \eqn{T * .75 - Q} where \eqn{T} is the number of timepoints in BOLD
#'  and \eqn{Q} is the number of group ICs. If \code{NULL}, \code{Q2_max} will
#'  be set to \eqn{T * .50 - Q}, rounded.
#'
#'  The defaults for both arguments is \code{"template"}, to use the same option
#'  used for estimation of the \code{template}.
#' @param covariates Numeric vector of covariates to take into account for model
#'  estimation. Names should give the name of each variable. The covariates must
#'  match those of the template. Default: \code{NULL} (no covariates).
#'  NOTE: Not implemented yet.
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI
#'  file paths. This is a character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{"template"} to use the same brainstructures present in the
#'  \code{template}).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI
#'  file paths or \code{"nifti"} objects. This is a brain map formatted as a
#'  binary array of the same spatial dimensions as the fMRI data, with
#'  \code{TRUE} corresponding to in-mask voxels.
#' @param time_inds Subset of fMRI BOLD volumes to include in analysis, as a
#'  vector of numeric integers. If \code{NULL} (default), use all volumes.
#'  Subsetting is performed before any centering, scaling, detrending, and
#'  denoising. If multiple \code{BOLD} are provided, \code{time_inds} can be a
#'  list of the same length, with each element being a vector of numeric
#'  integers indicating the timepoints to include for each BOLD.
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
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{"template"} to use the same brainstructures present in the
#'  \code{template}). Variance is calculated on the original data, before
#'  any normalization. Set to \code{0} to avoid removing locations due to
#'  low variance.
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
#'  estimation, but may result in more accurate results. Ignored for FC template
#'  ICA
#' @param method_FC Variational Bayes (VB) method for FC template ICA model:
#'  \code{"VB1"} (default) uses a conjugate Inverse-Wishart prior for the cor(A);
#'  \code{"VB2"} draws samples from p(cor(A)) to emulate the population distribution
#'  using a combination of Cholesky, SVD, and random pivoting.
#'  \code{"none"} Uses standard template ICA without FC prior
#' @param maxiter Maximum number of EM or VB iterations. Default: \code{100}.
#' @param miniter Minimum number of EM or VB iterations. Default: \code{3}.
#' @param epsilon Smallest proportion change between iterations. Default:
#'  \code{.001}.
# @param eps_inter Intermediate values of epsilon at which to save results (used
#  to assess benefit of more stringent convergence rules). Default:
#  \code{10e-2} to \code{10e-5}. These values should be in decreasing order
#  (larger to smaller error) and all values should be between zero and
#  \code{epsilon}.
#' @param kappa_init Starting value for kappa. Default: \code{0.2}.
#' @param usePar Parallelize the computation? Default: \code{TRUE}. Can be the
#' number of cores to use or \code{TRUE}, which will use the number available minus two.
#' @param PW Prewhiten to account for residual autocorrelation? Default: \code{FALSE}.
#' Only affects FC template ICA models.
# @param tESS_correction Take into account effective sample size due to temporal autocorrelation?
# Default: \code{FALSE}. Only affects FC template ICA models.
# @param sESS_correction To take into account effective sample size due to spatial correlation,
# provide list of mesh vertices, faces, and data indices (in that order). Default: \code{FALSE}.
# Only affects FC template ICA models.
#' @param verbose If \code{TRUE}, display progress of algorithm
# @param common_smoothness If \code{TRUE}. use the common smoothness version
#  of the spatial template ICA model, which assumes that all IC's have the same
#  smoothness parameter, \eqn{\kappa}
#'
#' @return A (spatial) template ICA object, which is a list containing:
#'  \code{subjICmean}, the \eqn{V \times L} estimated independent components
#'  \strong{S}; \code{subjICse}, the standard errors of \strong{S}; the
#'  \code{mask} of locations without template values due to too many low
#'  variance or missing values; the \code{nuisance} design matrix or matrices if
#'  applicable; and the function \code{params} such as the type of scaling and
#'  detrending performed.
#'
#'  If \code{BOLD} represented CIFTI or NIFTI data, \code{subjICmean} and
#'  \code{subjICse} will be formatted as \code{"xifti"} or \code{"nifti"}
#'  objects, respectively.
#'
#' @export
#'
# @importFrom INLA inla inla.spde.result inla.pardiso.check inla.setOption
#' @importFrom fMRItools infer_format_ifti_vec unmask_mat unvec_vol is_1 is_posNum dual_reg
#' @importFrom fMRIscrub flags_to_nuis_spikes
#' @importFrom stats optim
#' @importFrom matrixStats rowVars
#' @importFrom Matrix bandSparse Matrix
#'
#' @examples
#' \dontrun{
#'  tm <- estimate_template(cii1_fnames, cii2_fnames, gICA_fname)
#'  templateICA(newcii_fname, tm, spatial_model=TRUE, resamp_res=2000)
#' }
templateICA <- function(
  BOLD, template,
  tvar_method=c("non-negative", "unbiased"),
  #tinds=NULL,
  scale=c("template", "global", "local", "none"),
  scale_sm_surfL=NULL,
  scale_sm_surfR=NULL,
  scale_sm_FWHM="template",
  nuisance=NULL,
  scrub=NULL,
  TR=NULL, hpf="template",
  GSR="template",
  Q2="template",
  Q2_max="template",
  covariates=NULL,
  brainstructures="template",
  mask=NULL,
  time_inds=NULL,
  varTol="template",
  spatial_model=NULL,
  resamp_res=NULL,
  rm_mwall=TRUE,
  reduce_dim=FALSE,
  method_FC=c("VB1", "VB2", "none"),
  maxiter=100,
  miniter=3,
  epsilon=0.001,
  #eps_inter=NULL,
  kappa_init=0.2,
  #common_smoothness=TRUE,
  usePar=TRUE,
  PW=FALSE,
  verbose=TRUE){

  t0 <- Sys.time()

  # Check arguments ------------------------------------------------------------

  # Read in template and set arguments designated to match the template.
  if (is.character(template)) { template <- readRDS(template) }
  TFORMAT <- class(template)[grepl("template", class(template))]
  if (length(TFORMAT) != 1) { stop("`template` is not a template.") }
  TFORMAT <- toupper(gsub("template.", "", TFORMAT, fixed=TRUE))

  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("template", "global", "local", "none"))
  if (scale == "template") { scale <- template$params$scale }
  if (scale_sm_FWHM == "template") {
    scale_sm_FWHM <- template$params$scale_sm_FWHM
  }
  if (is.null(TR)) { TR <- "from_xifti_metadata" }
  stopifnot(length(TR)==1) # be explicit, diff TR for multi-BOLD is not allowed
  if (TR == "template") { TR <- template$params$TR }
  stopifnot(length(hpf)==1) # be explicit, diff TR for multi-BOLD is not allowed
  if (hpf == "template") { hpf <- template$params$hpf }
  if (GSR == "template") {
    GSR <- template$params$GSR
  }
  if (Q2 == "template") { Q2 <- template$params$Q2 }
  if (Q2_max == "template") { Q2_max <- template$params$Q2_max }
  if (!is.null(template$params$covariate_names)) {
    covariate_names <- template$params$covariate_names
    nC <- length(covariate_names)
    if (is.null(covariates)) {
      stop("These covariates were used during template estimation: ", 
        paste0(covariate_names, collapse=", "), ". They must also be provided ", 
        "to `templateICA` with the `covariates` argument.")
    }
    stopifnot(is.numeric(covariates) && is.vector(covariates))
    stopifnot(length(covariates) == length(covariate_names))
    if(!all(names(covariates) == covariate_names)) {
      stop("These covariates were used during template estimation: ", 
        paste0(covariate_names, collapse=", "), ". The same covariates must ", 
        "also be provided to `templateICA` with the `covariates` argument. ",
        "However, the names for `covariates` provided here differ.")
    }
  } else {
    if (!is.null(covariates)) { 
      stop("`covariates` is not `NULL`, yet none were used during template ",
      "estimation. Any covariates must match those used in the template.")
    }
    covariate_names <- NULL
    nC <- 0
  }
  brainstructures <- match.arg(
    brainstructures,
    c("template", "left", "right", "subcortical", "all"),
    several.ok=TRUE
  )
  if (length(brainstructures)==1 && brainstructures == "template") {
    brainstructures <- template$params$brainstructures
  }
  if (varTol == "template") { varTol <- template$params$varTol }

  # Remaining simple argument checks.
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
  stopifnot(is_1(scale_sm_FWHM, "numeric"))
  if (is.null(hpf)) { hpf <- 0 }
  if (TR!= "from_xifti_metadata") { stopifnot(fMRItools::is_posNum(TR)) }
  stopifnot(fMRItools::is_posNum(hpf, zero_ok=TRUE))
  stopifnot(is_1(GSR, "logical"))
  if (!is.null(Q2)) { stopifnot(is_posNum(Q2, zero_ok=TRUE)) } # Q2_max checked later.
  stopifnot(is_posNum(varTol))
  if (isFALSE(spatial_model)) { spatial_model <- NULL }
  if (!is.null(resamp_res)) {
    stopifnot(is_posNum(resamp_res) && round(resamp_res) == resamp_res)
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
  }
  stopifnot(is_1(rm_mwall, "logical"))
  stopifnot(is_1(reduce_dim, "logical"))
  method_FC <- match.arg(method_FC, c("VB1", "VB2", "none"))
  stopifnot(is_posNum(maxiter))
  stopifnot(is_posNum(miniter))
  stopifnot(miniter <= maxiter)
  stopifnot(is_posNum(epsilon))
  # if (!is.null(eps_inter)) {
  #   stopifnot(is.numeric(eps_inter) && all(diff(eps_inter) < 0))
  #   stopifnot(eps_inter[length(eps_inter)]>0 && eps_inter[1]>epsilon)
  # }
  if (!is.null(kappa_init)) { stopifnot(is_posNum(kappa_init)) }
  stopifnot(is_1(usePar, "logical") || is_1(usePar, "numeric"))
  stopifnot(is_1(verbose, "logical"))

  # `usePar`
  if (!isFALSE(usePar)) {
    check_parallel_packages()
    cores <- parallel::detectCores()
    if (isTRUE(usePar)) { nCores <- cores[1] - 2 } else { nCores <- usePar; usePar <- TRUE }
    if (nCores < 2) {
      usePar <- FALSE
    } else {
      #cluster <- parallel::makeCluster(nCores)
      doParallel::registerDoParallel(nCores)
    }
  }

  # `BOLD` ---------------------------------------------------------------------
  cat("\n")
  # Determine the format of `BOLD`.
  # [TO DO]: more elegant way? b/c list of xifti vs single xifti...
  format <- suppressWarnings(fMRItools::infer_format_ifti(BOLD))[1]
  if (is.na(format)) {
    format <- fMRItools::infer_format_ifti_vec(BOLD)[1]
  }
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
    if ("all" %in% brainstructures) { brainstructures <- c("left", "right", "subcortical") }
    do_left <- "left" %in% brainstructures
    do_right <- "right" %in% brainstructures
    do_sub <- "subcortical" %in% brainstructures

    if (format == "xifti") {
      if (!requireNamespace("ciftiTools", quietly = TRUE)) {
        stop("Package \"ciftiTools\" needed to work with xifti data. Please install it.", call. = FALSE)
      }

      for (bb in seq(nN)) {
        if (!do_left && !is.null(BOLD[[bb]]$data$cortex_left)) {
          BOLD[[bb]] <- ciftiTools::remove_xifti(BOLD[[bb]], "cortex_left")
        }
        if (!do_right && !is.null(BOLD[[bb]]$data$cortex_right)) {
          BOLD[[bb]] <- ciftiTools::remove_xifti(BOLD[[bb]], "cortex_right")
        }
        if (!do_sub && !is.null(BOLD[[bb]]$data$subcort)) {
          BOLD[[bb]] <- ciftiTools::remove_xifti(BOLD[[bb]], "subcortical")
        }
      }
    }
  }

  # Read in CIFTI, GIFTI, or NIFTI files.
  if (verbose && format %in% c("CIFTI", "GIFTI", "NIFTI", "RDS")) {
    cat("Reading BOLD data.\n")
  }

  # (Need to do now rather than later, so that CIFTI resolution info can be used.)
  if (format == "CIFTI") {
    for (bb in seq(nN)) {
      if (is.character(BOLD[[bb]])) {
        BOLD[[bb]] <- ciftiTools::read_cifti(
          BOLD[[bb]], resamp_res=resamp_res,
          brainstructures=brainstructures,
          verbose=FALSE
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

  # Check template format matches BOLD format.
  if (TFORMAT != FORMAT) {
    stop("The BOLD format is '", FORMAT, ",' but the template format is '", TFORMAT, ".'")
  }

  # Check that parameters match.
  if (is.null(template$params)) {
    # warning("Old template does not have `params` metadata, so the parameters can't be checked.")
    NULL
  } else {
    pmatch <- c(
      scale=scale,
      scale_sm_FWHM=scale_sm_FWHM,
      hpf=hpf,
      GSR=GSR,
      # Q2=Q2, Q2_max=Q2_max,
      varTol=varTol
    )
    for (pp in seq(length(pmatch))) {
      pname <- names(pmatch)[pp]
      if (pmatch[pname] != template$params[[pname]]) {
        warning(paste0(
          "The `", pname, "` parameter was ",
          template$params[[pname]], " for the template, but is ",
          pmatch[pname], " here. Generally, these should match. (Proceeding anyway.)\n"
        ))
      }
    }
  }

  #check for FC template
  do_FC <- FALSE; template_FC <- NULL
  if(method_FC == "VB1"){
    if('FC' %in% names(template$template)) {
      do_FC <- TRUE
      template_FC <- template$template$FC
      template$template$FC <- template$template$FC_Chol <- NULL
    } else {
      warning("FC information not available in `template`.  I will set `method_FC` to 'none' and perform standard template ICA.")
      method_FC <- "none"
    }
  } else if(method_FC == "VB2"){
    if('FC_Chol' %in% names(template$template)) {
      do_FC <- TRUE
      template_FC <- template$template$FC_Chol
      template$template$FC <- template$template$FC_Chol <- NULL
    } else {
      warning("Cholesky FC information not available in `template`.  I will set `method_FC` to 'none' and perform standard template ICA.")
      method_FC <- "none"
    }
  }

  # Get `nI`, `nL`, and `nV`.
  # Check brainstructures and resamp_res if CIFTI, where applicable.
  # Convert to CIFTI if GIFTI.
  # Check mask if NIFTI.
  xii1 <- NULL
  if (FORMAT == "CIFTI") {

    # Check `resamp_res`.
    if (!is.null(resamp_res)) {
      template <- resample_template(template, resamp_res=resamp_res)
    }

    xii1 <- template$dat_struct
    # if (!is.null(resamp_res)) {
    #   if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    #     stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    #   }
    #   #xii1 <- ciftiTools::resample_xifti(xii1, resamp_res = resamp_res)
    # }
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
          "Ignoring ", ifelse(length(bs_missing) > 1, "these brain structures", "this brain structure"),
          " in the template: ", paste(bs_missing, collapse=", "), ".\n"
        ))
      }
      template <- removebs_template(template, bs_missing)
      xii1 <- ciftiTools::remove_xifti(xii1, bs_missing)
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
    tds_dim <- dim(template$mask_input) # has an empty last dim
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
    var = template$template[[tvar_name]],
    sigma_sq0 = template$sigma_sq0
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
    if (!any(grepl('FAILURE',flag))) { INLA::inla.setOption(smtp='pardiso') }

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
          if (length(unique(res))==1) {
            res <- res[1]
          } else if (any(res==0)) {
            res <- res[res!=0]
          } else {
            # [TO DO]: handle the case of unequal resolutions? this won't work
            #   but is a placeholder.
            res <- max(res)
          }
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
  ## Vectorize and mask data; get dimensions -----------------------------------

  # Get each entry of `BOLD` as a data matrix or array.
  # If `FORMAT==CIFTI`, get `TR` if needed.
  if (FORMAT %in% c("CIFTI", "GIFTI")) {
    for (bb in seq(nN)) {
      if (ciftiTools::is.xifti(BOLD[[bb]])) {
        # Get `TR` if not have yet; or check that it matches the current `TR`
        if (!is.null(BOLD[[bb]]$meta$cifti$time_step)) {
          if (TR == "from_xifti_metadata") {
            if (verbose) {
              message(
                "Setting TR=", BOLD[[bb]]$meta$cifti$time_step,
                ifelse(nN==1, " based on BOLD.", paste0(" based on BOLD #", bb, "."))
              )
            }
            TR <- BOLD[[bb]]$meta$cifti$time_step
          } else if (TR != BOLD[[bb]]$meta$cifti$time_step) {
            warning(
              "The `time_step` metadata for BOLD #", bb, " is ",
              BOLD[[bb]]$meta$cifti$time_step,
              ", but TR=", TR, " based on either the provided `TR` argument",
              " or an earlier BOLD. Proceeding anyway with TR=", TR, "."
            )
          }
        }
        BOLD[[bb]] <- as.matrix(BOLD[[bb]])
      }
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
      time_inds_bb <- if (is.list(time_inds)) {
        time_inds[[bb]]
      } else {
        time_inds
      }
      stopifnot(is.numeric(time_inds_bb))
      if (!all(time_inds_bb %in% seq(nT))) {
        warning('Not all `time_inds` available.') # [TO DO]: improve
        time_inds_bb <- intersect(time_inds_bb, seq(nT))
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
    cat('Number of BOLD sessions:       ', nN, "\n")
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
    template$sigma_sq0 <- template$sigma_sq0[mask2]
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
  use_mask3 <- any(!mask3)

  if (use_mask3) {
    if (do_spatial) {
      stop('Not supported yet: flat or NA voxels in data, after applying template mask, with spatial model.')
    }

    # [NOTE] For same results, would have needed to also update "A" matrix
    #   (projection from mesh to data locations)

    if(verbose) {
      cat(paste0(
        'In total, excluding ', sum(!mask3),
        ' locations from analysis due to flat or NA values in BOLD.\n'
      ))
    }
    template_orig <- template
    BOLD <- lapply(BOLD, function(x){x[mask3,]})
    dBOLDs <- lapply(BOLD, dim); dBOLD <- dBOLDs[[1]]
    template$mean <- template$mean[mask3,]
    template$var <- template$var[mask3,]
    template$sigma_sq0 <- template$sigma_sq0[mask3]
    if (use_mask2) { mask2[mask2][!mask3] <- FALSE }
  }

  ## Nuisance regression and scrubbing -----------------------------------------
  if (verbose) { cat("\n") }
  if (verbose) { cat("Pre-processing BOLD data.\n") }

  if (is.list(nuisance)) { stopifnot(length(nuisance)==nN) }
  if (is.list(scrub)) { stopifnot(length(scrub)==nN) }

  add_to_nuis <- function(x, nuis) {
    if (is.null(nuis)) { x } else { cbind(x, nuis) }
  }

  nmat <- vector("list", nN)
  nDCT <- if (hpf==0) { NULL } else { vector("numeric", nN) } # could return this
  nT_pre <- nT
  for (nn in seq(nN)) {
    if (verbose && nN > 1) { cat(paste0("Session ", nn, ":")) }
    # Collect nuisance matrix columns.
    nmat[nn] <- list(NULL)
    ## `nuisance`
    nuisance_nn <- if (is.list(nuisance)) { nuisance[[nn]] } else { nuisance }
    if (!is.null(nuisance_nn)) {
      stopifnot(is.numeric(nuisance_nn) && is.matrix(nuisance_nn))
      stopifnot(nrow(nuisance_nn) == nT[nn])
      if (verbose && nN > 1) { cat("\t") }
      if (verbose) { cat("Using", ncol(nuisance_nn), "regressors from `nuisance`.\n") }
      nmat[[nn]] <- add_to_nuis(nuisance_nn, nmat[[nn]])
    }
    ## `scrub`
    scrub_nn <- NULL
    if (!is.null(scrub)) {
      scrub_nn <- if (is.list(scrub)) { scrub[[nn]] } else { scrub }
      if (is.logical(scrub_nn)) { scrub_nn <- which(scrub_nn) }
      if (length(scrub_nn) > 0) {
        scrub_nn_mat <- fMRIscrub::flags_to_nuis_spikes(scrub_nn, nT[nn])
        if (verbose && nN > 1) { cat("\t") }
        if (verbose) { cat("Scrubbing", ncol(scrub_nn_mat), "volumes.\n") }
        nmat[[nn]] <- add_to_nuis(scrub_nn_mat, nmat[[nn]])
      } else {
        scrub_nn <- NULL
      }
    }
    ## DCT
    if (hpf != 0) {
      if (TR=="from_xifti_metadata") {
        stop("`hpf!=0`, but `TR`` was neither provided nor able to be inferred from the data. Please provide `TR`.")
      }
      nDCT[nn] <- round(dct_convert(nT[nn], TR=TR, f=hpf))
      if (verbose && nN > 1) { cat("\t") }
      if (verbose) { cat("Using", nDCT[nn], "DCT bases.\n") }
      nmat[[nn]] <- add_to_nuis(dct_bases(nT[nn], nDCT[nn]), nmat[[nn]])
    }
    # Perform nuisance regression, if applicable.
    if (!is.null(nmat[[nn]])) {
      nmat[[nn]] <- cbind(1, nmat[[nn]])
      if (verbose && nN > 1) { cat("\t") }
      if (verbose) { cat("Doing nuisance regression with", ncol(nmat[[nn]]), "total regressors.\n") }
      BOLD[[nn]] <- nuisance_regression(BOLD[[nn]], nmat[[nn]])
    }
    # Drop scrubbed volumes, if applicable.
    if (!is.null(scrub_nn)) {
      BOLD[[nn]] <- BOLD[[nn]][,-scrub_nn]
      dBOLDs <- lapply(BOLD, dim)
      nT <- vapply(dBOLDs, function(x){x[ldB]}, 0)
      nTmin <- min(nT)
    }
  }
  if (sum(nT) != sum(nT_pre)) {
    if (verbose && nN > 1) { cat("\t") }
    if (verbose) { cat('Timepoints after scrubbing:    ', sum(nT), "\n") }
  }

  if (all(vapply(nmat, is.null, FALSE))) { nmat <- NULL }

  ## Center and scale `BOLD` ---------------------------------------------------
  if (verbose) {
    cat("Normalizing BOLD: centering location timecourses")
    if (GSR) { cat(", centering volumes (GSR)") }
    if (scale != "none") { cat(",", scale, "scaling") }
    cat(".\n")
  }

  if (!is.null(xii1) && scale=="local" && scale_sm_FWHM > 0) {
    xii1 <- ciftiTools::add_surf(xii1, surfL=scale_sm_surfL, surfR=scale_sm_surfR)
  }

  mask2and3 <- if (use_mask2) { mask2 } else { mask3 } # [TO DO] patch???

  BOLD <- lapply(BOLD, norm_BOLD,
    center_rows=TRUE, center_cols=GSR,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    scale_sm_xifti_mask=mask2and3,
    hpf=0
  )

  ## Estimate and subtract nuisance ICs ----------------------------------------
  if (verbose) { cat("Removing nuisance ICs.\n") }

  Q2_est <- vector("numeric", nN)
  for (nn in seq(nN)) {
    x <- rm_nuisIC(
      BOLD[[nn]], template_mean=template$mean, Q2=Q2, Q2_max=Q2_max,
      verbose=verbose, return_Q2=TRUE
    )
    BOLD[[nn]] <- x$BOLD
    Q2_est[nn] <- x$Q2
  }
  rm(x)

  ## Center and scale `BOLD` again to ensure mean zero and correct scaling -----
  if (verbose) {
    cat("Normalizing BOLD again: centering location timecourses")
    if (scale != "none") { cat(",", scale, "scaling") }
    cat(".\n")
  }


  BOLD <- lapply(BOLD, norm_BOLD,
    center_rows=TRUE, center_cols=FALSE,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    scale_sm_xifti_mask=mask2and3,
    hpf=0
  )

  ## Concatenate the data. -----------------------------------------------------
  if (verbose && length(nT) > 1) { cat("Concatenating sessions.\n") }

  BOLD <- do.call(cbind, BOLD)
  nT <- sum(nT)

  ## Scale by residual SD from template estimation -----------------------------
  if (verbose) { cat("Scaling by residual SD from template estimation.\n") }
  if(length(template$sigma_sq0) != nrow(BOLD)) stop('Check length of sigma_sq')
  rescale <- sqrt(template$sigma_sq0)/mean(sqrt(template$sigma_sq0), na.rm=TRUE)
  rescale <- matrix(rescale, nrow=length(template$sigma_sq0), ncol=ncol(BOLD))
  BOLD <- BOLD / rescale

  # Initialize with the dual regression-based estimate -------------------------
  if (verbose) { cat("Computing DR.\n") }

  BOLD_DR <- dual_reg(
    BOLD, template$mean, GSR=FALSE,
    scale=FALSE, hpf=0
  )

  t1 <- Sys.time() - t0

  # Bayesian Computation -------------------------------------------------------
  #Three algorithms to choose from:
  #1) Template ICA
  #2) FC Template ICA (EM or VB)
  #3) Spatial Template ICA (initialize with standard Template ICA)

  verbose0 <- verbose
  if (verbose) {
    if (do_spatial | do_FC) {
      cat("Initializing with standard Template ICA.\n")
      verbose <- FALSE
    }
    if (!do_spatial & !do_FC) { cat("Computing Template ICA.\n") }
  }

  if(do_spatial) {
    if(!reduce_dim) warning("Setting reduce_dim to TRUE for spatial template ICA")
    reduce_dim <- TRUE
  }

  if(do_FC) {
    if(reduce_dim) warning("Setting reduce_dim to FALSE for FC template ICA")
    reduce_dim <- TRUE #only temporary, for initilizing with template ICA, then will set to FALSE
  }


  ## 1) Template ICA -----------------------------------------------------------
  if (reduce_dim) {
    if (verbose) { cat("Reducing data dimensions.\n") }
    # Reduce data dimensions
    BOLD_PCA <- dim_reduce(BOLD, nL)
    err_var <- BOLD_PCA$sigma_sq # spw: need to run dim red to get this quantity
    BOLD2 <- BOLD_PCA$data_reduced
    H <- BOLD_PCA$H
    Hinv <- BOLD_PCA$H_inv
    # In original template ICA model nu^2 is separate
    #   for spatial template ICA it is part of C
    C_diag <- BOLD_PCA$C_diag
    if (do_spatial) { C_diag <- C_diag * (BOLD_PCA$sigma_sq) } #(nu^2)HH' in paper
    rm(BOLD_PCA)
    # Apply dimension reduction
    HA <- H %*% BOLD_DR$A
    theta0 <- list(A = HA)
    # #initialize residual variance --- no longer do this, because we use dimension reduction-based estimate
    # theta0$nu0_sq = dat_list$sigma_sq
    # if(verbose) paste0('nu0_sq = ',round(theta0$nu0_sq,1)))
  } else {
    # [TO DO]: what if just compute eigenvalues? faster, right?
    err_var <- dim_reduce(BOLD, nL)$sigma_sq
    BOLD2 <- BOLD
    theta0 <- list(A = BOLD_DR$A)
    C_diag <- rep(1, nT)
    H <- Hinv <- NULL
  }

  theta00 <- theta0
  theta00$nu0_sq <- err_var
  result <- EM_templateICA.independent(
    template_mean=template$mean,
    template_var=template$var,
    BOLD=BOLD2,
    theta0=theta00,
    C_diag=C_diag,
    H=H, Hinv=Hinv,
    maxiter=maxiter,
    epsilon=epsilon,
    reduce_dim=reduce_dim,
    usePar=usePar,
    verbose=verbose
  )
  if (reduce_dim) { result$A <- Hinv %*% result$theta_MLE$A }
  if (!reduce_dim) { result$A <- result$theta_MLE$A }
  class(result) <- 'tICA'
  #end of standard template ICA estimation

  t2 <- Sys.time() - t0

  t3 <- t4 <- 0 #these will be overwritten if FC-tICA or stICA is run

  ## 2) FC Template ICA --------------------------------------------------------
  if (do_FC) {

    verbose <- verbose0

    reduce_dim <- FALSE
    result_tICA <- result

    if (verbose) { cat("Estimating FC Template ICA Model\n") }
    prior_params = c(0.001, 0.001) #alpha, beta (uninformative) -- note that beta is scale parameter in IG but rate parameter in the Gamma

    #no parallelization implemented for VB1
    if(method_FC=='VB1') {
      usePar <- FALSE
      doParallel::stopImplicitCluster()
    }

    result <- VB_FCtemplateICA(
        template_mean = template$mean,
        template_var = template$var,
        template_FC = template_FC,
        method_FC = method_FC,
        prior_params, #for prior on tau^2
        BOLD=BOLD,
        A0 = result_tICA$A,
        S0 = result_tICA$subjICmean,
        S0_var = (result_tICA$subjICse)^2,
        miniter=miniter,
        maxiter=maxiter,
        epsilon=epsilon,
        #eps_inter=eps_inter,
        usePar = usePar,
        PW = PW,
        verbose=verbose
      )

    result$result_tICA <- result_tICA

    t3 <- Sys.time() - t0

  } # end of FC template ICA estimation

  ## 3) Spatial Template ICA ---------------------------------------------------
  if (do_spatial) {
    result_tICA <- result #this is the standard template ICA result
    theta0$kappa <- rep(kappa_init, nL)
    if(verbose) cat('ESTIMATING SPATIAL MODEL\n')
    t000 <- Sys.time()
    result <- EM_templateICA.spatial(template$mean,
                                        template$var,
                                        meshes,
                                        BOLD=BOLD2,
                                        theta0,
                                        C_diag,
                                        H=H, Hinv=Hinv,
                                        maxiter=maxiter,
                                        usePar=usePar,
                                        epsilon=epsilon,
                                        verbose=verbose)
    #common_smoothness=common_smoothness)
    print(Sys.time() - t000)

    #organize estimates and variances in matrix form
    result$subjICmean <- matrix(result$subjICmean, ncol=nL)
    result$subjICse <- sqrt(matrix(diag(result$subjICcov), ncol=nL))
    result$A <- Hinv %*% result$theta_MLE$A

    result$result_tICA <- result_tICA
    class(result) <- 'stICA'

    t4 <- Sys.time() - t0

  }

  ## Wrapping up ---------------------------------------------------------------
  if (usePar) { doParallel::stopImplicitCluster() }

  # Return DR estimates.
  result$result_DR <- BOLD_DR

  #This part problematic for spatial template ICA, but can bring back
  #for template ICA and FC template ICA.  When we check for bad locations,
  #can return an error only for spatial template ICA.

  #result$keep <- keep
  # #map estimates & templates back to original locations
  # if(sum(!keep)>0){
  #   #estimates
  #   subjICmean <- subjICse <- matrix(nrow=length(keep), ncol=L)
  #   subjICmean[keep,] <- result$subjICmean
  #   subjICse[keep,] <- result$subjICse
  #   result$subjICmean <- subjICmean
  #   result$subjICse <- subjICse
  #   #templates
  #   result$template_mean <- template_mean_orig
  #   result$template_var <- template_var_orig
  # }

  # Params
  tICA_params <- list(
    time_inds=time_inds, GSR=GSR,
    scale=scale, TR=TR, hpf=hpf,
    Q2=Q2, Q2_max=Q2_max, Q2_est=Q2_est,
    covariate_names=covariate_names,
    brainstructures=brainstructures,
    tvar_method=tvar_method,
    spatial_model=do_spatial,
    rm_mwall=rm_mwall,
    reduce_dim=reduce_dim,
    miniter=miniter,
    maxiter=maxiter,
    epsilon=epsilon,
    #eps_inter=eps_inter,
    kappa_init=kappa_init
  )

  # Format output.
  if (use_mask2) {
    result$subjICmean <- fMRItools::unmask_mat(result$subjICmean, mask2)
    result$subjICse <- fMRItools::unmask_mat(result$subjICse, mask2)
  } else if (use_mask3) {
    result$subjICmean <- fMRItools::unmask_mat(result$subjICmean, mask3)
    result$subjICse <- fMRItools::unmask_mat(result$subjICse, mask3)
  }

  if (FORMAT %in% c("CIFTI", "GIFTI") && !is.null(xii1)) {
    xiiL <- ciftiTools::select_xifti(xii1, rep(1, nL))
    xiiL <- ciftiTools::convert_to_dscalar(xiiL, names=paste("IC", IC_inds))
    result$subjICmean <- ciftiTools::newdata_xifti(xiiL, result$subjICmean)
    result$subjICse <- ciftiTools::newdata_xifti(xiiL, result$subjICse)

    if (FORMAT == "GIFTI") {
      # Apply `mask2`. # [TO DO] what?
      result$subjICmean <- ciftiTools::move_to_mwall(result$subjICmean)
      result$subjICse <- ciftiTools::move_to_mwall(result$subjICse)
      mask2 <- NULL
    }

    if (do_spatial) {
      result$result_tICA$subjICmean <- ciftiTools::newdata_xifti(xiiL, result$result_tICA$subjICmean)
      result$result_tICA$subjICse <- ciftiTools::newdata_xifti(xiiL, result$result_tICA$subjICse)
    }
    class(result) <- 'tICA.cifti'

  } else if (FORMAT == "NIFTI") {
    result$subjICmean <- RNifti::asNifti(
      fMRItools::unvec_vol(result$subjICmean, mask, fill=NA)
    )
    result$subjICse <- RNifti::asNifti(
      fMRItools::unvec_vol(result$subjICse, mask, fill=NA)
    )
    if (do_spatial) {
      result$result_tICA$subjICmean <- RNifti::asNifti(
        fMRItools::unvec_vol(result$result_tICA$subjICmean, mask, fill=NA)
      )
      result$result_tICA$subjICse <- RNifti::asNifti(
        fMRItools::unvec_vol(result$result_tICA$subjICse, mask, fill=NA)
      )
    }
    result$mask_nii <- mask
    # result$mask <- mask
    class(result) <- 'tICA.nifti'
  } else {
    class(result) <- 'tICA.matrix'
  }

  result$BOLD <- BOLD
  result$mask <- mask2
  result$nuisance <- nmat
  result$params <- tICA_params

  #record computation time of each algorithm
  result$comptime <- c(as.numeric(t1, units = "secs"),
                       as.numeric(t2, units = "secs"),
                       as.numeric(t3, units = "secs"),
                       as.numeric(t4, units = "secs"))

  names(result$comptime) <- c('DR','tICA','FC-tICA','stICA')

  result
}
