#' Perform group ICA based on CIFTI data
#'
#' @param cifti_fnames Vector of file paths of CIFTI-format fMRI timeseries
#'  (*.dtseries.nii) for which to calculate group ICA
#' @param subjects Use this argument if some subjects have more than one scan.
#'  This is a numeric or factor vector the same length as \code{cifti_fnames}. Scans
#'  for the same subject should have the same value. For example, if there are four
#'  subjects with two visits each, and the scans are ordered with the first subject's
#'  two scans first, then the second subject's two scans next, etc., then this argument
#'  should be \code{rep(seq(4), each=2)}. If there are three subjects and four scans,
#'  with the last two scans belonging to the same subject, this argument should be
#'  \code{c(1,2,3,3)}. If all subjects only have one scan, keep this argument as
#'  \code{NULL} (default).
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param num_PCs Maximum number of PCs to retain in initial subject-level PCA
#' @param num_ICs Number of independent components to identify.
#' @param max_rows_GPCA The maximum number of rows for the matrix on which group
#'  level PCA will be performed.  This matrix is the result of temporal concatenation
#'  and contains \code{length(cifti_fnames) * num_PCs} rows.
#' @param center_Bcols,scale,scale_sm_FWHM,detrend_DCT Center BOLD columns, scale by the
#'  standard deviation, and detrend voxel timecourses? See 
#'  \code{\link{norm_BOLD}}. Normalization is applied separately to each scan.
#'  Defaults: Center BOLD columns, scale by the global standard deviation, but
#'  do not detrend.
#' @param verbose If \code{TRUE}, display progress updates
#' @param out_fname (Optional) If not specified, a xifti object will be returned but
#'  the GICA maps will not be written to a CIFTI file.
#' @param surfL (Optional) File path to a surface GIFTI for the left cortex.
#'  If provided, will be part of xifti result object for visualization in R.
#'  Will also be used to perform any smoothing.
#' @param surfR (Optional) File path to a surface GIFTI for the right cortex.
#'  If provided, will be part of xifti result object for visualization in R.
#'  Will also be used to perform any smoothing.
#' @param smooth Smooth the CIFTI data prior to reading in each file? Default:
#'  \code{TRUE}. Use the following arguments to control the smoothing parameters.
#' @param smooth_surf_FWHM,smooth_vol_FWHM,smooth_zeroes_as_NA,smooth_subcortical_merged
#'  See \code{\link[ciftiTools]{smooth_cifti}}. The defaults here are the same.
#'  Note that \code{smooth_zeroes_as_NA} will control both of the corresponding
#'  surface and volume arguments to \code{\link[ciftiTools]{smooth_cifti}}.
#'  These arguments only apply if \code{smooth}.
#'
#' @importFrom fMRItools match_input sign_flip
#' @importFrom ica icaimax
#'
#' @return \code{out_fname} if a file was written, or the GICA as a \code{"xifti"} object
#'  if not.
#'
#' @export
#'
groupICA.cifti <- function(
  cifti_fnames, subjects=NULL,
  brainstructures=c("left","right"),
  num_PCs=100,
  num_ICs,
  max_rows_GPCA=10000,
  center_Bcols=FALSE, 
  scale=c("global", "local", "none"), scale_sm_FWHM=2,
  detrend_DCT=0,
  verbose=TRUE,
  out_fname=NULL,
  surfL=NULL,
  surfR=NULL,
  smooth=TRUE,
  smooth_surf_FWHM=5,
  smooth_vol_FWHM=5,
  smooth_zeroes_as_NA=FALSE,
  smooth_subcortical_merged=FALSE#,
  #nRuns = 5
  ){

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
  }

  if (!is.null(out_fname)) {
    if(!dir.exists(dirname(out_fname))) stop("`out_fname` directory does not exist.")
    out_fname <- as.character(out_fname)
    if (length(out_fname) > 1) {
      warning("Using first entry of `out_fname`.")
      out_fname <- out_fname[1]
    }
    if (!endsWith(out_fname, ".dscalar.nii")) {
      out_fname <- paste0(out_fname, ".dscalar.nii")
    }
  }

  if (!is.null(subjects)) {
    subjects <- as.factor(subjects)
    if (length(subjects) == length(unique(subjects))) {
      subjects <- NULL
    } else {
      stopifnot(length(subjects) == length(cifti_fnames))
    }
  } else {
    subjects <- as.factor(seq(length(cifti_fnames)))
  }

  if (missing(num_ICs)) { stop("`num_ICs` must be provided.") }

  notthere <- !file.exists(cifti_fnames)
  if (all(notthere)) {
    stop("All files in `cifti_fnames` do not exist.")
  } else if (any(notthere)) {
    warning(
      "There are ", sum(notthere),
      " files in `cifti_fnames` that do not exist. These will be excluded from the group ICA.\n"
    )
    cifti_fnames <- cifti_fnames[!notthere]
    subjects <- subjects[!notthere]
  }
  subjects <- droplevels(subjects)

  Nsub <- length(levels(subjects))
  Nscans <- length(subjects)
  any_multi <- Nscans != Nsub

  if (verbose) {
    cat(paste0("Number of subjects: ", Nsub, "\n"))
    cat(paste0("Number of total scans: ", Nscans, "\n"))
  }

  if (any_multi) {
    if (verbose) {
      cat("Scans per subject: ")
      sub_counts <- unique(table(subjects))
      if (length(sub_counts) < 2) {
        cat(sub_counts); cat("\n")
      } else {
        cat(paste0(min(sub_counts), "-", max(sub_counts), "\n"))
      }
    }
  }

  if (Nsub*num_PCs > max_rows_GPCA) {
    warning(
      "The number of subjects times the number of ",
      "subject-level PCs exceeds the limit of ", max_rows_GPCA,
      ". Some subjects will be excluded from analysis.\n"
    )
  }

  brainstructures <- fMRItools::match_input(
    brainstructures, c("left","right","subcortical","all"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }

  size_bigY <- min(max_rows_GPCA, (num_PCs*Nsub)) #for group-level PCA

  if (smooth) {
    # Warn the user when they don"t provide surfaces, which will require using
    #   the default ciftiTools surfaces (inflated) to smooth on.
    sm_defaultL <- ("left" %in% brainstructures) & is.null(surfL)
    sm_defaultR <- ("right" %in% brainstructures) & is.null(surfR)
    if (sm_defaultL || sm_defaultR) {
      qwarn <- c("left cortex", "right cortex", "cortex")[as.numeric(sm_defaultL + 2*sm_defaultR)]
      message("Using inflated surface to smooth the ", qwarn, " data.\n")
    }
  }

  # READ IN DATA AND PERFORM INITIAL DIMENSION REDUCTION
  next_col <- 1
  for (ii in seq(Nsub)) {

    # Get CIFTI files.
    if (verbose) cat(paste0("Reading data for subject ",ii," of ",Nsub, "\n"))
    fnames_ii <- cifti_fnames[as.numeric(subjects)==ii]

    # Smooth BOLD data.
    # The call returns the path to the smoothed cifti, which will be in a temp dir.
    if (smooth) {
      for (jj in seq(length(fnames_ii))) {
        fnames_ii[jj] <- ciftiTools::smooth_cifti(
          fnames_ii[jj],
          cifti_target_fname=file.path(tempdir(), basename(fnames_ii[jj])),
          surf_FWHM = smooth_surf_FWHM, vol_FWHM = smooth_vol_FWHM,
          surfL_fname=surfL, surfR_fname=surfR,
          subcortical_zeroes_as_NA = smooth_zeroes_as_NA,
          cortical_zeroes_as_NA = smooth_zeroes_as_NA,
          subcortical_merged = smooth_subcortical_merged
        )
      }
    }

    # Read in BOLD data.
    BOLD_ii <- lapply(fnames_ii, ciftiTools::read_cifti, brainstructures=brainstructures)
    # Normalize each scan (keep in `"xifti"` format for `merge_xifti` next).
    BOLD_ii <- lapply(BOLD_ii, function(x){
      ciftiTools::newdata_xifti(x, norm_BOLD(
        as.matrix(x), center_cols=center_Bcols, 
        scale=scale, scale_sm_xifti=ciftiTools::select_xifti(x, 1), scale_sm_FWHM=scale_sm_FWHM,
        detrend_DCT=detrend_DCT
      ))
    })

    # Concatenate and convert to matrix.
    # `merge_xifti` will check that voxels and vertices align;
    #   i.e. that the resolutions are the same.
    BOLD_ii <- as.matrix(ciftiTools::merge_xifti(xifti_list=BOLD_ii))

    # Check the total timepoints for this subject.
    #   Update the running total.
    ntime <- ncol(BOLD_ii)
    num_PCs_ii <- min(num_PCs, ntime-1)
    cols_ii <- next_col - 1 + seq(num_PCs_ii)
    # Quit if adding these next PCs would take the number of PCs past the limit.
    if (max(cols_ii) > size_bigY) { break }
    next_col <- next_col + seq(num_PCs_ii)

    # Check original and reduced dimensions.
    if (ii==1) {
      nvox <- nrow(BOLD_ii)
      bigY <- matrix(NA, nrow=nvox, ncol=size_bigY) #for initializing online PCA
    } else {
      if (nrow(BOLD_ii) != nvox) {
        stop('Data resolution for subject ',ii,' does not match that of first subject.')
      }
    }

    # Perform subject-level PCA and add rows to bigY
    if(verbose) cat(paste0('... performing dimension reduction\n'))
    # if (ntime < num_PCs) {
    #   warning(
    #     'For subject ',ii,
    #     ' there are fewer than `num_PCs` (',num_PCs,') time points. Skipping dimension reduction.'
    #   )
    # } #note: still use SVD to standardize variance
    # TIMER -----------------------------------------
    tick = Sys.time()
    bigY[,cols_ii] <- PCA(BOLD_ii, Q=num_PCs_ii, nV="Q")$v
    # TIMER -----------------------------------------
    # [TO DO]: if verbose?
    cat("Single subject dimension reduction timer:\n")
    print(Sys.time() - tick)
  }
  bigY <- bigY[,seq(next_col-1),drop=FALSE]

  # TIMER -----------------------------------------
  # start timer
  tick = Sys.time()

  # PERFORM GROUP-LEVEL PCA FOR DIMENSION REDUCTION
  if (verbose) { cat("Performing group-level PCA.\n") }
  bigY <- t(PCA(bigY, Q=num_ICs, nV="Q")$v)

  # PERFORM GROUP-LEVEL ICA
  if (verbose) { cat("Performing group-level ICA.\n") }
  GICA <- ica::icaimax(bigY, nc=num_ICs, center=TRUE)
  rm(bigY)
  GICA <- fMRItools::sign_flip(GICA)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # # Repeat GICA. Does not work with icaimax because there is no change, but works (badly) with fastica
  # for (nRun in 1:numRun){
  #   print(nRun)
  #   GICA <- icaimax(t(Vt_bigY), nc=num_ICs, center=TRUE)
  #   #GICA = fastICA(t(Vt_bigY), n.comp = num_ICs)
  #   icasig = GICA$S
  #   M = GICA$M # A in GIFT
  #   #M = GICA$A # for fastICA
  #   # For "averaging" multiple runs as donde in GIFT.
  #   # See icatb_calculateICA.m that can be found at
  #   # https://github.com/mandymejia/templateICA/blob/master/GIFT_GICA/GroupICATv4.0b/icatb/icatb_analysis_functions/icatb_calculateICA.m
  #   #if(nRun == 1) icasig2 = matrix(0, numRuns, size(icasig, 1), size(icasig, 2))
  #   if(nRun == 1){
  #     oldica = icasig
  #     icasig2 = matrix(data = 0, ncol = ncol(icasig), nrow = nrow(icasig))
  #     icasig2 = array(c(icasig2, 1:numRun), dim = c(nrow(icasig),ncol(icasig), numRun))
  #     M2 = matrix(data = 0, ncol = ncol(M), nrow = nrow(M))
  #     M2 = array(c(M2, 1:numRun), dim = c(nrow(M2),ncol(M2), numRun))
  #   }
  #   if(nRun > 1){
  #     # Match components
  #     # Identifies cross-correlation between current and all gicas
  #     print(sum(icasig == oldica))
  #     rho = matrix(0, ncol = ncol(icasig), nrow = ncol(icasig))
  #     for(k1 in 1:ncol(icasig)){
  #       for (k2 in 1:ncol(icasig)){
  #         rho[k1, k2] = corr_coeff(icasig[,k1], icasig2[,k2,1])
  #       }
  #     }
  #     Y = array(data=0, dim=ncol(icasig))
  #     I = Y
  #     Ys = Y
  #     for(k in 1:ncol(icasig)){
  #       Y[k] = max(abs(rho[,k]))
  #       I[k] = which.max(abs(rho[,k]))
  #       Ys[k] = rho[I(k), k] # signed correlation
  #       rho[I[k], k] = 0
  #     }
  #     # reorder and force to be positively correlated
  #     icasig = sign(rep(t(Ys), 1, nrow(icasig)))*icasig[,I]
  #     M = sign(rep(Ys, ncol(M), 1))*M[,I]
  #   }
  #   # Store icasig and M
  #   icasig2[,,nRun] = icasig
  #   M2[,,nRun] = M
  # }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # # Average Infomax runs

  # if(numRun > 1){
  #   #icasig = c(mean(icasig2))
  #   icasig = apply(icasig2, c(1,2), mean)
  #   #M = c(mean(M2))
  #   M = apply(M2, c(1,2), mean)
  #   W = pinv(M)
  # }

  # TIMER -----------------------------------------
  elapsed = Sys.time() - tick
  cat("Group-level `icaimax` timer: \n")
  # print(paste0("group-level infomax (x", numRun, ") timer:"))
  print(elapsed)

  # Format GICA as a xifti object and write out a cifti

  xii_GICA <- ciftiTools::read_cifti(
    fnames_ii[1], brainstructures=brainstructures,
    surfL_fname=surfL, surfR_fname=surfR,
  )
  xii_GICA <- ciftiTools::select_xifti(xii_GICA, rep(1, ncol(GICA$M)))
  xii_GICA <- ciftiTools::convert_xifti(xii_GICA, "dscalar")
  xii_GICA <- ciftiTools::newdata_xifti(xii_GICA, GICA$M)

  if(!is.null(out_fname)){
    return(ciftiTools::write_cifti(xii_GICA, out_fname, verbose=verbose))
  } else {
    return(xii_GICA)
  }
}
