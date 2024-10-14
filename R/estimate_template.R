#' Estimate template from DR
#'
#' Estimate variance decomposition and templates from DR estimates.
#'
#' @param DR the test/retest dual regression estimates, as an array with
#'  dimensions \eqn{M \times N \times (L \times V)}, where \eqn{M} is the number
#'  of visits (2), \eqn{N} is the number of subjects, \eqn{L} is the number of
#'  IC networks, and \eqn{V} is the number of data locations.
#'
#'  (\eqn{L} and \eqn{V} are collapsed because they are treated equivalently
#'  in the context of calculating the variance decomposition and templates).
#' @param LV A length-two integer vector giving the dimensions \eqn{L} and
#'  \eqn{V} to reshape the result. Default: \code{NULL} (do not reshape the
#'  result).
#'
#' @return List of two elements: the templates and the variance decomposition.
#'
#'  There are two version of the variance template: \code{varUB} gives the
#'  unbiased variance estimate, and \code{varNN} gives the upwardly-biased
#'  non-negative variance estimate. Values in \code{varUB} will need to be
#'  clamped above zero before using in \code{templateICA}.
#'
#' @importFrom fMRItools var_decomp
#' @export
estimate_template_from_DR <- function(
  DR, LV=NULL){

  # Check arguments.
  stopifnot(length(dim(DR)) == 3)
  nM <- dim(DR)[1]  # visits
  nN <- dim(DR)[2]  # subjects
  nLV <- dim(DR)[3] # locations & networks
  if (!is.null(LV)) {
    stopifnot(is.numeric(nLV) && all(nLV > 0) && all(nLV == round(nLV)))
    stopifnot(prod(LV) == nLV)
  }

  # Variance decomposition
  vd <- var_decomp(DR)

  # Template calculation
  # Below true for M==2. Double check correct for M > 3? (Not used currently.)
  MSB_divM <- (vd$SSB / (nN-1)) / nM
  MSE_divM <- (vd$SSR / ((nM-1)*(nN-1))) / nM
  template <- list(
    mean = vd$grand_mean,
    varUB = MSB_divM - MSE_divM,
    varNN = MSB_divM
  )

  # Format `vd`
  vd$nM <- vd$nS <- vd$grand_mean <- NULL # Get rid of redundant entries

  ## Format `template`: clamp var est above zero.
  # template$varUB <- pmax(0, template$varUB)

  # Format as matrix.
  if (!is.null(LV)) {
    template <- lapply(template, function(x){ matrix(x, nrow=LV[1], ncol=LV[2]) })
    vd <- lapply(vd, function(x){ matrix(x, nrow=LV[1], ncol=LV[2]) })
  }

  # Return
  list(template=template, var_decomp=vd)
}

#' Estimate template from DR estimates (when there are two measurements)
#'
#' Legacy version of \code{\link{estimate_template_from_DR}}
#'
#' @param DR1,DR2 the test and retest dual regression estimates (\eqn{N \times L \times V})
#'
#' @return List of two elements: the mean and variance templates
#' @keywords internal
estimate_template_from_DR_two <- function(DR1, DR2){

  # Check arguments.
  stopifnot(length(dim(DR1)) == length(dim(DR2)))
  stopifnot(all(dim(DR1) == dim(DR2)))
  N <- dim(DR1)[1]

  template <- list(mean=NULL, var=NULL)

  # Mean.
  template$mean <- t(colMeans(DR1 + DR2, na.rm=TRUE) / 2)

  # Variance.
  SSB <- 2 * colSums(((DR1 + DR2)/2 - rep(t(template$mean), each=N))^2, na.rm=TRUE)
  template$var_nn <- t(SSB / (N-1)) / 2 # MSB/2
  # Unbiased.
  # 1. Fastest method.
  var_noise <- t( (1/2) * apply(DR1 - DR2, c(2,3), var, na.rm=TRUE) )
  template$var_ub <- template$var_nn - var_noise/2

  # # 2. Previous, equivalent calculation.
  # var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
  # var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
  # var_tot <- t((var_tot1 + var_tot2)/2)
  # # noise (within-subject) variance
  # DR_diff <- DR1 - DR2;
  # var_noise <- t((1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE))
  # # signal (between-subject) variance
  # template$var <- var_tot - var_noise
  #
  # # 3. Another equivalent calculation.
  # template$var <- t(apply(
  #   abind::abind(DR1, DR2, along=1),
  #   seq(2, 3),
  #   function(q){ cov(q[seq(N)], q[seq(N+1, 2*N)], use="complete.obs") }
  # ))

  # Make negative estimates equal to zero.
  template$var_ub[template$var_ub < 0] <- 0

  template
}

#' Estimate FC template
#'
#' @param FC0 The FC estimates from \code{\link{estimate_template}}.
#' @param nu_adjust Factor by which to adjust estimate of nu.  Values < 1 will
#' inflate the template variance to avoid an over-informative prior on FC.
#' @importFrom matrixStats colVars
#' @export
estimate_template_FC <- function(FC0, nu_adjust=1){

  nL <- dim(FC0)[3]
  stopifnot(nL == dim(FC0)[4])

  FC1 <- FC0[1,,,]; FC2 <- FC0[2,,,]
  FCavg <- (FC1 + FC2)/2
  mean_FC <- apply(FCavg, c(2,3), mean, na.rm=TRUE)
  var_FC_between <- apply(FCavg, c(2,3), var, na.rm=TRUE) #this may be an overestimate but that's ok
  # mean_FC <- (colMeans(FC1, na.rm=TRUE) + colMeans(FC2, na.rm=TRUE)) / 2
  # var_FC_tot  <- (apply(FC1, c(2, 3), var, na.rm=TRUE) + apply(FC2, c(2, 3), var, na.rm=TRUE))/2
  # #var_FC_within  <- 1/2*(apply(FC1-FC2, c(2, 3), var, na.rm=TRUE))
  # var_FC_between <- var_FC_tot # var_FC_within #to avoid under-estimating the variance, given that within-subject variance in FC is often high
  # var_FC_between[var_FC_between < 0] <- NA

  nu_est <- estimate_nu(var_FC_between, mean_FC)
  nu_est <- max(nL+4, nu_est*nu_adjust)

  list(nu = nu_est,
       psi = mean_FC*(nu_est - nL - 1),
       mean_empirical = mean_FC,
       var_empirical = var_FC_between)

}

#' Cholesky-based FC sampling
#'
#' @param Chol_vals Matrix of Cholesky factorizations (upper triangular values) for all template sessions (nN*nM x nChol)
#' @param p Pivot/reordering applied to FC matrices prior to Cholesky factorization
#' @param M Number of samples to draw
#' @param chol_diag Indices of diagonal upper triangular elements
#' @param chol_offdiag Indices of off-diagonal upper triangular elements
#' @param Chol_mat_blank A nLxnL matrix indexing the upper triangular elements
#' @importFrom stats rnorm
Chol_samp_fun <- function(Chol_vals, p, M, chol_diag, chol_offdiag, Chol_mat_blank){

  nUT <- ncol(Chol_vals)
  nL <- length(chol_diag)

  logit <- function(p){log(p/(1-p))}
  expit <- function(x){exp(x)/(1+exp(x))}
  fishZ <- function(r){(1/2)*log((1+r)/(1-r))}
  fishZinv <- function(z){(exp(2*z)-1)/(exp(2*z)+1)}

  #grab upper triangle of each Cholesky matrix and logit/Fisher transform values
  Chol_vals[,1] <- NA #ignore first element since it is always 1
  Chol_vals[,chol_diag] <- logit(Chol_vals[,chol_diag]) #logit-transform diagonal values
  Chol_vals[,chol_offdiag] <- fishZ(Chol_vals[,chol_offdiag]) #z-transform off-diagonal values
  Chol_vals <- Chol_vals[,-1] #sess x (nUT-1) matrix (first element always equals 1)

  #perform SVD
  chol_svd <- svd(scale(Chol_vals, scale=FALSE)) # X = u %*% diag(d) %*% t(v)
  K <- length(chol_svd$d) #keep all the factors
  U <- chol_svd$u
  V <- chol_svd$v
  D <- diag(chol_svd$d)

  #take samples of U
  sd <- sqrt(mean(apply(U, 2, var))) #estimate SD of the scores (same across PCs)
  U_samp <- matrix(rnorm(M*K, mean=0, sd = sd), nrow=M, ncol=K)

  #2. transform back to Cholesky elememnts
  UDV_samp <- U_samp %*% D %*% t(V)
  Chol_mean <- colMeans(Chol_vals) #mean that was was removed prior to PCA
  Chol_samp <- UDV_samp + matrix(Chol_mean, nrow=M, ncol=nUT-1, byrow=TRUE) #add back in mean
  Chol_samp <- cbind(1, Chol_samp) #add back in first element (always = 1)
  Chol_samp[,chol_diag] <- expit(Chol_samp[,chol_diag]) #apply inverse-logit transformation to diagonal elements
  Chol_samp[,chol_offdiag] <- fishZinv(Chol_samp[,chol_offdiag]) #apply inverse-Fisher transformation to off-diagonal elements
  Chol_samp[,1] <- 1 #correct the first element to equal 1 again

  #3. enforce correlation scale
  for(col in 2:nL){
    inds_col <- Chol_mat_blank[1:col,col]
    #for each sample, compute the sum of squares for the current column
    SS_col <- rowSums((Chol_samp[,inds_col])^2)
    Chol_samp[,inds_col] <- Chol_samp[,inds_col] / sqrt(SS_col)
  }

  # #compute the mean and variance over samples
  # Chol_samp_mean <- UT2mat(colMeans(Chol_samp))
  # Chol_samp_var <- UT2mat(colVars(Chol_samp))
  # Chol_samp_mean[lower.tri(Chol_samp_mean)] <- NA
  # Chol_samp_var[lower.tri(Chol_samp_var)] <- NA

  #reconstruct corresponding FC matrices, compute mean and variance
  #these are provided so the user can assess how well the samples emulate the FC edge-wise pop. mean and variance
  Chol_samp_list <- apply(Chol_samp, 1, list)
  FC_samp_list <- lapply(Chol_samp_list, function(R){
    R <- R[[1]]
    R_reorder <- UT2mat(R)[,order(p)] #reverse-pivot columns of the Cholesky UT matrix
    crossprod(R_reorder) #t(R_reorder) %*% R_reorder
  })

  result <- list(Chol_samp = Chol_samp,
                 # Chol_samp_mean = Chol_samp_mean,
                 # Chol_samp_var = Chol_samp_var,
                 FC_samp_list = FC_samp_list,
                 #FC_samp_mean = FC_samp_mean,
                 #FC_samp_var = FC_samp_var,
                 chol_svd = chol_svd)
  return(result)
}

#' Transform upper-triangular elements to matrix form
#'
#' @param x Vector of upper triangular values
#' @param diag Are diagonal values included in x?  Default: \code{TRUE}.
UT2mat <- function(x, diag=TRUE){

  #determine V based on M (solution to quadratic formula since M = V*(V+1)/2 (diag=TRUE) or V*(V-1)/2 (diag=FALSE))
  nUT <- length(x)
  if(diag) {
    V <- (-1+sqrt(8*nUT+1))/2 #this is for when diag=TRUE
    if(round(V) != V) stop('Length of x must equal V(V+1)/2 for some integer V when diag=TRUE')  #this is for when diag=TRUE
  }
  if(!diag) {
    V <- (1+sqrt(8*nUT+1))/2 #this is for when diag=FALSE
    if(round(V) != V) stop('Length of x must equal V(V-1)/2 for some integer V when diag=FALSE')  #this is for when diag=TRUE
  }

  mat <- matrix(0, nrow=V, ncol=V)
  if(diag) mat[upper.tri(mat, diag=TRUE)] <- x
  if(!diag) mat[upper.tri(mat, diag=FALSE)] <- x
  return(mat)
}


# estimate_template_FC_chol <- function(FC0_chol){
#
#   # FC0_chol is 2 x n_subjects x Q*(Q+1)/2
#
#   FC_chol_avg <- (FC0_chol[1,,] + FC0_chol[2,,])/2
#   mean_FC_chol <- apply(FC_chol_avg, 2, mean, na.rm=TRUE)
#   var_FC_chol1 <- apply(FC0_chol[1,,], 2, var, na.rm=TRUE)
#   var_FC_chol2 <- apply(FC0_chol[2,,], 2, var, na.rm=TRUE)
#   var_FC_chol <- (var_FC_chol1 + var_FC_chol2)/2
#
#   list(mean_empirical = mean_FC_chol,
#        var_empirical = var_FC_chol)
#
# }

#' Estimate template
#'
#' Estimate template for Template ICA based on fMRI data
#'
#' All fMRI data (entries in \code{BOLD} and \code{BOLD2}, and \code{GICA}) must
#'  be in the same spatial resolution.
#'
#' @param BOLD,BOLD2 Vector of subject-level fMRI data in one of the following
#'  formats: CIFTI file paths, \code{"xifti"} objects, GIFTI file paths,
#'  \code{"gifti"} objects, NIFTI file paths, \code{"nifti"} objects,
#'  or \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data
#'  locations and \eqn{T} is the number of timepoints.
#
#  If GIFTI or \code{"gifti"}, the input can also be a length two list,
#  where the first list element is a length \eqn{N} vector for the left
#  hemisphere and the second list element is a length \eqn{N} vector for the
#  right hemisphere.
#'
#'  If \code{BOLD2} is provided it must be in the same format as \code{BOLD};
#'  \code{BOLD} will be the test data and \code{BOLD2} will be the retest data.
#'  \code{BOLD2} should be the same length as \code{BOLD} and have the same
#'  subjects in the same order. If \code{BOLD2} is not provided, \code{BOLD}
#'  will be split in half; the first half will be the test data and the second
#'  half will be the retest data.
#' @param GICA Group ICA maps in a format compatible with \code{BOLD}. Can also
#'  be a (vectorized) numeric matrix (\eqn{V \times Q}) no matter the format of
#'  \code{BOLD}. Its columns will be centered.
#'
#'  New: can also be a parcellation in CIFTI format (other formats to be
#'  implemented in the future). The parcellation should have the same locations
#'  and one column, with integer values indicating the parcel to which each
#'  location belongs to. Each parcel is modeled as a brain map; instead of the
#'  first step of dual regression, the medial timecourse of each parcel is used.
#'
#' @param mask Required if the entries of \code{BOLD} are NIFTI
#'  file paths or \code{"nifti"} objects, optional for other formats. For
#'  NIFTI, this is a brain map formatted as a binary array of the same spatial
#'  dimensions as the fMRI data, with \code{TRUE} corresponding to in-mask
#'  voxels. For other formats, a logical vector.
#' @param inds Numeric indices of the group ICs to include in the template. If
#'  \code{NULL}, use all group ICs (default).
#'
#'  If \code{inds} is provided, the ICs not included will be removed after
#'  calculating dual regression, not before. This is because removing the ICs
#'  prior to dual regression would leave unmodelled signals in the data, which
#'  could bias the templates.
#' @inheritParams scale_Param
#' @param scale_sm_surfL,scale_sm_surfR,scale_sm_FWHM Only applies if
#'  \code{scale=="local"} and \code{BOLD} represents surface data (CIFTI or
#'  GIFTI). To smooth the standard deviation estimates used for local scaling,
#'  provide the surface geometries along which to smooth as GIFTI geometry files
#'  or \code{"surf"} objects, as well as the smoothing FWHM (default: \code{2}).
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
#' @inheritParams TR_param
#' @inheritParams hpf_param
#' @inheritParams GSR_Param
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI
#'  file paths. This is a character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.
#' @param resamp_res Only applies if the entries of \code{BOLD} are CIFTI file paths.
#'  Resample the data upon reading it in? Default: \code{NULL} (no resampling).
#' @param mask Required if \code{BOLD} are NIFTI file paths or \code{"nifti"} 
#'  objects, and optional for other formats. For NIFTI data, this is a brain map
#'  formatted as a logical array of the same spatial dimensions as the fMRI 
#'  data, with \code{TRUE} corresponding to in-mask voxels. For other data, this
#'  is a logical vector with the same length as the number of locations in
#'  \code{GICA}, with \code{TRUE} corresponding to in-mask locations.
#' @param keep_S Keep the DR estimates of S? If \code{FALSE} (default), do not save
#'  the DR estimates and only return the templates. If \code{TRUE}, the DR
#'  estimates of S are returned too. If a single file path, save the DR estimates as
#'  an RDS file at that location rather than returning them.
#   [TO DO] If a list of two vectors of file paths with the same lengths as
#   \code{BOLD}, save the DR estimates as individual files at these locations in
#   the appropriate format (CIFTI, NIFTI, or RDS files, depending on \code{BOLD}).
#' @param keep_FC Keep the DR estimates of the FC cor(A)? If \code{FALSE} (default), do not save
#'  the DR estimates and only return the templates. If \code{TRUE}, the DR
#'  estimates of cor(A) and its Cholesky factor are returned too.
#'  If a single file path, save the DR estimates as
#'  an RDS file at that location rather than returning them.
#   [TO DO] If a list of two vectors of file paths with the same lengths as
#   \code{BOLD}, save the DR estimates as individual files at these locations in
#   the appropriate format (CIFTI, NIFTI, or RDS files, depending on \code{BOLD}).
#' @param Q2,Q2_max Obtain dual regression estimates after denoising? Denoising is
#'  based on modeling and removing nuisance ICs. It may result in a cleaner
#'  estimate for smaller datasets, but it may be unnecessary (and time-consuming)
#'  for larger datasets.
#'
#'  Set \code{Q2} to control denoising: use a positive integer to specify the
#'  number of nuisance ICs, \code{NULL} to have the number of nuisance ICs
#'  estimated by PESEL, or zero (default) to skip denoising.
#'
#'  If \code{is.null(Q2)}, use \code{Q2_max} to specify the maximum number of
#'  nuisance ICs that should be estimated by PESEL. \code{Q2_max} must be less
#'  than \eqn{T * .75 - Q} where \eqn{T} is the minimum number of timepoints in
#'  each fMRI scan and \eqn{Q} is the number of group ICs. If \code{NULL}
#'  (default), \code{Q2_max} will be set to \eqn{T * .50 - Q}, rounded.
#' @param FC Include the functional connectivity template? Default: \code{TRUE}.
#' @param FC_nPivots Number of pivots to use in Cholesky-based FC template
#' estimation.  Set to zero to skip Cholesky-based FC template estimation. Default: 100.
#' @param FC_nSamp Number of FC matrix samples to generate across all pivots. This
#' should be a multiple of FC_nPivots.
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{1e-6}. Variance is calculated on the original data, before
#'  any normalization.
#' @param maskTol For computing the dual regression results for each subject:
#'  tolerance for number of locations masked out due to low
#'  variance or missing values. If more than this many locations are masked out,
#'  a subject is skipped without calculating dual regression. \code{maskTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of locations can be masked out.
#'
#'  If \code{BOLD2} is provided, masks are calculated for both scans and then
#'  the intersection of the masks is used, for each subject.
#' @param missingTol For computing the variance decomposition across all subjects:
#'  tolerance for number of subjects masked out due to low variance or missing
#'  values at a given location. If more than this many subjects are masked out,
#'  the location's value will be \code{NA} in the templates. \code{missingTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of subjects can be masked out
#'  at a given location.
#' @param usePar,wb_path Parallelize the DR computations over subjects? Default:
#'  \code{FALSE}. Can be the number of cores to use or \code{TRUE}, which will
#'  use the number on the PC minus two. If the input data is in CIFTI format, the
#'  \code{wb_path} must also be provided.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#'
#' @importFrom stats cov quantile
#' @importFrom fMRItools is_1 is_integer is_posNum colCenter unmask_mat infer_format_ifti_vec all_binary
#' @importFrom abind abind
#'
#' @return A list: the \code{template} and \code{var_decomp} with entries in
#'  matrix format; the \code{mask} of locations without template values due to
#'  too many low variance or missing values; the function \code{params} such as
#'  the type of scaling and detrending performed; the \code{dat_struct} which can be
#'  used to convert \code{template} and \code{var_decomp} to \code{"xifti"} or
#'  \code{"nifti"} objects if the \code{BOLD} format was CIFTI or NIFTI data;
#'  and DR results if \code{isTRUE(keep_S)} and/or \code{isTRUE(keep_FC)}.
#'
#'  Use \code{summary} to print a description of the template results, and
#'  for CIFTI-format data use \code{plot} to plot the template mean and variance
#'  estimates. Use \code{\link{export_template}} to save the templates to
#'  individual RDS, CIFTI, or NIFTI files (depending on the \code{BOLD} format).
#' @export
#'
#' @examples
#' nT <- 21
#' nV <- 140
#' nQ <- 6
#' mU <- matrix(rnorm(nV*nQ), nrow=nV)
#' mS <- mU %*% diag(seq(nQ, 1)) %*% matrix(rnorm(nQ*nT), nrow=nQ)
#' BOLD <- list(B1=mS, B2=mS, B3=mS)
#' BOLD <- lapply(BOLD, function(x){x + rnorm(nV*nT, sd=.05)})
#' GICA <- mU
#' estimate_template(BOLD=BOLD, GICA=mU, FC_nSamp=2000)
#'
#' \dontrun{
#'  estimate_template(
#'    run1_cifti_fnames, run2_cifti_fnames,
#'    gICA_cifti_fname, brainstructures="all",
#'    scale="global", TR=0.71, Q2=NULL, varTol=10
#'  )
#' }
estimate_template <- function(
  BOLD, BOLD2=NULL,
  GICA,
  mask=NULL,
  inds=NULL,
  scale=c("local", "global", "none"),
  scale_sm_surfL=NULL,
  scale_sm_surfR=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01,
  GSR=FALSE,
  Q2=0, Q2_max=NULL,
  brainstructures="all", resamp_res=NULL,
  keep_S=FALSE, keep_FC=FALSE,
  FC=TRUE,
  FC_nPivots=100,
  FC_nSamp=50000,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  usePar=FALSE, wb_path=NULL,
  verbose=TRUE) {

  # Check arguments ------------------------------------------------------------

  # Simple argument checks.
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("local", "global", "none"))
  stopifnot(fMRItools::is_1(scale_sm_FWHM, "numeric"))
  if (is.null(hpf)) { hpf <- 0 }
  if (is.null(TR)) {
    if (hpf==.01) {
      message("Setting `hpf=0` because `TR` was not provided. Either provide `TR` or set `hpf=0` to disable this message.")
      hpf <- 0
    } else if (hpf!=0) {
      stop("Cannot apply `hpf` because `TR` was not provided. Either provide `TR` or set `hpf=0`.")
    }
  } else {
    stopifnot(fMRItools::is_posNum(TR))
    stopifnot(fMRItools::is_posNum(hpf, zero_ok=TRUE))
  }
  stopifnot(fMRItools::is_1(GSR, "logical"))
  if (!is.null(Q2)) { # Q2_max checked later.
    stopifnot(fMRItools::is_integer(Q2) && (Q2 >= 0))
  }
  stopifnot(fMRItools::is_1(FC, "logical"))
  stopifnot(fMRItools::is_1(varTol, "numeric"))
  if (varTol < 0) { cat("Setting `varTol=0`."); varTol <- 0 }
  stopifnot(fMRItools::is_posNum(maskTol, zero_ok=TRUE))
  stopifnot(fMRItools::is_posNum(missingTol, zero_ok=TRUE))
  stopifnot(fMRItools::is_1(verbose, "logical"))
  real_retest <- !is.null(BOLD2)

  # `keep_S`
  if (is.logical(keep_S)) {
    stopifnot(length(keep_S)==1)
  } else {
    if (is.character(keep_S)) {
      stopifnot(length(keep_S)==1)
      if (!dir.exists(dirname(keep_S))) { stop('Directory part of `keep_S` does not exist.') }
      if (!endsWith(keep_S, ".rds")) { keep_S <- paste0(keep_S, ".rds") }
    } else if (is.list(keep_S)) {
      stop("Not supported: `keep_S` must be `TRUE`, `FALSE`, or a single file path.")
      # [TO DO]
      # if (length(keep_S) != 2) {
      #   stop("If `keep_S` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      # }
      # if (length(keep_S[[1]]) != nN || length(keep_S[[2]]) != nN) {
      #   stop("If `keep_S` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      # }
      # if (!all(dir.exists(dirname(do.call(c, keep_S))))) { stop('At least one directory part of `keep_S` does not exist.') }
    }
  }

  # `keep_FC`
  if (is.logical(keep_FC)) {
    stopifnot(length(keep_FC)==1)
  } else {
    if (is.character(keep_FC)) {
      stopifnot(length(keep_FC)==1)
      if (!dir.exists(dirname(keep_FC))) { stop('Directory part of `keep_FC` does not exist.') }
      if (!endsWith(keep_FC, ".rds")) { keep_FC <- paste0(keep_FC, ".rds") }
    } else if (is.list(keep_FC)) {
      stop("Not supported: `keep_FC` must be `TRUE`, `FALSE`, or a single file path.")
      # [TO DO]
      # if (length(keep_FC) != 2) {
      #   stop("If `keep_FC` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      # }
      # if (length(keep_FC[[1]]) != nN || length(keep_FC[[2]]) != nN) {
      #   stop("If `keep_FC` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      # }
      # if (!all(dir.exists(dirname(do.call(c, keep_FC))))) { stop('At least one directory part of `keep_FC` does not exist.') }
    }
  }

  # `usePar`
  if (!isFALSE(usePar)) {

    check_parallel_packages()
    cores <- parallel::detectCores()
    if (isTRUE(usePar)) { nCores <- cores[1] - 2 } else { nCores <- usePar; usePar <- TRUE }
    if (nCores < 2) {
      usePar <- FALSE
    } else {
      cluster <- parallel::makeCluster(nCores, outfile="")
      doParallel::registerDoParallel(cluster)
    }
  }

  # `BOLD` and `BOLD2` ---------------------------------------------------------
  # Determine the format of `BOLD` and `BOLD2`.
  format <- infer_format_ifti_vec(BOLD)[1]
  if (real_retest) {
    format2 <- infer_format_ifti_vec(BOLD2)[1]
    if (format2 != format) {
      stop("`BOLD` format is ", format, ", but `BOLD2` format is ", format2, ".")
    }
  }
  FORMAT <- get_FORMAT(format)
  FORMAT_extn <- switch(FORMAT,
    CIFTI=".dscalar.nii",
    GIFTI=".func.gii",
    NIFTI=".nii",
    MATRIX=".rds"
  )

  check_req_ifti_pkg(FORMAT)

  # Ensure `BOLD2` is the same length.
  if (real_retest) {
    if (length(BOLD) != length(BOLD2)) {
      stop("`BOLD2` represents corresponding retest data for `BOLD`, so it must have the same length as `BOLD`.")
    }
  }

  # If BOLD (and BOLD2) is a CIFTI, GIFTI, NIFTI, or RDS file, check that the file paths exist.
  if (format %in% c("CIFTI", "GIFTI", "NIFTI", "RDS")) {
    missing_BOLD <- !file.exists(BOLD)
    if (all(missing_BOLD)) stop('The files in `BOLD` do not exist.')
    if (real_retest) {
      missing_BOLD2 <- !file.exists(BOLD2)
      if (all(missing_BOLD2)) stop('The files in `BOLD2` do not exist.')
      # determine pairs with at least one missing scan
      missing_BOLD <- missing_BOLD | missing_BOLD2
      rm(missing_BOLD2)
      if (all(missing_BOLD)) stop('Files in `BOLD` and/or `BOLD2` are missing such that no complete pair of data exists.')
    }
    if (any(missing_BOLD)) {
      if (real_retest) {
        warning('There are ', sum(missing_BOLD), ' pairs of `BOLD` and `BOLD2` with at least one non-existent scan. These pairs will be excluded from template estimation.')
        BOLD <- BOLD[!missing_BOLD]
        BOLD2 <- BOLD[!missing_BOLD]
      } else {
        warning('There are ', sum(missing_BOLD), ' scans in `BOLD` that do not exist. These scans will be excluded from template estimation.')
        BOLD <- BOLD[!missing_BOLD]
      }
    }
  }
  nN <- length(BOLD)

  # Check `scale_sm_FWHM`
  if (scale_sm_FWHM !=0 && FORMAT %in% c("NIFTI", "MATRIX")) {
    if (scale_sm_FWHM==2) {
      cat("Setting `scale_sm_FWHM == 0`.\n")
    } else {
      if (FORMAT == "NIFTI") {
        # [TO DO] make this available
        warning( "Setting `scale_sm_FWHM == 0` (Scale smoothing not available for volumetric data.).\n")
      } else {
        warning( "Setting `scale_sm_FWHM == 0` (Scale smoothing not available for data matrices: use CIFTI/GIFTI files.).\n")
      }
    }
    scale_sm_FWHM <- 0
  }

  # [TO DO]: Mysteriously, RDS format is not working with parallel.
  if (usePar && format=="RDS") { stop("Parallel computation not working with RDS file input. Working on this!") }

  # `GICA` ---------------------------------------------------------------------
  # Convert `GICA` to a numeric data matrix or array.
  GICA_parc <- FALSE; GICA_parc_table <- NULL
  if (FORMAT == "CIFTI") {
    if (is.character(GICA)) { GICA <- ciftiTools::read_cifti(GICA, brainstructures=brainstructures, resamp_res=resamp_res) }
    if (ciftiTools::is.xifti(GICA, messages=FALSE)) {
      #for parcellation (instead of ICA)
      if ((ncol(GICA) == 1) && (all(as.matrix(GICA) == round(as.matrix(GICA))))) {
        GICA_parc <- TRUE
        GICA_parc_vals <- sort(unique(c(as.matrix(GICA))))
        GICA_parc_table <- GICA$meta$cifti$labels[[1]]
        if (is.null(GICA_parc_table)) {
          GICA_parc_table <- data.frame(
            Key = GICA_parc_vals,
            Red = 1,
            Green = 1,
            Blue = 1,
            Alpha = 1,
            row.names = paste("ParcelKey", GICA_parc_vals)
          )
        } else {
          # Drop empty levels
          GICA_parc_table <- GICA_parc_table[GICA_parc_table$Key %in% GICA_parc_vals,]
          stopifnot(nrow(GICA_parc_table) > 1)
        }
      }
      xii1 <- ciftiTools::select_xifti(GICA, 1) # for formatting output
      GICA <- as.matrix(GICA)
    } else {
      # Get `xii1` from first data entry.
      xii1 <- BOLD[[1]]
      if (is.character(xii1)) {
        xii1 <- ciftiTools::read_cifti(xii1, brainstructures=brainstructures, resamp_res=resamp_res, idx=1)
      }
      xii1 <- ciftiTools::convert_xifti(ciftiTools::select_xifti(xii1, 1), "dscalar")
    }
    stopifnot(is.matrix(GICA))
  } else if (FORMAT == "GIFTI") {
    if (is.character(GICA)) { GICA <- gifti::readgii(GICA) }
    ghemi <- GICA$file_meta["AnatomicalStructurePrimary"]
    if (!(ghemi %in% c("CortexLeft", "CortexRight"))) {
      stop("AnatomicalStructurePrimary metadata missing or invalid for GICA.")
    }
    ghemi <- switch(ghemi, CortexLeft="left", CortexRight="right")
    GICA <- do.call(cbind, GICA$data)
  } else if (FORMAT == "NIFTI") {
    if (is.character(GICA)) { GICA <- RNifti::readNifti(GICA) }
    stopifnot(length(dim(GICA)) > 1)
  } else if (FORMAT == "MATRIX") {
    if (is.character(GICA)) { GICA <- readRDS(GICA) }
    stopifnot(is.matrix(GICA))
  }

  # `inds`.
  if (GICA_parc) {
    nQ <- nrow(GICA_parc_table)
    if (!is.null(inds)) {
      if (!all(inds %in% GICA_parc_table$Key)) stop('Invalid entries in inds argument.')
      nL <- length(inds)
    } else {
      inds <- GICA_parc_table$Key
      nL <- nQ
    }
    inds <- match(inds, GICA_parc_table$Key)
  } else {
    nQ <- dim(GICA)[length(dim(GICA))]
    if (!is.null(inds)) {
      if (!all(inds %in% seq(nQ))) stop('Invalid entries in inds argument.')
      nL <- length(inds)
    } else {
      inds <- seq(nQ)
      nL <- nQ
    }
  }

  # [TO DO]: NA in GICA?

  # [TO DO]: Check that FC_nPivots is a positive integer or is equal to zero.

  # [TO DO]: Check that FC_nSamp is a multiple of FC_nPivots

  # `mask` ---------------------------------------------------------------------
  # Get `mask` as a logical array (NIFTI) or vector (everything else).
  # For NIFTI, append NIFTI header from GICA to `mask`.
  # Apply mask to `GICA`, and if NIFTI, vectorize `GICA`.
  if (FORMAT == "NIFTI") {
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- RNifti::readNifti(mask); mask <- array(as.logical(mask), dim=dim(mask)) }
    if (dim(mask)[length(dim(mask))] == 1) { mask <- array(mask, dim=dim(mask)[length(dim(mask))-1]) }
    if (is.numeric(mask)) {
      cat("Coercing `mask` to a logical array.\n")
      if (!fMRItools::all_binary(mask)) {
        cat("Warning: values other than 0 or 1 in mask.\n")
      }
      mask <- array(as.logical(mask), dim=dim(mask))
    }
    nI <- dim(mask); nV <- sum(mask)
    stopifnot(length(dim(GICA)) %in% c(2, length(nI)+1))
    if (length(dim(GICA)) == length(nI)+1) {
      if (length(dim(GICA)) != 2) {
        stopifnot(all(dim(GICA)[seq(length(dim(GICA))-1)] == nI))
      }
      # Append NIFTI header.
      mask <- RNifti::asNifti(array(mask, dim=c(dim(mask), 1)), reference=GICA)
      # Vectorize `GICA`.
      if (all(dim(GICA)[seq(length(dim(GICA))-1)] == nI)) {
        GICA <- if (GICA_parc) {
          matrix(GICA[as.logical(mask)], ncol=1)
        } else {
          matrix(GICA[rep(as.logical(mask), nQ)], ncol=nQ)
        }
        stopifnot(nrow(GICA) == nV)
      }
    }
  } else { #For non-NIFTI data, mask is not required but can be provided
    if (!is.null(mask)) {
      if (FORMAT == "GIFTI") {
        mask <- as.logical(do.call(cbind, gifti::read_gifti(mask)$data))
      } else if (FORMAT == "GIFTI2") {
        mask <- as.list(mask)
        if (length(mask) != 2) {
          stop("`mask` should be a length-2 list of GIFTI data: left hemisphere first, right hemisphere second.")
        }
        for (ii in 1:2) {
          mask[[ii]] <-  as.logical(do.call(cbind, gifti::read_gifti(mask[[ii]])$data))
        }
        mask <- do.call(c, mask)
      } else if (FORMAT == "MATRIX") {
        stopifnot(is.vector(mask))
        stopifnot(is.numeric(mask) || is.logical(mask))
      }
      if (is.numeric(mask)) {
        cat("Coercing `mask` to a logical vector.\n")
        if (!all_binary(mask)) {
          cat("Warning: values other than 0 or 1 in mask.\n")
        }
        mask <- as.logical(mask)
      }
      nI <- length(mask); nV <- sum(mask)
      # Check `GICA` and `mask` dimensions match.
      stopifnot(nrow(GICA) == nI)
      # Apply mask to GICA.
      GICA <- GICA[mask,,drop=FALSE]
    } else { #end if(!is.null(mask))
      nI <- nV <- nrow(GICA)
    }
  } #end else (not NIFTI format)

  # Center `GICA` columns.
  center_Gcols <- TRUE
  if (center_Gcols && (!GICA_parc)) { GICA <- fMRItools::colCenter(GICA) }

  # Print summary of data ------------------------------------------------------
  format2 <- if (format == "data") { "numeric matrix" } else { format }
  if (verbose) {
    cat("Data input format:             ", format2, "\n")
    cat("Image dimensions:              ", paste(nI, collapse=" x "), "\n")
    cat("Initial in-mask locations:     ", nV, "\n")
    if (FORMAT == "GIFTI") {
      cat("Cortex Hemisphere:             ", ghemi, "\n")
    }
    if (GICA_parc) {
      cat('Number of original parcels:    ', nQ, "\n")
      cat('Number of template parcels:    ', nL, "\n")
    } else {
      cat('Number of original group ICs:  ', nQ, "\n")
      cat('Number of template ICs:        ', nL, "\n")
    }
    cat('Number of training subjects:   ', nN, "\n")
    if (FC) { cat('\nIncluding Cholesky-based FC template with ',
      FC_nPivots,' random pivots.\n') }
  }

  # Process each scan ----------------------------------------------------------
  nM <- 2

  # Initialize Cholesky pivots for Chol-based FC template ---------------------
  if (FC) {
    if(FC_nPivots > 0){
      pivots <- Chol_samp <- Chol_svd <- vector('list', length=FC_nPivots)
      FC_nSamp2 <- round(FC_nSamp/FC_nPivots) #number of samples per pivot
      for(pp in 1:FC_nPivots){
        pivots[[pp]] <- sample(1:nL, nL, replace = FALSE)
      }
    }
  } #end setup for FC template estimation


  if (usePar) {
    check_parallel_packages()

    # Loop over subjects.
    `%dopar%` <- foreach::`%dopar%`
    q <- foreach::foreach(ii = seq(nN), .packages=c("ciftiTools", "fMRItools", "templateICAr")) %dopar% {
      if (FORMAT=="CIFTI" || FORMAT=="GIFTI") {
        # Load the workbench.
        if (is.null(wb_path)) {
          stop("`wb_path` is required for parallel computation.")
        }
        requireNamespace("ciftiTools")
        ciftiTools::ciftiTools.setOption("wb_path", wb_path)
      }

      # Initialize output.
      out <- list(DR=array(NA, dim=c(nM, 1, nL, nV)))
      if (FC) {
        out$FC <- array(NA, dim=c(nM, 1, nL, nL))
       } #end setup for FC template estimation
      out$sigma_sq <- array(NA, dim=c(nM, 1, nV))

      # Dual regression.
      if(verbose) { cat(paste0(
        '\nSubject ', ii,' of ', nN, ".\n"
      )) }
      if (real_retest) { B2 <- BOLD2[[ii]] } else { B2 <- NULL }
      DR_ii <- try(dual_reg2(
        BOLD[[ii]], BOLD2=B2,
        format=format,
        GICA=GICA, GICA_parc_table=GICA_parc_table,
        mask=mask,
        keepA=FC,
        GSR=GSR,
        scale=scale,
        scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
        scale_sm_FWHM=scale_sm_FWHM,
        TR=TR, hpf=hpf,
        Q2=Q2, Q2_max=Q2_max,
        brainstructures=brainstructures, resamp_res=resamp_res,
        varTol=varTol, maskTol=maskTol,
        verbose=verbose
      ))

      # Add results if this subject was not skipped.
      # (Subjects are skipped if too many locations are masked out.)
      if (is.null(DR_ii)) {
        if(verbose) { cat(paste0(
          '\tSubject ', ii,' was skipped (too many masked locations).\n'
        )) }
      } else if (inherits(DR_ii, "try-error")) {
        cat(paste0(
          '\tSubject ', ii,' was skipped (error).\n'
        ))
        if (ii==1) { message(DR_ii); stop("Error on first subject. Check data?") }
      } else {
        out$DR[1,,,] <- DR_ii$test$S[inds,]
        out$DR[2,,,] <- DR_ii$retest$S[inds,]
        if(FC) {
          out$FC[1,,,] <- cov(DR_ii$test$A[,inds])
          out$FC[2,,,] <- cov(DR_ii$retest$A[,inds])
          #out$FC_chol[1,,] <- chol(out$FC[1,,,])[upper.tri(out$FC[1,,,], diag=TRUE)]
          #out$FC_chol[2,,] <- chol(out$FC[2,,,])[upper.tri(out$FC[2,,,], diag=TRUE)]
        }
        #save residual variance for rescaling
        out$sigma_sq[1,,] <- DR_ii$test$sigma_sq
        out$sigma_sq[2,,] <- DR_ii$retest$sigma_sq
      }
      out
    }

    # Aggregate.
    DR0 <- abind::abind(lapply(q, `[[`, "DR"), along=2)
    if (FC) {
      FC0 <- abind::abind(lapply(q, `[[`, "FC"), along=2)
      #FC0_chol <- abind::abind(lapply(q, `[[`, "FC_chol"), along=2)
    }
    sigma_sq0 <- abind::abind(lapply(q, `[[`, "sigma_sq"), along=2)


    doParallel::stopImplicitCluster()

  } else {
    # Initialize output.
    DR0 <- array(NA, dim=c(nM, nN, nL, nV)) # measurements by subjects by components by locations
    if(FC) {
      FC0 <- array(NA, dim=c(nM, nN, nL, nL)) # for functional connectivity template
      #FC0_chol <- array(NA, dim=c(nM, nN, nL*(nL+1)/2))
    }
    sigma_sq0 <- array(NA, dim=c(nM, nN, nV))

    for (ii in seq(nN)) {
      if(verbose) { cat(paste0(
        '\nSubject ', ii,' of ', nN, '.\n'
      )) }
      if (real_retest) { B2 <- BOLD2[[ii]] } else { B2 <- NULL }

      DR_ii <- try(dual_reg2(
        BOLD[[ii]], BOLD2=B2,
        format=format,
        GICA=GICA, GICA_parc_table=GICA_parc_table,
        mask=mask,
        keepA=FC,
        GSR=GSR,
        scale=scale,
        scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
        scale_sm_FWHM=scale_sm_FWHM,
        TR=TR, hpf=hpf,
        Q2=Q2, Q2_max=Q2_max,
        brainstructures=brainstructures, resamp_res=resamp_res,
        varTol=varTol, maskTol=maskTol,
        verbose=verbose
      ))

      # Add results if this subject was not skipped.
      # (Subjects are skipped if error or too many locations are masked out.)
      if (is.null(DR_ii)) {
        if(verbose) { cat(paste0(
          '\tSubject ', ii,' was skipped (too many masked locations).\n'
        )) }
      } else if (inherits(DR_ii, "try-error")) {
        cat(paste0(
          '\tSubject ', ii,' was skipped (error).\n'
        ))
        if (ii==1) { message(DR_ii); stop("Error on first subject. Check data?") }
      } else {
        DR0[1,ii,,] <- DR_ii$test$S[inds,]
        DR0[2,ii,,] <- DR_ii$retest$S[inds,]
        if(FC) {
          FC0[1,ii,,] <- cov(DR_ii$test$A[,inds])
          FC0[2,ii,,] <- cov(DR_ii$retest$A[,inds])
          if(!all(round(diag(FC0[1,ii,,]),6) == 1)) stop('var(A) should be 1 but it is not')
          if(!all(round(diag(FC0[2,ii,,]),6) == 1)) stop('var(A) should be 1 but it is not')
          #FC0_chol[1,ii,] <- chol(FC0[1,ii,,])[upper.tri(FC0[1,ii,,], diag=TRUE)]
          #FC0_chol[2,ii,] <- chol(FC0[2,ii,,])[upper.tri(FC0[2,ii,,], diag=TRUE)]
        }
        sigma_sq0[1,ii,] <- DR_ii$test$sigma_sq
        sigma_sq0[2,ii,] <- DR_ii$retest$sigma_sq
      }
    }
    # Conserve memory.
    rm(DR_ii)
  }

  # Aggregate results, and compute templates. ----------------------------------

  # Mask out locations for which too many subjects do not have data
  #   because of missing values or low variance.
  nVm <- nV # `nVm` will be the number of locations after masking with `mask2`.
  if (missingTol < 1) { missingTol <- missingTol * nN }
  # Assume is.na(DR0[mm,,,]) is the same for all mm
  # Assume is.na(DR0[,,cc,]) is the same for all cc
  NA_counts <- colSums(is.na(DR0[1,,1,]))
  mask2 <- NA_counts < missingTol
  use_mask2 <- !all(mask2)
  if (use_mask2) {
    if (all(!mask2)) { stop(
      "No locations meet the minimum number of subjects with data (",
      round(missingTol, digits=3), "). Check `maskTol` and `missingTol`."
    ) }
    DR0 <- DR0[,,,mask2,drop=FALSE]
    nVm <- sum(mask2)
    sigma_sq0 <- sigma_sq0[,,mask2]
  }
  # Note that `NA` values may still exist in `DR0`.

  #Average `sigma_sq0` across subjects and test/retest.
  sigma_sq0 <- colMeans(sigma_sq0, dims=2, na.rm=TRUE)

  # Vectorize components and locations
  DR0 <- array(DR0, dim=c(nM, nN, nL*nVm))

  if (verbose) { cat("\nCalculating spatial IC template.\n") }
  # Estimate the mean and variance templates.
  # Also obtain the variance decomposition.
  x <- estimate_template_from_DR(DR0, c(nL, nVm))
  template <- lapply(x$template, t)
  var_decomp <- lapply(x$var_decomp, t)
  rm(x)

  #rescale mean and variance of S to standardize residual var
  #rescaling residuals by \sigma_v implies that s_v rescaled in the same way
  #this is the same as rescaling mean(s_v) and var(s_v)
  rescale <- sqrt(sigma_sq0)/mean(sqrt(sigma_sq0), na.rm=TRUE)
  rescale <- matrix(rescale, nrow=length(sigma_sq0), ncol=nL)
  template$mean <- template$mean / rescale #scale mean(S)
  template[2:3] <- lapply(template[2:3], function(x) return(x / (rescale^2)) ) #scale var(S)
  var_decomp <- lapply(var_decomp, function(x) return(x / (rescale^2) ) ) #scale var(S)

  # Unmask the data matrices (relative to `mask2`, not `mask`).
  if (use_mask2) {
    template <- lapply(template, fMRItools::unmask_mat, mask=mask2)
    var_decomp <- lapply(var_decomp, fMRItools::unmask_mat, mask=mask2)
  }



  # Estimate FC template
  if(FC){

    if (verbose) { cat("\nCalculating parametric FC template.\n") }

    template$FC <- estimate_template_FC(FC0) #estimate IW parameters

    #for Cholesky-based FC template
    FC_samp_list <- NULL #collect FC samples to compute mean/var at end
    if(FC_nPivots > 0){

      if (verbose) { cat("\nGenerating samples for Cholesky-based FC template.\n") }

      #do Cholesky-based FC template sampling
      FC_samp_logdet <- NULL
      FC_samp_cholinv <- vector('list', length = FC_nPivots)

      #setup
      nUT <- nL*(nL+1)/2 #number of upper triangular elements
      Chol_mat_blank <- matrix(0, nL, nL)
      Chol_mat_blank[upper.tri(Chol_mat_blank, diag=TRUE)] <- 1:nUT #indices of UT elements
      chol_diag <- diag(Chol_mat_blank) #which Cholesky UT elements are on the diagonal
      chol_offdiag <- Chol_mat_blank[upper.tri(Chol_mat_blank, diag=FALSE)] #which Cholesky UT elements are on the diagonal

      #for each pivot: (steps 1-4 performed by Chol_samp function)
      #1. perform Cholesky on each session, transform values to real line
      #2. vectorize and perform SVD (save D, V and sd(U) for each pivot to facilitate additional sampling)
      #3. draw FC_nSamp2 samples, construct Cholesky elements, and reverse transform
      #4. rescale values to correspond to a correlation matrix
      #5. for each sample, compute log determinant and inverse
          # a) log|X| = 2*sum(log(diag(L))), where L is the upper or lower triangular Cholesky matrix
          # b) a'X^(-1)a can be written b'b, where b = R_p^(-T)*P*a, where R_p is the upper triangular Cholesky matrix
      for(pp in 1:FC_nPivots){

        #perform Cholesky decomposition on each matrix in FC0 --> dim(FC0) = c(nM, nN, nL, nL)
        Chol_p <- apply(FC0, 1:2, function(x, pivot){ # dim = nChol x nM x nN
          xp <- x[pivot,pivot]
          chol_xp <- chol(xp)
          chol_xp[upper.tri(chol_xp, diag = TRUE)]
        }, pivot = pivots[[pp]])
        #rbind across sessions to form a matrix
        Chol_mat_p <- rbind(t(Chol_p[,1,]), t(Chol_p[,2,])) # dim = nM*nN x nChol

        #take samples
        Chol_samp_pp <- Chol_samp_fun(Chol_mat_p, p=pivots[[pp]], M=FC_nSamp2,
                                 chol_diag, chol_offdiag, Chol_mat_blank) #returns a nSamp2 x nChol matrix
        Chol_samp[[pp]] <- Chol_samp_pp$Chol_samp
        Chol_svd[[pp]] <- Chol_samp_pp$chol_svd
        FC_samp_list <- c(FC_samp_list, Chol_samp_pp$FC_samp_list)

        # 5(a) compute the log(det(FC)) for each pivoted Cholesky sample
        logdet_p <- 2 * rowSums(log(Chol_samp[[pp]][,chol_diag])) #log(det(X)) = 2*sum(log(diag(R)))
        FC_samp_logdet <- c(FC_samp_logdet, logdet_p)

        # 5(b) compute the inverse of each pivoted Cholesky sample
        # If chol_inv is the inverse of the pivoted UT cholesky matrix,
        # then chol_inv[op,] %*% t(chol_inv[op,]) gives inv(FC), where op = order(pivot)
        op <- order(pivots[[pp]])
        FC_samp_cholinv[[pp]] <- t(apply(Chol_samp[[pp]], 1, function(x){
          x_mat <- UT2mat(x) #form into a matrix (upper triangular)
          x_mat_inv <- backsolve(x_mat, diag(nL)) #take the inverse of an UT matrix fast
          #x_mat_inv_reo <- x_mat_inv[op,] #if we reverse the pivot -- no longer upper triangular!  so we will need to do this later.
          x_mat_inv[upper.tri(x_mat_inv, diag=TRUE)] #the inverse of an UT matrix is also upper triangular
        }, simplify=TRUE)) #this is extremely fast

      } #end loop over pivots

      #compute mean across FC samples
      M_all <- FC_nPivots*FC_nSamp2 #total number of samples
      FC_samp_mean <- Reduce("+", FC_samp_list)/M_all #mean of FC
      FC_samp_var <- Reduce("+", lapply(FC_samp_list, function(x){ (x - FC_samp_mean)^2 }))/M_all #var of FC

      # Compute the maximum eigenvalue of each FC^(-1) (same as 1 / min eigenvalue of each FC sample)
      FC_samp_maxeig <- sapply(FC_samp_list, function(x){
        vals <- eigen(x, only.values = TRUE)$values
        return(1/min(vals))
      })

      template$FC_Chol <- list(Chol_samp = Chol_samp, #pivoted Cholesky factors for every sample
                               FC_samp_logdet = FC_samp_logdet, #log determinant values for every sample
                               FC_samp_cholinv = FC_samp_cholinv, #pivoted Cholesky inverses for every sample
                               FC_samp_maxeig = FC_samp_maxeig, #maximum eigenvalue of inverse FC samples
                               FC_samp_mean = FC_samp_mean, #mean of FC samples
                               FC_samp_var = FC_samp_var, #var of FC samples
                               Chol_svd = Chol_svd,
                               pivots = pivots) #need to use these along with FC_samp_cholinv to determine inv(FC)
    } #end Cholesky-based FC template estimation
  }

  # [TO DO]: replace with fMRIscrub::unmask_mat or a vector version
  unmask_vec <- function(vec, mask) {
    vec2 <- rep(NA, length(mask))
    vec2[mask] <- vec
    vec2
  }

  # Format result ---------------------------------------------------
  if (use_mask2) {
    sigma_sq0 <- unmask_vec(sigma_sq0, mask2)
  }

  # Keep DR estimate of S
  if (!isFALSE(keep_S)) {
    DR0 <- array(DR0, dim=c(nM, nN, nL, nVm)) # Undo vectorize
    if (use_mask2) {
      DR0temp <- array(NA, dim=c(dim(DR0)[seq(3)], length(mask2)))
      DR0temp[,,,mask2] <- DR0
      DR0 <- DR0temp
    }
    if (is.character(keep_S)) {
      saveRDS(DR0, keep_S)
      keep_S <- FALSE # no longer need it.
    } else if (!isTRUE(keep_S)) {
      warning("`keep_S` should be `TRUE`, `FALSE`, or a file path. Using `FALSE`.")
      keep_S <- FALSE
    }
  }
  if (!keep_S) { rm(DR0) }

  # Keep DR estimate of FC
  if (FC) {
    if (!isFALSE(keep_FC)) {
      if (is.character(keep_FC)) {
        saveRDS(FC0, keep_FC) # in this case we don't save the Cholesky factors
        keep_FC <- FALSE # no longer need it.
      } else if (!isTRUE(keep_FC)) {
        warning("`keep_FC` should be `TRUE`, `FALSE`, or a file path. Using `FALSE`.")
        keep_FC <- FALSE
      }
    }
    if (!keep_FC) rm(FC0)
  }

  tparams <- list(
    FC=FC, FC_nPivots=FC_nPivots, FC_nSamp=FC_nSamp,
    num_subjects=nN, num_visits=nM,
    inds=inds,
    GSR=GSR, scale=scale,
    scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf,
    Q2=Q2, Q2_max=Q2_max,
    brainstructures=brainstructures, resamp_res=resamp_res,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    pseudo_retest=!real_retest
  )

  # If masking was performed, and if verbose, report the NA count in templates.
  # If no masking was performed, remove unnecessary metadata.
  if (use_mask2) {
    if (verbose) {
      cat("Template locations removed (filled with `NA`) due to varTol & missingTol:", sum(!mask2), "\n")
    }
  } else {
    mask2 <- NULL
  }

  if (GICA_parc && FORMAT=="CIFTI") {
    xii1 <- ciftiTools::convert_xifti(xii1, "dscalar")
  }

  dat_struct <- switch(FORMAT,
    CIFTI = ciftiTools::newdata_xifti(xii1, 0),
    GIFTI = list(hemisphere=ghemi),
    MATRIX = NULL
  )

  # Results list.
  result <- list(
    template=template,
    var_decomp=var_decomp,
    sigma_sq0=sigma_sq0,
    mask_input=mask,
    mask=mask2,
    params=tparams,
    dat_struct=dat_struct,
    GICA_parc_table=GICA_parc_table
  )

  # Add DR if applicable.
  if (keep_S) { result$S <- DR0 }
  if (FC && keep_FC) { result$FC <- FC0 }

  # Return results.
  class(result) <- paste0("template.", tolower(FORMAT))
  result
}

#' @rdname estimate_template
#'
estimate_template.cifti <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("local", "global", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01,
  GSR=FALSE,
  Q2=0, Q2_max=NULL,
  brainstructures="all", resamp_res=resamp_res,
  keep_S=FALSE, keep_FC=FALSE,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  usePar=FALSE, wb_path=NULL,
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    scale=scale, scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
    scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf,
    GSR=GSR,
    Q2=Q2, Q2_max=Q2_max,
    brainstructures=brainstructures, resamp_res=resamp_res,
    keep_S=keep_S,
    FC=FC,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    usePar=usePar, wb_path=wb_path,
    verbose=verbose
  )
}

#' @rdname estimate_template
#'
estimate_template.gifti <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("local", "global", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01,
  GSR=FALSE,
  Q2=0, Q2_max=NULL,
  brainstructures="all",
  keep_S=FALSE,keep_FC=FALSE,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  usePar=FALSE, wb_path=NULL,
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    scale=scale, scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
    scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf,
    GSR=GSR,
    Q2=Q2, Q2_max=Q2_max,
    brainstructures=brainstructures,
    keep_S=keep_S,
    FC=FC,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    usePar=usePar, wb_path=wb_path,
    verbose=verbose
  )
}

#' @rdname estimate_template
#'
estimate_template.nifti <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("local", "global", "none"),
  TR=NULL, hpf=.01,
  GSR=FALSE,
  Q2=0, Q2_max=NULL,
  mask=NULL,
  keep_S=FALSE,keep_FC=FALSE,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  usePar=FALSE, wb_path=NULL,
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    scale=scale,
    TR=TR, hpf=hpf,
    GSR=GSR,
    Q2=Q2, Q2_max=Q2_max,
    mask=mask,
    keep_S=keep_S,
    FC=FC,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    usePar=usePar, wb_path=wb_path,
    verbose=verbose
  )
}

