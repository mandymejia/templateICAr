#' Perform group ICA based on CIFTI data
#'
#' @param cifti_fnames Vector of file paths of CIFTI-format fMRI timeseries
#'  (*.dtseries.nii) for template estimation
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param num_PCs Maximum number of PCs to retain in initial subject-level PCA
#' @param num_ICs Number of independent components to identify.
#' @param max_rows_GPCA The maximum number of rows for the matrix on which group
#' level PCA will be performed.  This matrix is the result of temporal concatenation
#' and contains (number of subjects)*(number of subject-level PCs) rows.
#' @param verbose If \code{TRUE}, display progress updates
#' @param out_fname (Optional) If not specified, a xifti object will be returned but
#' the GICA maps will not be written to a CIFTI file.
#' @param surfL (Optional) File path to a surface GIFTI for the left cortex.
#' If provided, will be part of xifti result object for visualization in R.
#' @param surfR (Optional) File path to a surface GIFTI for the right cortex.
#' If provided, will be part of xifti result object for visualization in R.
#'
#' @importFrom ciftiTools read_cifti write_cifti make_surf
#' @importFrom ica icaimax
#'
#' @return xifti object containing group ICA maps
#'
#' @export
#'
groupICA.cifti <- function(
  cifti_fnames,
  brainstructures=c("left","right"),
  num_PCs=100,
  num_ICs,
  max_rows_GPCA=10000,
  verbose=TRUE,
  out_fname=NULL,
  surfL=NULL,
  surfR=NULL){

  if(!is.null(out_fname)){
    if(!dir.exists(dirname(out_fname))) stop('directory part of out_fname does not exist')
  }

  notthere <- !file.exists(cifti_fnames)
  if(sum(notthere) == length(cifti_fnames)) stop('The files in cifti_fnames do not exist.')
  if(sum(notthere) > 0) warning(paste0('There are ', sum(notthere), ' files in cifti_fnames that do not exist. These will be excluded from the group ICA.'))
  cifti_fnames <- cifti_fnames[!notthere]

  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }

  N <- length(cifti_fnames)
  if(verbose){ cat(paste0('\n Number of subjects: ',N)) }

  if(N*num_PCs > max_rows_GPCA) warning(paste0('The number of subjects times the number of
                                        subject-level PCs exceeds the limit of ',max_rows_GPCA,
                                        '. Some subjects will be excluded from analysis.'))
  size_bigY <- min(max_rows_GPCA, (num_PCs*N)) #for group-level PCA

  # READ IN DATA AND PERFORM INITIAL DIMENSION REDUCTION
  for(ii in 1:N){

    #read in BOLD
    if(verbose) cat(paste0('\n Reading data for subject ',ii,' of ',N))
    fname_ii <- cifti_fnames[ii]
    if(!file.exists(fname_ii)) {
      missing_data <- c(missing_data, fname_ii)
      if(verbose) cat(paste0('\n Data not available for file:', fname_ii))
      next
    }

    BOLD_ii <- do.call(rbind, read_cifti(fname_ii, brainstructures=brainstructures)$data)
    nvox <- nrow(BOLD_ii)
    ntime <- ncol(BOLD_ii)

    if(ii==1){
      #for temporal concatenation
      #bigY <- matrix(0, nrow=(num_PCs*N), ncol=nvox)
      bigY <- matrix(0, nrow=size_bigY, ncol=nvox) #for initializing online PCA
    } else {
      if(nvox != ncol(bigY)) stop(paste0('Number of brain locations in subject ',ii,' does not match first subject.'))
    }
    rows_ii <- (1:num_PCs) + (ii-1)*num_PCs

    #stop adding subjects once we reach the nrow limit
    if(max(rows_ii) > size_bigY){
      break
    } else {
      #perform subject-level PCA and add rows to bigY
      if(verbose) cat(paste0('... performing dimension reduction'))
      if(ntime < num_PCs) warning(paste0('In the file ',cifti_fnames[ii],' there are fewer than num_PCs (',num_PCs,') time points. Skipping dimension reduction.')) #note: still use SVD to standardize variance
      num_PCs_ii <- min(num_PCs, ntime-1)
      Y_ctr <- scale(t(BOLD_ii), scale=FALSE)
      YYt <- tcrossprod(Y_ctr) #1 min
      Y_svd <- svd(YYt, nv = 0, nu=num_PCs_ii) #5 sec
      U_ii <- Y_svd$u
      D_ii <- sqrt(Y_svd$d[1:num_PCs_ii])
      V_ii <- diag(1/D_ii) %*% t(U_ii) %*% Y_ctr

      rows_ii <- rows_ii[1:num_PCs_ii] #this will result in rows of all zeros if num_PCs_ii < num_PCs. Exclude at the end.
      bigY[rows_ii,] <- V_ii
    }
  }

  # PERFORM GROUP-LEVEL PCA FOR DIMENSION REDUCTION
  zero_rows <- (rowSums(abs(bigY)) == 0)
  bigY <- bigY[!zero_rows,]
  bigY_ctr <- scale(bigY, scale=FALSE)
  bigYYt <- tcrossprod(bigY_ctr) #40 min
  bigY_svd <- svd(bigYYt, nv = 0, nu=num_ICs) #45 min
  U_bigY <- bigY_svd$u
  D_bigY <- sqrt(bigY_svd$d[1:num_ICs])
  V_bigY <- diag(1/D_bigY) %*% t(U_bigY) %*% bigY_ctr

  # PERFORM GROUP-LEVEL ICA
  GICA <- icaimax(t(V_bigY), nc=num_ICs, center=TRUE)

  #fix the mean of each component to be positive
  sign_flip <- which(colMeans(GICA$S) < 0)
  GICA$S[,sign_flip] <- (-1)*GICA$S[,sign_flip]

  # Format GICA as a xifti object and write out a cifti

  xifti_GICA <- read_cifti(fname_ii, brainstructures=brainstructures) #as a template
  nleft <- nrow(xifti_GICA$data$cortex_left)
  nright <- nrow(xifti_GICA$data$cortex_right)
  nsub <- nrow(xifti_GICA$data$subcort)
  if ("left" %in% brainstructures) xifti_GICA$data$cortex_left <- GICA$S[1:nleft,]  #[flat_bs_mask == "left",, drop=FALSE]
  if ("right" %in% brainstructures) xifti_GICA$data$cortex_right <- GICA$S[nleft+(1:nright),] #[flat_bs_mask == "right",, drop=FALSE]
  if ("subcortical" %in% brainstructures) xifti_GICA$data$subcort <- GICA$S[nleft+nright+(1:nsub),] #[flat_bs_mask == "subcortical",, drop=FALSE]

  #format as a dscalar
  xifti_GICA$meta$cifti$intent <- 3006
  xifti_GICA$meta$cifti[c('time_start','time_step','time_unit')] <- NULL
  xifti_GICA$meta$cifti$names <- paste0('IC ',1:num_ICs)

  #add surfaces, if provided
  if(!is.null(surfL)) {
    surfL <- make_surf(surf = '~/Google Drive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii')
    xifti_GICA <- add_surf(xifti = xifti_GICA, surfL = surfL)
  }
  if(!is.null(surfR)) {
    surfR <- make_surf(surf = '~/Google Drive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii')
    xifti_GICA <- add_surf(xifti = xifti_GICA, surfR = surfR)
  }

  # FOR VISUALIZATION
  # for(k in 1:25){
  #   print(k)
  #   view_xifti_surface(xifti_GICA, idx=k, fname = paste0('~/Desktop/GICA_',k))
  # }


  if(!is.null(out_fname)){
    write_cifti(xifti_GICA, paste0(out_fname, '.dscalar.nii'), verbose=verbose)
  }

  return(xifti_GICA)
}

