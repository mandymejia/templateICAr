# Build --> Install and Restart

# Setup ------------------------------------------------------------------------
# ciftiTools
library(ciftiTools)
print(packageVersion("ciftiTools"))
ciftiTools.setOption("wb_path", "~/../Desktop/fMRI/workbench")

# templateICAr
library(templateICAr)
# roxygen2::roxygenize("../../templateICAr")
print(packageVersion("templateICAr"))

library(RNifti)
library(rgl)

# file paths
data_dir <- "tests/testWithData/data"
subjects <- c(100307, 100408, 100610)
cii_fnames <- c(
  paste0(data_dir, "/", subjects, "_rfMRI_REST1_LR_Atlas.dtseries.nii"),
  paste0(data_dir, "/", subjects, "_rfMRI_REST2_LR_Atlas.dtseries.nii")
)
nii_fnames <- c(
  paste0(data_dir, "/", subjects, "_rfMRI_REST1_LR.nii.gz"),
  paste0(data_dir, "/", subjects, "_rfMRI_REST2_LR.nii.gz")
)

# CIFTI ------------------------------------------------------------------------

# Load CIFTI group IC
cgIC_fname <- file.path(data_dir, "melodic_IC_100.4k.dscalar.nii")
cgIC <- read_cifti(cgIC_fname)
xii1 <- select_xifti(cgIC, 1) * 0

tm <- estimate_template(
  cii_fnames[seq(4)], GICA=cgIC_fname, scale=FALSE, keep_DR=TRUE
)
tm
plot(tm)
rgl.close(); rgl.close()
tICA <- templateICA(cii_fnames[5], tm, scale=FALSE, maxiter=7)
plot(tICA)
rgl.close()
plot(activations(tICA))
rgl.close()

tm <- estimate_template(
  cii_fnames[seq(3)], cii_fnames[seq(4, 6)],
  GICA=cgIC_fname, scale="local", detrend_DCT=3,
  normA=TRUE, brainstructures="right", varTol=1, verbose=FALSE
)
tm
plot(tm, "var")
rgl.close()

tm2 <- estimate_template(
  cii_fnames[seq(3)], cii_fnames[seq(4, 6)],
  GICA=cgIC_fname, scale="local", detrend_DCT=3, scale_FWHM=20,
  normA=TRUE, brainstructures="right", varTol=1, verbose=FALSE
)

cii <- lapply(cii_fnames[seq(4)], read_xifti, brainstructures="left")
cii[[1]]$data$cortex_left[3,17] <- NA
cii[[1]]$data$cortex_left[11,5] <- NA
cii[[2]]$data$cortex_left[11,seq(10)] <- NA
cii[[3]]$data$cortex_left[11,] <- NA
cii[[1]]$data$cortex_left[78,5] <- NA
cii[[2]]$data$cortex_left[78,seq(10)] <- NA
cii[[3]]$data$cortex_left[78,] <- NA
tm <- estimate_template(
  cii, GICA=read_cifti(cgIC_fname, brainstructures="left"),
  scale="global", inds=c(1,4,7,11)
)
tm
rm(cii)
plot(tm, idx=3)
rgl.close(); rgl.close()

# CIFTI pseudo retest vs data true retest: should get same results.
tm2 <- estimate_template(
  lapply(cii, function(x){as.matrix(x)[,seq(600)]}),
  lapply(cii, function(x){as.matrix(x)[,seq(601,1200)]}),
  GICA=as.matrix(read_cifti(cgIC_fname, brainstructures="left")),
  scale="global", inds=c(1,4,7,11),
)
stopifnot(
  max(abs(do.call(c, tm$var_decomp) - do.call(c, tm2$var_decomp)), na.rm=TRUE) < 1e-8
)
rm(tm2)

tICA <- templateICA(cii_fnames[2], tm, brainstructures="left")
tICA
plot(tICA)
rgl.close()
# plot(activations(tICA))

tICA <- templateICA(
  cii_fnames[3], tm, brainstructures="left",
  tvar_method="unbiased", Q2=0, reduce_dim=FALSE, usePar=TRUE
)
tICA
plot(tICA)
rgl.close()
# temp:
tICA$subjICmean$data$cortex_left[is.na(tICA$subjICmean$data$cortex_left)] <- 0
tICA$subjICse$data$cortex_left[is.na(tICA$subjICse$data$cortex_left)] <- 5
tICA$mask <- rep(TRUE, length(tICA$mask))
# -----
plot(activations(tICA))
rgl.close()

# NIFTI ------------------------------------------------------------------------

rm(cgIC, xii1)

# Load NIFTI group IC
ngIC_fname <- file.path(data_dir, "melodic_IC_sum.nii.gz")
ngIC <- readNifti(ngIC_fname)
nmask <- apply(ngIC!=0, seq(3), all)
mask_fname <- file.path(data_dir, "mask.nii.gz")
RNifti::writeNifti(nmask, mask_fname)

# mask erosion?
tm <- estimate_template(
  nii_fnames[seq(4)], GICA=ngIC_fname, scale=FALSE,
  keep_DR=TRUE, mask=mask_fname, varTol = 500, maskTol=.3, missingTol=.9
)
tm
tICA <- templateICA(nii_fnames[2], tm, scale=FALSE, maxiter=7, mask=mask_fname)
tICA
activations(tICA)
