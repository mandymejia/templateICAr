# Build --> Install and Restart

# Setup ------------------------------------------------------------------------
# ciftiTools
library(ciftiTools)
print(packageVersion("ciftiTools"))
ciftiTools.setOption("wb_path", "~/Desktop/workbench")

# templateICAr
library(templateICAr)
# roxygen2::roxygenize("../../templateICAr")
print(packageVersion("templateICAr"))

library(RNifti)
library(gifti)
library(rgl)

# file paths
data_dir <- "tests/testWithData/data"
subjects <- c(100307, 100408, 100610)
cii_fnames <- c(
  paste0(data_dir, "/", subjects, "_rfMRI_REST1_LR_Atlas.dtseries.nii"),
  paste0(data_dir, "/", subjects, "_rfMRI_REST2_LR_Atlas.dtseries.nii")
)
giiL_fnames <- gsub("dtseries.nii", "sep.L.func.gii", cii_fnames, fixed=TRUE)
giiL_ROI_fnames <- gsub("dtseries.nii", "sep.ROI_L.func.gii", cii_fnames, fixed=TRUE)
nii_fnames <- gsub("dtseries.nii", "nii.gz", cii_fnames, fixed=TRUE)
rds_fnames <-gsub("dtseries.nii", "rds", cii_fnames, fixed=TRUE)

GICA_fname <- c(
  cii = file.path(data_dir, "melodic_IC_100.4k.dscalar.nii"),
  gii = file.path(data_dir, "melodic_IC_100.4k.sep.L.func.gii"),
  nii = file.path(data_dir, "melodic_IC_sum.nii.gz"),
  rds = file.path(data_dir, "melodic_IC_100.4k.rds")
)

xii1 <- select_xifti(read_cifti(GICA_fname["cii"]), 1) * 0

# `estimate_template`: check for same result w/ different file types -----------
### Test 1: basic
tm_cii <- estimate_template(
  cii_fnames[seq(5)], brainstructures="left", GICA = GICA_fname["cii"],
  keep_DR=TRUE, FC=TRUE
)
tm_gii <- estimate_template(
  giiL_fnames[seq(5)], GICA = GICA_fname["gii"],
  keep_DR=TRUE, FC=TRUE
)
tm_rds <- estimate_template(
  rds_fnames[seq(5)], GICA = GICA_fname["rds"],
  keep_DR=TRUE, FC=TRUE
)
testthat::expect_equal(
  lapply(tm_cii$template[seq(3)], fMRItools:::unmask_mat, tm_cii$dat_struct$meta$cortex$medial_wall_mask$left),
  tm_gii$template[seq(3)]
)
testthat::expect_equal(tm_cii$template, tm_rds$template)
##### Misc follow-up
tm_cii; tm_gii; tm_rds
plot(tm_cii); plot(tm_gii)

### Test 2: with various parameters changed
tm_cii <- estimate_template(
  cii_fnames[seq(3)], cii_fnames[seq(4,6)], GICA = GICA_fname["cii"],
  inds=c(2,7,11,90), scale="local", scale_sm_FWHM=5, detrend_DCT=4, normA=TRUE,
  maskTol=.9, brainstructures="left", wb_path="~/Desktop/workbench",
  usePar=TRUE, FC=TRUE, varTol=10000
)
tm_gii <- estimate_template(
  giiL_fnames[seq(3)], giiL_fnames[seq(4,6)], GICA = GICA_fname["gii"],
  inds=c(2,7,11,90), scale="local", scale_sm_FWHM=5, detrend_DCT=4, normA=TRUE,
  maskTol=.9, wb_path="~/Desktop/workbench",
  usePar=TRUE, FC=TRUE, varTol=10000
)
testthat::expect_equal(
  lapply(tm_cii$template[seq(3)], fMRItools:::unmask_mat, tm_cii$dat_struct$meta$cortex$medial_wall_mask$left),
  tm_gii$template[seq(3)]
)
tm_gii <- estimate_template(
  giiL_fnames[seq(3)], giiL_fnames[seq(4,6)], GICA = GICA_fname["gii"],
  inds=seq(5), scale="none", detrend_DCT=4, Q2=5,
  maskTol=.9, wb_path="~/Desktop/workbench",
  usePar=TRUE
)
tm_rds <- estimate_template(
  lapply(rds_fnames[seq(4,6)], readRDS), lapply(rds_fnames[seq(3)], readRDS), GICA = GICA_fname["rds"],
  inds=seq(5), scale="none", detrend_DCT=4, Q2=5,
  maskTol=.9, ,usePar=TRUE
)
testthat::expect_equal(
  tm_gii$template,
  lapply(tm_rds$template, fMRItools:::unmask_mat, tm_cii$dat_struct$meta$cortex$medial_wall_mask$left),
)

rgl.close(); rgl.close(); rgl.close(); rgl.close()

# CIFTI ------------------------------------------------------------------------
tm <- estimate_template(
  cii_fnames[seq(4)], GICA=cgIC_fname, scale=FALSE, keep_DR=TRUE#, FC=TRUE
)
tm
plot(tm)
rgl.close(); rgl.close()

cii <- read_cifti(cii_fnames[5])
cii$data$cortex_left[33,] <- mean(cii$data$cortex_left[33,])
tICA <- templateICA(cii, tm, scale=FALSE, maxiter=7, Q2=0)
plot(tICA)
rgl.close()
actICA <- activations(tICA)
actICA_fname <- paste0(tempfile(), ".dlabel.nii")
write_cifti(actICA, actICA_fname)
plot(actICA)
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
cii[[4]]$data$cortex_left[,] <- NA
tm <- estimate_template(
  cii, GICA=read_cifti(cgIC_fname, brainstructures="left"),
  scale="global", inds=c(1,4,7,11), maskTol = .5, missingTol=.5
)
tm
rm(cii)
plot(tm, idx=3)
rgl.close(); rgl.close()
cii <- read_cifti(cii_fnames[5], brainstructures="left")
cii$data$cortex_left[33,] <- mean(cii$data$cortex_left[33,])
tICA <- templateICA(cii, tm, brainstructures="left", scale="global", maxiter=7, Q2=0, spatial_model = TRUE)

cii <- lapply(cii_fnames[seq(4)], read_xifti, brainstructures="right")
cii0 <- lapply(cii, as.matrix)
cii0f <- paste0(c(tempfile(), tempfile(), tempfile(), tempfile()), ".rds")
tm <- estimate_template(cii0f, GICA=as.matrix(read_cifti(cgIC_fname, brainstructures="right")))


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
tICA <- templateICA(
  nii_fnames[2], tm, scale=FALSE,
  maxiter=1, mask=mask_fname, Q2=0
)
tICA
activations(tICA)

# GIFTI ------------------------------------------------------------------------
