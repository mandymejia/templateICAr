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
data_dir <- "data_notInPackage"
subjects <- c(100307, 100408, 100610)
cii_fnames <- c(
  paste0(data_dir, "/", subjects, "_rfMRI_REST1_LR_Atlas.dtseries.nii"),
  paste0(data_dir, "/", subjects, "_rfMRI_REST2_LR_Atlas.dtseries.nii")
)
giiL_fnames <- gsub("dtseries.nii", "sep.L.func.gii", cii_fnames, fixed=TRUE)
giiL_ROI_fnames <- gsub("dtseries.nii", "sep.ROI_L.func.gii", cii_fnames, fixed=TRUE)
nii_fnames <- gsub("_Atlas.dtseries.nii", ".nii.gz", cii_fnames, fixed=TRUE)
rds_fnames <-gsub("dtseries.nii", "rds", cii_fnames, fixed=TRUE)

GICA_fname <- c(
  cii = file.path(data_dir, "melodic_IC_100.4k.dscalar.nii"),
  gii = file.path(data_dir, "melodic_IC_100.4k.sep.L.func.gii"),
  nii = file.path(data_dir, "melodic_IC_sum.nii.gz"),
  rds = file.path(data_dir, "melodic_IC_100.4k.rds")
)

xii1 <- select_xifti(read_cifti(GICA_fname["cii"]), 1) * 0

# Quick little check of the three main functions, w/ CIFTI ---------------------
tm_cii <- estimate_template(
  cii_fnames[seq(3)], GICA = GICA_fname["cii"], TR=.72, FC=FALSE,
  brainstructures=c("left", "right")
)
tICA_cii <- templateICA(
  cii_fnames[4], tm_cii, brainstructures="left", maxiter=5, TR="template", resamp_res=2000
)
actICA_cii <- activations(tICA_cii)
actICA_cii <- activations(tICA_cii, z=c(0, .1, 3, 11))

# `estimate_template`: check for same result w/ different file types -----------
### Test 1: basic
tm_cii <- estimate_template(
  cii_fnames[seq(5)], brainstructures="left", GICA = GICA_fname["cii"],
  keep_DR=TRUE, FC=TRUE, TR=.72, scale_sm_FWHM=0
)
tm_gii <- estimate_template(
  giiL_fnames[seq(5)], GICA = GICA_fname["gii"],
  keep_DR=TRUE, FC=TRUE, TR=.72, scale_sm_FWHM=0
)
tm_rds <- estimate_template(
  rds_fnames[seq(5)], GICA = GICA_fname["rds"],
  keep_DR=TRUE, FC=TRUE, TR=.72, scale_sm_FWHM=0
)
testthat::expect_equal(
  lapply(tm_cii$template[seq(3)], fMRItools::unmask_mat, tm_cii$dat_struct$meta$cortex$medial_wall_mask$left),
  tm_gii$template[seq(3)]
)
testthat::expect_equal(tm_cii$template, tm_rds$template)
##### Misc follow-up
tm_cii; tm_gii; tm_rds
plot(tm_cii); plot(tm_gii)

### Test 2: with various parameters changed
tm_cii <- estimate_template(
  cii_fnames[seq(3)], cii_fnames[seq(4,6)], GICA = GICA_fname["cii"],
  inds=c(2,7,11,90), scale="global", scale_sm_FWHM=5,
  maskTol=.9, brainstructures="left", wb_path="~/Desktop/workbench",
  usePar=TRUE, FC=TRUE, varTol=10000
)
tm_gii <- estimate_template(
  giiL_fnames[seq(3)], giiL_fnames[seq(4,6)], GICA = GICA_fname["gii"],
  inds=c(2,7,11,90), scale="global", scale_sm_FWHM=5,
  maskTol=.9, wb_path="~/Desktop/workbench",
  usePar=TRUE, FC=TRUE, varTol=10000
)
testthat::expect_equal(
  lapply(tm_cii$template[seq(3)], fMRItools::unmask_mat, tm_cii$dat_struct$meta$cortex$medial_wall_mask$left),
  tm_gii$template[seq(3)]
)
tm_gii <- estimate_template(
  giiL_fnames[seq(3)], giiL_fnames[seq(4,6)], GICA = GICA_fname["gii"],
  inds=seq(5), scale="none", Q2=5,
  maskTol=.9, wb_path="~/Desktop/workbench",
  usePar=TRUE
)
tm_rds <- estimate_template(
  lapply(rds_fnames[seq(4,6)], readRDS), lapply(rds_fnames[seq(3)], readRDS), GICA = GICA_fname["rds"],
  inds=seq(5), scale="none",  Q2=5,
  maskTol=.9, usePar=TRUE
)
testthat::expect_equal(
  tm_gii$template,
  lapply(tm_rds$template, fMRItools::unmask_mat, tm_cii$dat_struct$meta$cortex$medial_wall_mask$left),
)

close3d(); close3d(); close3d(); close3d()

# `export_template` and `templateICA`: check for same result w/ different file types -----------------
tm_cii <- estimate_template(
  cii_fnames[seq(3)], brainstructures="left", GICA = GICA_fname["cii"], inds=seq(3),
  keep_DR=TRUE, scale="global"
)
tm_gii <- estimate_template(
  giiL_fnames[seq(3)], GICA = GICA_fname["gii"], inds=seq(3),
  keep_DR=TRUE, scale="global"
)
tm_rds <- estimate_template(
  rds_fnames[seq(3)], GICA = GICA_fname["rds"], inds=seq(3),
  keep_DR=TRUE, scale="global"
)

# `export_template`
out_fname=export_template(tm_cii, tempfile())
tm_cii2 <- list(read_cifti(out_fname[1]), read_cifti(out_fname[2]), readRDS(out_fname[3]))
out_fname=export_template(tm_gii, tempfile())
tm_gii2 <- list(readgii(out_fname[1]), readgii(out_fname[2]), readRDS(out_fname[3]))
out_fname=export_template(tm_rds, tempfile())
tm_rds2 <- lapply(out_fname, readRDS)

tm_cii <- estimate_template(
  cii_fnames[seq(3)], brainstructures="left", GICA = GICA_fname["cii"], inds=seq(3),
  keep_DR=TRUE, FC=FALSE
)

# `templateICA`
tICA_cii <- templateICA(cii_fnames[4], brainstructures="left", tm_cii, maxiter=20, Q2=0, TR=.72)
tICA_gii <- templateICA(giiL_fnames[4], tm_gii, Q2=0, maxiter=20, TR=.72)
tICA_rds <- templateICA(rds_fnames[4], tm_rds, Q2=0, maxiter=20, TR=.72)
tICA_cii; tICA_gii; tICA_rds
testthat::expect_equal(tICA_cii$theta_MLE, tICA_rds$theta_MLE)
testthat::expect_equal(tICA_gii$A, tICA_rds$A)
actICA_rds <- activations(tICA_rds)
actICA_cii <- activations(tICA_cii)
plot(activations(tICA_cii)); plot(activations(tICA_gii))
close3d(); close3d()

gamma <- 2
gamma_scaled <- gamma*sqrt(matrixStats::colVars(tICA_cii$template_mean))
act2 <- activations(tICA_cii, u=gamma_scaled)
act3 <- activations(tICA_cii, z=gamma)
testthat::expect_equal(act2[names(act2)!="z"], act3[names(act3)!="z"])

# CIFTI ------------------------------------------------------------------------
tm <- estimate_template(
  cii_fnames[seq(4)], GICA=GICA_fname["cii"], scale=FALSE, keep_DR=TRUE#, FC=TRUE
)
tm
plot(tm)
close3d(); close3d()

cii <- read_cifti(cii_fnames[5])
cii$data$cortex_left[33,] <- mean(cii$data$cortex_left[33,])
tICA <- templateICA(cii, tm, scale=FALSE, maxiter=7, Q2=0)
plot(tICA)
close3d()
actICA <- activations(tICA)
actICA_fname <- paste0(tempfile(), ".dlabel.nii")
# [TO DO]: ciftiTools 12.0. get rid of below.
actICA$active$data$cortex_left[50,] <- 0
write_cifti(actICA$active, actICA_fname)
actICA2 <- read_cifti(actICA_fname)
plot(actICA); plot(actICA2)
close3d(); close3d()

# LEFT OFF HERE.
tm <- estimate_template(
  cii_fnames[seq(3)], cii_fnames[seq(4, 6)],
  GICA=GICA_fname["cii"], scale="local",
  brainstructures="right", varTol=1, verbose=FALSE
)
tm
plot(tm, "var")
close3d()

tm2 <- estimate_template(
  cii_fnames[seq(3)], cii_fnames[seq(4, 6)],
  GICA=GICA_fname["cii"], scale="local", scale_sm_FWHM=20,
  brainstructures="right", varTol=1, verbose=FALSE
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
  cii, GICA=read_cifti(GICA_fname["cii"], brainstructures="left"),
  scale="global", inds=c(1,4,7,11), maskTol = .5, missingTol=.5
)
tm
rm(cii)
plot(tm, idx=3)
close3d(); close3d()
cii <- read_cifti(cii_fnames[5], brainstructures="left")
#squarem1 error ... ?
#templateICA(cii, tm, brainstructures="left", scale="global", maxiter=7, Q2=0, spatial_model = TRUE)
cii$data$cortex_left[33,] <- mean(cii$data$cortex_left[33,])
tICA <- testthat::expect_error( # Not supported yet: flat or NA voxels in data, after applying template mask, with spatial model.
  templateICA(cii, tm, brainstructures="left", scale="global", maxiter=7, Q2=0, spatial_model = TRUE, TR=.72)
)

cii <- lapply(cii_fnames[seq(4)], read_xifti, brainstructures="right")
cii0 <- lapply(cii, as.matrix)
cii0f <- paste0(c(tempfile(), tempfile(), tempfile(), tempfile()), ".rds")
for (ii in seq(4)) { saveRDS(cii0[[ii]], cii0f[ii]) }
tm <- estimate_template(
  cii0f,
  GICA=as.matrix(read_cifti(GICA_fname["cii"], brainstructures="right")),
  scale="global", inds=c(1,4,7,11)
)


# CIFTI pseudo retest vs data true retest: should get same results.
tm2 <- estimate_template(
  lapply(cii, function(x){as.matrix(x)[,seq(600)]}),
  lapply(cii, function(x){as.matrix(x)[,seq(601,1200)]}),
  GICA=as.matrix(read_cifti(GICA_fname["cii"], brainstructures="right")),
  scale="global", inds=c(1,4,7,11),
)
stopifnot(
  max(abs(do.call(c, tm$var_decomp) - do.call(c, tm2$var_decomp)), na.rm=TRUE) < 1e-8
)
rm(tm2)

# NIFTI ------------------------------------------------------------------------
rm(xii1)

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
  miniter=1, maxiter=1, mask=mask_fname, Q2=0, TR=.72
)
tICA
activations(tICA)
close3d()
