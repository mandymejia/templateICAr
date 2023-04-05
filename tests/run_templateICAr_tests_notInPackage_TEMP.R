# Build --> Install and Restart

# Setup ------------------------------------------------------------------------
# ciftiTools
library(ciftiTools)
print(packageVersion("ciftiTools"))
ciftiTools.setOption("wb_path", "~/Desktop/workbench")

# templateICAr
roxygen2::roxygenize("../../templateICAr")
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

# `estimate_template`: check for same result w/ different file types -----------
### Test 1: basic
tm_cii <- estimate_template(
  cii_fnames[seq(5)], brainstructures="left", GICA = GICA_fname["cii"],
  keep_DR=TRUE, FC=FALSE
)

# `templateICA`
tICA_cii <- templateICA(cii_fnames[4], brainstructures="left", tm_cii, maxiter=20, Q2=0)
