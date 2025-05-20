# Build --> Install and Restart

# [Edit these]
## path to your Connectome Workbench
my_wb <- "~/Applications"
## path to test data
dir_data <- "tests/data"
## path to results from tests
dir_results <- file.path(dir_data, "results")

library(testthat)
library(ciftiTools)
if (interactive()) { ciftiTools.setOption('wb_path', my_wb) }
library(templateICAr)

tests_dir <- "testthat"
if (!endsWith(getwd(), "tests")) { tests_dir <- file.path("tests", tests_dir) }
source(file.path(tests_dir, "test-misc.R"))
