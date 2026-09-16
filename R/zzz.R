.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Note: templateICAr has been superseded by BayesBrainMap,\n",
    "which includes model improvements and renamed functions.\n",
    "New projects should use BayesBrainMap! Install with:\n",
    "> install.packages('BayesBrainMap\')"
  )
}
