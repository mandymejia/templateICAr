## Test environments

* Windows x86_64-w64-mingw32/x64, R 4.2.2
* Mac x86_64-apple-darwin17.0, R 4.1.1

## R CMD check results

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded

*This R extension was previously on CRAN before it was archived.*

## Downstream dependencies

None.

## Tests

Passes all the tests in `tests/testthat.R`

# Resubmission 1 

  Please add \value to .Rd files regarding exported methods and explain
  the functions results in the documentation. Please write about the
  structure of the output (class) and also what the output means. (If a
  function does not return a value, please document that too, e.g.
  \value{No return value, called for side effects} or similar)
  Missing Rd-tags in up to 12 .Rd files, e.g.:
        Gibbs_AS_posteriorCPP.Rd: \value
        resample_template.Rd: \value
        summary.template.cifti.Rd: \value
        summary.template.gifti.Rd: \value
        summary.template.matrix.Rd: \value
        summary.template.nifti.Rd: \value
        ...

Done! All exported functions should have \value now.

  \dontrun{} should only be used if the example really cannot be executed
  (e.g. because of missing additional software, missing API keys, ...) by
  the user. That's why wrapping examples in \dontrun{} adds the comment
  ("# Not run:") as a warning for the user. Does not seem necessary.
  Please replace \dontrun with \donttest.

  Please unwrap the examples if they are executable in < 5 sec, or replace
  dontrun{} with \donttest{}.

The instances of \dontrun indeed cannot be executed, because they require 
intermediate results that would be too large to host within the package.
The one for `estimate_template` shows the syntax for a typical use case, which
takes as input multiple data file paths. And all other instances of \dontrun are
for functions that require as input the template S3 object created by
`estimate_template`. Since these examples illustrate how to use the function for
the most typical use cases, but they cannot be run due to referencing objects
that don't exist (and cannot be created easily on-the-fly), we think \dontrun
is more appropriate than \donttest, and more appropriate than removing the
examples outright.
