## Test environments

* Mac x86_64-apple-darwin17.0, R 4.4.0

## R CMD check results

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded

## Downstream dependencies

None.

## Previous submission

Found the following Rd file(s) with Rd \link{} targets missing package
anchors:
  bdiag_m.Rd: dgCMatrix-class
  bdiag_m2.Rd: dgCMatrix-class

Found the following HTML validation problems:
bdiag_m.html:40:9 (bdiag_m.Rd:12): Error: <mat1> is not recognized!
bdiag_m.html:40:9 (bdiag_m.Rd:12): Warning: discarding unexpected <mat1>
bdiag_m.html:40:17 (bdiag_m.Rd:12): Error: <mat2> is not recognized!
bdiag_m.html:40:17 (bdiag_m.Rd:12): Warning: discarding unexpected <mat2>

* These `Rd` file problems have been fixed.

Examples with CPU (user + system) or elapsed time > 10s
                   user system elapsed
estimate_template 15.35   0.32   15.66

* This example has been revised to use smaller data so it finishes faster. 

## Tests

Passes all the tests in `tests/testthat.R`
