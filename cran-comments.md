## Test environments

* Windows x86_64-w64-mingw32/x64, R 4.2.2
* Mac x86_64-apple-darwin17.0, R 4.3.1

## R CMD check results

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded

## Downstream dependencies

None.

## Tests

Passes all the tests in `tests/testthat.R`

## Previous submission

  Package CITATION file contains call(s) to old-style personList() or
  as.personList().  Please use c() on person objects instead.

Fixed!

  checkRd: (-1) estimate_template.Rd:238: Lost braces; missing escapes or markup?
     238 | the type of scaling and detrending performed; the {dat_struct} which can be
         |                                                   ^
  checkRd: (-1) make_mask.Rd:15: Lost braces; missing escapes or markup?
      15 | {varTol}.}
         | ^

Fixed!