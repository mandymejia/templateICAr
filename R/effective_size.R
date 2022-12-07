# Look at the effective size of the first (S) and sixth (FC) result element
# and calculate the standard error
load("~/Downloads/testsubj1.RData")
names(result0)
names(result1) # This is the one with the MCMC samples

library(coda)
str(result1[[1]])
str(result1[[6]])

# Effective sample size for A
neff_A <- apply(result1[[1]],1:2, coda::effectiveSize)
apply(neff_A,2,summary)
#               [,1]      [,2]      [,3]      [,4]      [,5]
# Min.     367.9891  435.6564  408.9003  451.8461  453.7103
# 1st Qu.  850.0107  833.0011  826.8523  859.5236  828.1696
# Median   950.0000  950.0000  950.0000  950.0000  950.0000
# Mean     916.3597  898.8461  890.9201  923.1200  897.8728
# 3rd Qu.  950.0000  950.0000  950.0000  950.0000  950.0000
# Max.    1960.1537 1658.2725 1570.4217 1686.9252 2712.3933
# Sometimes you can get an effective sample size greater than the number of
# samples when you have negative autocorrelation in the sampling paths
# https://stats.stackexchange.com/questions/296059/effective-sample-size-greater-than-actual-sample-size

effic_A <- neff_A / 950 # There were 950 samples taken from the posterior
apply(effic_A,2,summary)
#               [,1]      [,2]      [,3]      [,4]      [,5]
# Min.    0.3873569 0.4585857 0.4304213 0.4756275 0.4775898
# 1st Qu. 0.8947482 0.8768432 0.8703709 0.9047617 0.8717575
# Median  1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
# Mean    0.9645892 0.9461538 0.9378106 0.9717053 0.9451293
# 3rd Qu. 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
# Max.    2.0633196 1.7455500 1.6530755 1.7757107 2.8551508

# Effective sample size for correlation
neff_corr <- apply(result1[[6]],1,coda::effectiveSize)
neff_corr_mat <- matrix(neff_corr, 5,5)
neff_corr_mat / 950 # This is the efficiency
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.0000000 0.2624071 0.2049072 0.2756693 0.2028445
# [2,] 0.2624071 0.0000000 0.2966938 0.2323800 0.2777222
# [3,] 0.2049072 0.2966938 0.0000000 0.2666900 0.2219476
# [4,] 0.2756693 0.2323800 0.2666900 0.0000000 0.2477732
# [5,] 0.2028445 0.2777222 0.2219476 0.2477732 0.0000000
# These are a bit low, but I'm not surprized here. Variance measures usually
# mix poorly in MCMC.

# Standard error of MCMC results using the effective sample size
# For A
sd_A <- apply(result1[[1]],1:2,sd)
se_A <- sd_A / sqrt(neff_A)
apply(se_A,2,summary)
#                 [,1]        [,2]        [,3]        [,4]        [,5]
# Min.    0.0009781418 0.001235335 0.001159358 0.001062199 0.000802049
# 1st Qu. 0.0014316851 0.001613056 0.001542167 0.001433226 0.001364720
# Median  0.0014734142 0.001675624 0.001604174 0.001488161 0.001412692
# Mean    0.0015076353 0.001709824 0.001630352 0.001531343 0.001436178
# 3rd Qu. 0.0015568984 0.001776887 0.001687830 0.001588480 0.001492277
# Max.    0.0024409376 0.002583843 0.002480831 0.002246293 0.002084771
# These look pretty small, which seems like a good thing. I would expect low
# numbers here after the scaling.

# For the correlation
sd_corr <- apply(result1[[6]],1,sd)
sd_corr_mat <- matrix(sd_corr,5,5)
(se_corr <- sd_corr_mat[upper.tri(sd_corr_mat)] / neff_corr_mat[upper.tri(neff_corr_mat)])
# [1] 3.833044e-05 6.169135e-05 4.451739e-05 3.586047e-05 5.526549e-05 4.777773e-05 7.201388e-05 5.278028e-05 6.910744e-05 5.452063e-05
# Again, very small values, as expected.
