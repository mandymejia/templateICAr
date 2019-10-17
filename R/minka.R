minka = function(X,alpha=0,beta=0,eps = 1e-12){

  #defaulting alpha and beta to 0 can be dangerous. please don't do it, please
  #probably better to default to 1 #TO DO: ask Andy about this
  if(alpha<=0){alpha=0;beta=0}

  #X is a matrix with observations in rows
  N = dim(X)[1]
  d = dim(X)[2]

  if(N > d){
    covX = cov(X)*(N-1)/N #S/N
    svd_covX = svd(covX)
    lambda = svd_covX$d
    #U=svd_covX$u
  } else {
    XXt = X %*% t(X) / N
    svd_XXt = svd(XXt)
    lambda = svd_XXt$d
  }
  lambda = lambda[lambda>eps]
  rankX =length(lambda)


  #replace d with rank, this should be done for model to make sense if X is rank deficient
  replace_d_with_rankX = FALSE
  if(replace_d_with_rankX==TRUE){D=rankX}else{D=d}
  #we will use D where using rankX could make more sense than using dimX
  maxK = rankX  - 1

  #k values
  k = c(1:maxK)
  #do everything in the log of the densities
  #try to do all values of k at once when possible
  #n = N+alpha-1
  n = N+ifelse(alpha==0,0,alpha-1)

  #mult correction
  multiplicity_correction = k*log(2)
  #print((multiplicity_correction))
  # power of 2pi on density for x
  pow_2pi_likelihood = -(N*D)/2*log(2*pi)
  #print((pow_2pi_likelihood))
  #N^(-d/2) from integrating out the m
  m_int_const = -D/2*log(N)+D/2*log(2*pi)
  #print((m_int_const))
  #prior const for U
  i = c(1:maxK)
  prior_const_U = -k*log(2)+cumsum(lgamma((d-i+1)/2)-(d-i+1/2)*log(pi))
  #print((prior_const_U))
    #cumsum is doing the sum (product) over the i=1:k
    #the d is the ambient dimension if X is rank deficient, not the rank
  #prior const for v and L
  prior_const_L = if(alpha==0){0}else{(alpha/2*log(beta/2)-lgamma(alpha/2))*k}
  #print((prior_const_L))
  prior_const_v = if(alpha==0){0}else{alpha*(D-k)/2*log(beta*(D-k)/2)-lgamma(alpha*(D-k)/2)}
  #print((prior_const_v))
  #det L part
  det_L_part = -n/2*cumsum(log(lambda[1:maxK]))
  #print((det_L_part))
  #v_hat_part
  v_hat = 1/(n*(D-k))*(cumsum(N*lambda[rankX:2]+beta))[maxK:1]
  v_hat_part = -n*(D-k)/2*log(v_hat)
  #print((v_hat_part))
  #exponential part
  exp_part = -n*D/2
  #print((exp_part))
  #2pi power k = dim lambda, 1 = dim v, d*k-k-choose(k,2)  = dim U
  two_pi_part_posterior_approx = (k+1+d*k-k-choose(k,2))/2*log(2*pi)
  #print((two_pi_part_posterior_approx))
  #AL part
  AL_part = -k/2*log(n/2)
  #print((AL_part))
  #Av part
  Av_part = -1/2*log(n*(D-k)/2)
  #print((Av_part))
  #Az part
  #part 1 is the difference in lambda part
  Az_part_1 = rep(0,maxK)
  for(i in 1:maxK){
    j_index = (i+1):rankX
    Az_part_1[i] = sum(log(lambda[i]-lambda[j_index]))
  }
  Az_part_1 = -1/2*cumsum(Az_part_1)
  #print((Az_part_1))
  #part 2 is for the inverse of the lambda hats
  #see eq 70 in updated paper
  Az_part_2 = rep(0,maxK)
  lambda_hat_leqk = (lambda[1:maxK]*N+beta)/(n) #note, name of variable applies when the index is <= k for a given k
  log_diff_inv_lambda_hat_leqk = rep(0,maxK)
  for(j in 2:maxK){
    i_index = 1:(j-1)
    log_diff_inv_lambda_hat_leqk[j] = sum(log(lambda_hat_leqk[j]^(-1)-lambda_hat_leqk[i_index]^(-1)))
  }
  Az_part_2 = -1/2*cumsum(log_diff_inv_lambda_hat_leqk)
  #print((Az_part_2))
  for(ell in 1:maxK){
    Az_part_2[ell] = Az_part_2[ell] - (d-ell)/2*sum(log(v_hat[ell]^(-1)-lambda_hat_leqk[1:ell]^(-1)))
  }
  Az_part_3 = -1/2*(d*k-(k*(k+1)/2))*log(N)
  #print((Az_part_3))

  return(list(Az_part_1, Az_part_2, Az_part_3) )
}


#   pow_2pi_likelihood_p0 = -(N*D)/2*log(2*pi)
#   m_int_const_p0 = -D/2*log(N)+D/2*log(2*pi)
#   prior_const_v_p0 = if(alpha==0){1}else{alpha*D/2*log(beta*D/2)-lgamma(alpha*D/2)}
#   v_hat_p0 = 1/(n*D)*(sum(N*lambda+beta))
#   v_hat_part_p0 = -n*D/2*log(v_hat_p0)
#   exp_part_p0 = -n*D/2
#   two_pi_part_posterior_approx_p0 = 1/2*log(2*pi)
#   Av_part_p0 = -1/2*log(n*D/2)
#   p0 = pow_2pi_likelihood_p0+
#     m_int_const_p0+prior_const_v_p0+
#     v_hat_part_p0+
#     exp_part_p0+
#     two_pi_part_posterior_approx_p0+
#     Av_part_p0
#
#   p = multiplicity_correction+
#     pow_2pi_likelihood+
#     m_int_const+
#     prior_const_U+
#     prior_const_L+
#     prior_const_v+
#     det_L_part+
#     v_hat_part+
#     exp_part+
#     two_pi_part_posterior_approx+
#     AL_part+
#     Av_part+
#     Az_part_1+
#     Az_part_2+
#     Az_part_3
#   p = c(p0,p)
#   names(p) = paste("rank",c(0:maxK),sep="=")
#   return(
#     p
#   )
# }
#
# # N = 100
# # D = 35
# # mu = rt(D,3)*4
# # k = 10
# # sd=2
# # H = matrix(rnorm(D*k),D,k)
# # w = matrix(rnorm(N*k),k,N)
# # lrm = t(H%*%w)
# # X = matrix(rnorm(N*D,sd=sd),N,D)+matrix(mu,N,D,byrow=TRUE)+lrm
# # out = minka(X,1,1)
# probs = out-max(out)
# probs = exp(probs)
# probs = probs/sum(probs)
# probs[probs==max(probs)]
# head(sort(probs,decreasing = TRUE))


