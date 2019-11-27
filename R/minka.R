#' Automated PCA-based Dimension Reduction based on Minka (2008)
#'
#' @param data Nxd data matrix for N observations of d variables.
#' @param alpha Parameter controlling sharpness of prior on eigenvalues and noise variance. Values close to zero lead to uninformative priors.
#' @param beta Parameter controlling sharpness of prior on eigenvalues and noise variance. Values close to zero lead to uninformative priors.
#' @param eps Value below which sample eigenvalues will be considered equal to zero.
#'
#' @return A list containing the estimated dimensionality, noise variance, and the vector of log probabilities for each dimensionality.
#' @export
#'
minka = function(data, alpha = 0, beta = 0, eps = 1e-10){

  if(alpha == 0){alpha = 0; beta = 0}

  N = dim(data)[1] #144 vs 145
  d = dim(data)[2]

  if(N < d){
    data0 = scale(data, scale=FALSE)
    lambda = svd(data0 %*% t(data0) / N)$d
  } else {
    data0 = scale(data, scale=FALSE)
    lambda = svd(t(data0) %*% data0 / N)$d
  }

  lambda = lambda[lambda > eps]
  rank_of_data = length(lambda)
  # if(rank_of_data < d) lambda[(rank_of_data+1):d] = 0
  # true_rank = rank_of_data
  # rank_of_data = d

  #logediff(i) = sum_{j>i} log(lambda(i) - lambda(j))
  logediff = rep(0, rank_of_data) #only consider values of k up to rank of data
  for (i in 1:(rank_of_data-1)){
    j = (i+1):rank_of_data
    logediff[i] = sum(log(lambda[i] - lambda[j])) + (d - rank_of_data)*log(lambda[i]) #accounting for missing eigenvalues
  }
  cumsum_logediff = cumsum(logediff)

  n = N - 1 + alpha
  lambda_hat = (lambda + beta/N)*N/n
  inv_lambda = 1/lambda_hat

  # loginvediff = rep(0, rank_of_data) #only consider values of k up to rank of data
  # for (i in 1:(rank_of_data-1)){
  #   j = (i+1):rank_of_data
  #   logediff[i] = sum(log(inv_lambda[j] - inv_lambda[i])) + (d - rank_of_data)*log(lambda[i]) #accounting for missing eigenvalues
  # }
  # cumsum_logediff = cumsum(logediff)



  r = matrix(inv_lambda, nrow=rank_of_data, ncol=rank_of_data, byrow=TRUE)
  c = matrix(inv_lambda, nrow=rank_of_data, ncol=rank_of_data)
  invediff = c - r
  invediff[invediff <= 0] = 1
  invediff = Matrix(log(invediff)) #store as sparse matrix because more than half zeros

  #rowSums_invediff_big <- rep(N*rowSums(invediff),N)
  #row_invediff <- cumsum(rowSums_invediff_big)

  # r = repmat(inv_lambda, 1, rank_of_data)
  # c = t(repmat(inv_lambda, 1, rank_of_data))
  # num = length(r) * length(c)
  # invediff = matrix(c(1:num), nrow = length(r))
  # for(i in 1:length(r)){
  #   for(j in 1:length(c)){
  #     invediff[i,j] = c[j] - r[i]
  #   }
  # }
  # #invediff = repmat(inv_lambda, 1, rank_of_data) - t(repmat(inv_lambda, 1, rank_of_data))
  # invediff[invediff <= 0] = 1
  #invediff = log(invediff)

  cumsum_invediff = apply(invediff, 2, cumsum)
  row_invediff = apply(cumsum_invediff, 1, sum)

  log_lambda = log(lambda_hat)
  cumsum_loglambda = cumsum(log_lambda)

  cumsum_lambdahat = cumsum(lambda_hat)

  kmax = rank_of_data - 1
  ks = c(1:kmax)

  p = rep(NA, length(ks)) #likelihood

  z = log(2) + (d-ks+1)/2*log(pi) - lgamma((d - ks + 1)/2)
  cumsum_z = cumsum(z)

  for(i in 1:length(ks)){
    k = ks[i]
    v = (cumsum_lambdahat[length(cumsum_lambdahat)] - cumsum_lambdahat[k])/(d-k)
    p[i] = -n/2*cumsum_loglambda[k] + (-n*(d-k)/2)*log(v)
    p[i] = p[i] - cumsum_z[k] - k/2*log(n)
    #h = logdet(A_Z)
    h = row_invediff[k] + cumsum_logediff[k]
    h = h + (d-k)*sum(log(1/v - inv_lambda[1:k]))
    m = d*k-k*(k+1)/2
    h = h + m*log(N)
    p[i] = p[i] + (m+k)/2*log(2*pi) - h/2
    p[i] = p[i] + 1.5*k*log(2)
    p[i] = p[i] - 0.5*log(d-k)
    if (alpha > 0){
      ck = alpha*(d-k)/2*log(beta*(d-k)/2) - lgamma(alpha*(d-k)/2) + k*(alpha/2*log(beta/2) - lgamma(alpha/2))
      p[i] = p[i] - n*d/2 - 0.5*log(n) + ck
    }
  }
  pmax = max(p, na.rm=TRUE)
  i = which.max(p)
  k = ks[i]

  #compare with p0
  v0 = cumsum_lambdahat[length(cumsum_lambdahat)]/length(cumsum_lambdahat)
  p0 = -n*d/2*log(v0) - 0.5*log(d)

  if(alpha > 0){
    p0 = p0 - n*d/2 - 0.5*log(n) + alpha*d/2*log(beta*d/2) - lgamma(alpha*d/2)
  }
  if(p0 >= pmax){
    k = 0
  }

  #get error variance
  v_hat = (cumsum_lambdahat[length(cumsum_lambdahat)] - cumsum_lambdahat[k])/(d-k)

  return(list(k, v_hat, p))
}

