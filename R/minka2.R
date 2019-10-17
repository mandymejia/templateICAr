# minka = function(X,alpha=0,beta=0,eps = 1e-12,rwr = FALSE){
#   #rwr = replace with rank
#   #defaulting alpha and beta to 0 can be dangerous. please don't do it, please
#   #probably better to default to 1
#   if(alpha<=0){alpha=0;beta=0}
#
#   #X is a matrix with observations in rows
#   N = dim(X)[1]
#   d = dim(X)[2]
#
#   #Check whether full rank (if d < N), return error otherwise?
#   #If X has already been centered, it will have rank d-1
#
#   #center X across observations
#   #X <- sweep(X,2,apply(X,2,mean),"-"))/sqrt(N))
#   X <- (X - matrix(colMeans(X), nrow=N, ncol=d, byrow=TRUE))/sqrt(N) #using matrix operations (fastest)
#
#   #perform SVD on XXt or XtX (whichever smaller)
#   if(N < d) {XX = tcrossprod(X)} else {XX=crossprod(X)} #transpose if necessarily ill formed
#   svd_XX = svd(XX)
#   lambda = svd_XX$d
#
#   #augment with zero eigenvalues if fewer than d estimated with SVD (N < d)
#   if(length(lambda) < d){
#     lambda = c(lambda,rep(0,d-length(lambda)))
#   }
#   lambda_min_zero = min(c(Inf,which(lambda/lambda[1]<=eps))) #index of first zero eigenvalue
#   #print(c(lambda_min_zero,length(lambda)))
#   #set eigenvalues below eps to zero
#   if(lambda_min_zero < Inf){
#     lambda[lambda_min_zero:length(lambda)] = rep(0, length(lambda) - lambda_min_zero + 1)
#     rankX = lambda_min_zero-1
#   }else{
#     rankX = length(lambda)
#   }
#
#   #replace d with rank, this should be done for model to make sense if X is rank deficient (not it N<d from sample size)
#   #replace_d_with_rankX = TRUE
#   #if(replace_d_with_rankX==rwr){D=rankX}else{D=d}
#   if(rwr){D=rankX}else{D=d}
#   #print(c(D,rankX))
#   #we will use D where using rankX could make more sense than using dimX
#   maxK = rankX-1 #rankX  - 1 #??? - estimable things?
#   lambda = lambda[1:D]
#   #k values
#   k = c(1:maxK)
#   #do everything in the log of the densities
#   #try to do all values of k at once when possible
#   #n = N+alpha-1
#   n = N+ifelse(alpha==0,0,alpha-1)
#
#   #mult correction
#   multiplicity_correction = k*log(2) #in front of c_k
#   #print("(multiplicity_correction)")
#   #print((multiplicity_correction))
#   #power of 2pi on density for x
#   pow_2pi_likelihood = -(N*D)/2*log(2*pi) #part of c_k
#   #print("(pow_2pi_likelihood)")
#   #print((pow_2pi_likelihood))
#   #N^(-d/2) from integrating out the m
#   m_int_const = -D/2*log(N)+D/2*log(2*pi) #part of c_k
#   #print("(m_int_const)")
#   #print((m_int_const))
#   #prior const for U
#   prior_const_U = -k*log(2)+cumsum(lgamma((D-k+1)/2)-(D-k+1/2)*log(pi)) #part of c_k - think about this
#   #print("(prior_const_U)")
#   #print((prior_const_U))
#     #cumsum is doing the sum (product) over the i=1:k
#     #the d is the ambient dimension if X is rank deficient, not the rank
#   #prior const for v and L
#   prior_const_L = if(alpha==0){1}else{(alpha/2*log(beta/2)-lgamma(alpha/2))*k}#part of c_k
#   #print("(prior_const_L)")
#   #print((prior_const_L))
#   prior_const_v = if(alpha==0){1}else{alpha*(D-k)/2*log(beta*(D-k)/2)-lgamma(alpha*(D-k)/2)}#part of c_k
#   #print("(prior_const_L)")
#   #print((prior_const_v))
#   #det L part
#   L_hat = 1/n*(N*lambda+beta)
#   det_L_part = -n/2*cumsum(log(L_hat[k]))
#   #print("(det_L_part)")
#   #print((det_L_part))
#   #v_hat_part
#   v_hat = 1/(n*(D-k))*N*(cumsum(lambda[rankX:2])+ifelse(rankX<D,sum(lambda[(rankX+1):D]),0))[maxK:1]+beta/n
#   #print(v_hat)
#   v_hat_part = -n*(D-k)/2*log(v_hat) #notsure
#   #print("(v_hat_part)")
#   #print((v_hat_part))
#   #exponential part
#   exp_part = -n*D/2 #notsure
#   #print("(exp_part)")
#   #print((exp_part))
#   #2pi power k = dim lambda, 1 = dim v, d*k-k-choose(k,2)  = dim U
#   two_pi_part_posterior_approx = (k+1+D*k-k-choose(k,2))/2*log(2*pi)
#   #print("(two_pi_part_posterior_approx)")
#   #print((two_pi_part_posterior_approx))
#   #AL part
#   AL_part = -k/2*log(n/2)
#   #print("(AL_part)")
#   #print((AL_part))
#   #Av part
#   Av_part = -1/2*log(n*(D-k)/2)
#   #print("(Av_part)")
#   #print((Av_part))
#   #Az part
#   #part 1 is the difference in lambda part
#   Az_part_1 = rep(0,maxK)
#   for(i in k){
#     j_index = (i+1):D
#     Az_part_1[i] = sum(log(lambda[i]-lambda[j_index]))
#   }
#   Az_part_1 = -1/2*cumsum(Az_part_1)
#   #print("(Az_part_1)")
#   #print((Az_part_1))
#   #part 2 is for the inverse of the lambda hats
#   #see eq 70 in updated paper
#   Az_part_2 = rep(0,maxK)
#   lambda_hat_leqk = L_hat #note, name of variable applies when the index is <= k for a given k
#   log_diff_inv_lambda_hat_leqk = rep(0,maxK)
#   for(j in (k+1)){
#     ##print(c(j,j-1))
#     i_index = 1:(j-1)
#     log_diff_inv_lambda_hat_leqk[j-1] = sum(log(lambda_hat_leqk[j]^(-1)-lambda_hat_leqk[i_index]^(-1)))
#   }
#   Az_part_2 = -1/2*cumsum(log_diff_inv_lambda_hat_leqk)
#   #print("(Az_part_2)")
#   #print((Az_part_2))
#   for(ell in k){
#     Az_part_2[ell] = Az_part_2[ell] - (D-ell)/2*sum(log(v_hat[ell]^(-1)-lambda_hat_leqk[1:ell]^(-1)))
#   }
#   #print("(Az_part_2)")
#   #print((Az_part_2))
#   Az_part_3 = -1/2*(D*k-(k*(k+1)/2))*log(N)
#   #print("(Az_part_3)")
#   #print((Az_part_3))
#
#
#   pow_2pi_likelihood_p0 = -(N*D)/2*log(2*pi)
#   m_int_const_p0 = -D/2*log(N)+D/2*log(2*pi)
#   prior_const_v_p0 = if(alpha==0){1}else{alpha*D/2*log(beta*D/2)-lgamma(alpha*D/2)}
#   v_hat_p0 = 1/(n*D)*(sum(N*lambda+beta))
#   v_hat_part_p0 = -n*D/2*log(v_hat_p0)
#   exp_part_p0 = -n*D/2
#   two_pi_part_posterior_approx_p0 = 1/2*log(2*pi)
#   Av_part_p0 = -1/2*log(n*D/2)
#   p0 = pow_2pi_likelihood_p0+m_int_const_p0+prior_const_v_p0+v_hat_part_p0+exp_part_p0+two_pi_part_posterior_approx_p0+Av_part_p0
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
#   max_estimable = maxK#max(which(is.finite(v_hat_part)))
#   p = c(p0,p[1:max_estimable])
#   names(p) = paste("rank",c(0:max_estimable),sep="=")
#   p[p==Inf] = -Inf
#   v_hat = c(mean(lambda),v_hat[1:max_estimable])
#   names(v_hat) = names(p)
#   return(
#     list(p=p,v_hat = v_hat, lambda = lambda)
#   )
# }
#
# plot.minka <- function(minka){
#
# }
#
# # N = 3000
# # D = 20
# # test_low = TRUE #whether to add collinear columns to X
# # M = 500 # # of extra columns to add
# # mu = rt(D,3)*4
# # k = 13
# # sd=2
# # H = matrix(rnorm(D*k),D,k)
# # w = matrix(rnorm(N*k),k,N)
# # lrm = t(H%*%w)
# # X = matrix(rnorm(N*D,sd=sd),N,D)+matrix(mu,N,D,byrow=TRUE)+lrm
# # if(test_low){A = matrix(rnorm(D*D),D,D); B = cbind(A,A[,sample(D,M,replace=TRUE)]+A[,sample(D,M,replace=TRUE)]);XX=cbind(X%*%A,X%*%B);X=X%*%A}else{XX=X;}
# # svd(X)$d
# # svd(XX)$d
# # out = minka(XX,0,0,1e-12,rwr=F)
# probs = out$p-max(out$p)
# probs = exp(probs)
# probs = probs/sum(probs)
# probs[probs==max(probs)]
# plot(log(probs))
# pp = sort(probs,decreasing = TRUE, index.return=TRUE)
# round(pp$x,4)[1:20]
# plot(out$v_hat)
# round(out$v_hat[pp$ix],4)[1:20]
