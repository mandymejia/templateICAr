#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

//' Update parameters (M-step of the EM)
//'
//' @param template_mean a matrix with dimensions V x Q giving the mean value
//'   of the independent components
//' @param template_var a matrix with dimensions V x Q giving the variance
//'   of the independent components
//' @param template_FC a list with two elements: psi, a Q x Q matrix, and nu, a
//'   scalar. These two values are the parameters of the Wishart prior on G
//' @param G a Q x Q matrix of the prior covariance of A
//' @param prior_params a length 2 vector with the prior parameters for tau_v
//' @param BOLD a V x T matrix of BOLD values
//' @param Y_sq_sum a length V vector with the sum of squared BOLD values at
//'   each data location
//' @param post_sums a list of posterior summary statistics including the named
//'   summaries \code{AS_sq_sum}, \code{yAS_sum}, \code{A_sum}, and \code{AtA_sum}
//' @param sigma2_alpha a scalar multiplier for the prior variance of alpha
//' @param verbose a boolean. Should messages be generated and output?
//' @return A list with quantities tau_sq, alpha, and G
//' @export
// [[Rcpp::export]]
Rcpp::List UpdateTheta_FCtemplateICAcpp(Eigen::MatrixXd template_mean,
                                        Eigen::MatrixXd template_var,
                                        Rcpp::List template_FC,
                                        Eigen::MatrixXd G,
                                        Eigen::VectorXd prior_params,
                                        Eigen::MatrixXd BOLD,
                                        Eigen::VectorXd Y_sq_sum,
                                        Rcpp::List post_sums,
                                        double sigma2_alpha,
                                        bool verbose){
  // Bring in data and initialize quantities
  int nvox = BOLD.rows();
  int ntime = BOLD.cols();
  int nICs = template_mean.cols();
  Eigen::VectorXd tau_sq_new(nvox);
  double alpha_tau = prior_params(0);
  double beta_tau = prior_params(1);
  double tau_sq_num, tau_sq_den;
  Eigen::Map<Eigen::VectorXd> AS_sq(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(post_sums["AS_sq_sum"]));
  Eigen::Map<Eigen::VectorXd> yAS(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(post_sums["yAS_sum"]));
  Eigen::Map<Eigen::VectorXd> A(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(post_sums["A_sum"]));
  Eigen::Map<Eigen::MatrixXd> AtA(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(post_sums["AtA_sum"]));
  Eigen::Map<Eigen::MatrixXd> psi0(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(template_FC["psi"]));
  // Eigen::VectorXd A = Eigen::VectorXd (post_sums["A_sum"]);
  // Eigen::MatrixXd AtA = Eigen::MatrixXd (post_sums["AtA_sum"]);
  Eigen::VectorXd alpha_new(nICs);
  double nu0 = double (template_FC["nu"]);
  // Eigen::MatrixXd psi0 = Eigen::MatrixXd (template_FC["psi"]);
  // UPDATE TAU^2 (ERROR VARIANCE)
  tau_sq_den = ntime + 2*alpha_tau;
  tau_sq_den += 2;
  for(int v=0;v < nvox; v++){
    tau_sq_num = AS_sq(v) - 2*yAS(v);
    tau_sq_num += Y_sq_sum(v);
    tau_sq_num += 2*beta_tau;
    tau_sq_new(v) = tau_sq_num / tau_sq_den;
    // tau_sq_new(v) = 1/(ntime + 2*alpha_tau + 2) * (AS_sq(v) - 2*yAS(v) + 2*beta_tau);
  }
  double tau_mean = tau_sq_new.mean();
  for(int v=0;v<nvox;v++){
    tau_sq_new(v) = tau_mean;
  }
  // UPDATE ALPHA (TEMPORAL INTERCEPT)
  Eigen::MatrixXd alpha_part = Eigen::MatrixXd::Identity(nICs,nICs);
  Eigen::MatrixXd G_inv = G.inverse();
  alpha_part /= sigma2_alpha;
  alpha_part += ntime * G_inv;
  Eigen::MatrixXd alpha_part_inv = alpha_part.inverse();
  alpha_part_inv *= G_inv;
  alpha_new = alpha_part_inv * A;
  // UPDATE G (TEMPORAL COVARIANCE, I.E. FUNCTIONAL CONNECTIVITY)
  Eigen::MatrixXd alpha_alpha_t = alpha_new * alpha_new.transpose();
  Eigen::MatrixXd tmp = AtA - 2 * A * alpha_new.transpose() + ntime*alpha_alpha_t + psi0;
  Eigen::MatrixXd G_new = tmp / (ntime + nu0 + nICs + 1);
  // RETURN NEW PARAMETER ESTIMATES
  Rcpp::NumericVector tau_sq_newX(wrap(tau_sq_new));
  Rcpp::NumericVector alpha_newX(wrap(alpha_new));
  Rcpp::NumericMatrix G_newX(wrap(G_new));
  return Rcpp::List::create(Rcpp::Named("tau_sq") = tau_sq_newX,
                            Rcpp::Named("alpha") = alpha_newX,
                            Rcpp::Named("G") = G_newX);
}

//' Use a Gibbs sampler for the A and S variables (E-step of the EM)
//'
//' @param nsamp the number of posterior samples to output after burn-in
//' @param nburn the number of posterior samples to throw away before saving
//' @param template_mean a matrix with dimensions V x Q giving the mean value
//'   of the independent components
//' @param template_var a matrix with dimensions V x Q giving the variance
//'   of the independent components
//' @param S a matrix with dimensions V x Q of subject independent components
//' @param G a Q x Q matrix of the prior covariance of A
//' @param tau_v a length V vector with noise variance for each data location
//' @param Y a matrix with dimensions V x T of observed BOLD data
//' @param alpha a length Q vector of the prior mean of all rows of A
//' @param final a boolean. Should posterior samples be returned instead of
//'   summary measures?
//' @param return_samp a boolean. Should posterior samples be returned?
//' @return List with estimates for A, S, and possibly other quantities
//' @export
// [[Rcpp::export]]
Rcpp::List Gibbs_AS_posteriorCPP(const int nsamp, const int nburn,
                                 const Eigen::MatrixXd template_mean,
                                 const Eigen::MatrixXd template_var,
                                 Eigen::MatrixXd S,
                                 const Eigen::MatrixXd G,
                                 const Eigen::VectorXd tau_v,
                                 const Eigen::MatrixXd Y,
                                 const Eigen::VectorXd alpha,
                                 bool final,
                                 bool return_samp) {
  // Find dimension quantities
  int niter = nsamp + nburn;
  int V = tau_v.size();
  int Q = template_mean.cols();
  int ntime = Y.cols();
  // This is just to reduce confusion about the input dimensions of S
  if(S.cols() == V) S = S.transpose();
  // Initialize output values
  Rcpp::List output;
  int TQ = ntime * Q;
  int VQ = V * Q;
  Eigen::MatrixXd A_samp, S_samp, AtA_sum, S_post = Eigen::MatrixXd::Zero(V,Q);
  Eigen::VectorXd A_sum, yAS_sum, AS_sq;
  if(final) return_samp = true;
  if(return_samp){
    // In the case of final, just output posterior samples of A and S
    A_samp = Eigen::MatrixXd::Zero(TQ,nsamp);
    S_samp = Eigen::MatrixXd::Zero(VQ,nsamp);
  }
  if(!final){
    // In the case of not final, output posterior summaries
    A_sum = Eigen::VectorXd::Zero(Q);
    AtA_sum = Eigen::MatrixXd::Zero(Q,Q);
    yAS_sum = Eigen::VectorXd::Zero(V);
    AS_sq = Eigen::VectorXd::Zero(V);
  }
  // Calculate static quantities
  double mean_tau = tau_v.mean();
  Eigen::VectorXd tau_inv(V);
  for(int v=0;v<V;v++) tau_inv(v) = 1 / mean_tau;
  // Eigen::VectorXd tau_inv = tau_v.cwiseInverse();
  Eigen::MatrixXd G_tau_inv = tau_inv.asDiagonal();
  Eigen::MatrixXd G_inv = G.inverse();
  Eigen::VectorXd alphaGinv = G_inv * alpha;
  Eigen::MatrixXd YtG_tau_inv = Y.transpose() * G_tau_inv;
  // Initialize intermediate objects
  Eigen::MatrixXd sig_inv_A = Eigen::MatrixXd::Zero(Q,Q);
  Eigen::MatrixXd chol_sig_A(Q,Q), chol_sig_S(Q,Q), A(ntime,Q);
  Eigen::MatrixXd YGS(ntime,Q), AtA(Q,Q), sig_inv_S(Q,Q), SGti(Q,V);
  Eigen::MatrixXd AtS(ntime, V), YtAtS(ntime,V);
  Eigen::VectorXd mu_at(Q), mu_sv(Q), tVarRow(Q), At(Q), Sv(Q);
  Eigen::VectorXd ygs_alphaGinv(Q), AtYvtempVarMean(Q), tVarMean(Q);
  Eigen::MatrixXd G_sv_inv(V,Q);
  // Initialize G_sv_inv so it doesn't happen in the loop
  for(int v=0;v<V;v++) {
    for(int q=0;q<Q;q++) {
      G_sv_inv(v,q) = 1 / template_var(v,q);
    }
  }
  NumericVector Z(Q);
  Eigen::LLT<Eigen::MatrixXd> chol_sig_inv_A, chol_sig_inv_S;
  // Start the Gibbs sampler
  for(int i=0;i < niter; i++) {
    // Update A
    SGti = S.transpose() * G_tau_inv;
    sig_inv_A = SGti * S;
    sig_inv_A += G_inv;
    chol_sig_inv_A.compute(sig_inv_A);
    chol_sig_A = chol_sig_inv_A.matrixL().toDenseMatrix().inverse();
    YGS = YtG_tau_inv * S;
    for(int t = 0;t < ntime; t++) {
      ygs_alphaGinv = YGS.row(t) + alphaGinv.transpose();
      mu_at = chol_sig_inv_A.solve(ygs_alphaGinv);
      Z = rnorm(Q);
      Eigen::Map<Eigen::VectorXd> ZZ = as<Eigen::Map<Eigen::VectorXd> >(Z);
      At = chol_sig_A * ZZ;
      At += mu_at;
      A.row(t) = At;
    }
    // Update S
    AtA = A.transpose() * A;
    for(int v=0; v < V; v++) {
      sig_inv_S = AtA;
      sig_inv_S /= tau_v(v);
      sig_inv_S += G_sv_inv.row(v).asDiagonal();
      chol_sig_inv_S.compute(sig_inv_S);
      chol_sig_S = chol_sig_inv_S.matrixL().toDenseMatrix().inverse();
      AtYvtempVarMean = A.transpose() * Y.row(v).transpose();
      AtYvtempVarMean /= tau_v(v);
      tVarMean = G_sv_inv.row(v).asDiagonal() * template_mean.row(v).transpose();
      AtYvtempVarMean += tVarMean;
      mu_sv = chol_sig_inv_S.solve(AtYvtempVarMean);
      Z = rnorm(Q);
      Eigen::Map<Eigen::VectorXd> ZZ = as<Eigen::Map<Eigen::VectorXd> >(Z);
      Sv = chol_sig_S * ZZ;
      Sv += mu_sv;
      S.row(v) = Sv;
    }
    // Store summaries
    if(i > (nburn-1)) {
      if(!final) {
        AtS = A * S.transpose();
        S_post += S;
        // AtA = A * A.transpose();
        YtAtS = Y.transpose().cwiseProduct(AtS);
        for(int t = 0; t < ntime; t++) {
          A_sum += A.row(t);
          yAS_sum += YtAtS.row(t);
          AS_sq += AtS.row(t).cwiseAbs2();
        }
        AtA_sum += AtA;
      }
      if(return_samp) {
        for(int t = 0;t<ntime;t++){
          for(int q=0;q<Q;q++){
            int Aidx = t + q*ntime;
            A_samp(Aidx,i - nburn) = A(t,q);
          }
        }
        for(int v=0;v<V;v++) {
          for(int q=0;q<Q;q++) {
            int Sidx = v + q*V;
            S_samp(Sidx,i - nburn) = S(v,q);
          }
        }
      }
    }
  }
  if(!final) {
    A_sum /= nsamp;
    AtA_sum /= nsamp;
    yAS_sum /= nsamp;
    AS_sq /= nsamp;
    S_post /= nsamp;
    Rcpp::NumericVector A_sumX(wrap(A_sum));
    SEXP aat = Rcpp::wrap(AtA_sum);
    SEXP s_post = Rcpp::wrap(S_post);
    Rcpp::NumericMatrix AtA_sumX(aat);
    Rcpp::NumericVector yAS_sumX(wrap(yAS_sum));
    Rcpp::NumericVector AS_sqX(wrap(AS_sq));
    Rcpp::NumericMatrix S_postX(s_post);
    if(!return_samp) {
      output = Rcpp::List::create(Rcpp::Named("A_sum") = A_sumX,
                                  Rcpp::Named("AtA_sum") = AtA_sumX,
                                  Rcpp::Named("yAS_sum") = yAS_sumX,
                                  Rcpp::Named("AS_sq_sum") = AS_sqX,
                                  Rcpp::Named("S_post") = S_post);
    } else {
      output = Rcpp::List::create(Rcpp::Named("A_sum") = A_sumX,
                                  Rcpp::Named("AtA_sum") = AtA_sumX,
                                  Rcpp::Named("yAS_sum") = yAS_sumX,
                                  Rcpp::Named("AS_sq_sum") = AS_sqX,
                                  Rcpp::Named("S_post") = S_post,
                                  Rcpp::Named("A_samp") = A_samp,
                                  Rcpp::Named("S_samp") = S_samp);
    }
  }
  if(final) {
    output = Rcpp::List::create(Rcpp::Named("A_samp") = A_samp,
                                Rcpp::Named("S_samp") = S_samp);
  }

  return output;
}
