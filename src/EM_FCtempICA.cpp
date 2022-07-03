#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

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

// [[Rcpp::export]]
Rcpp::List Gibbs_AS_posteriorCPP(const int nsamp, const int nburn,
                                 const Eigen::MatrixXd template_mean,
                                 const Eigen::MatrixXd template_var,
                                 Eigen::MatrixXd S, Eigen::MatrixXd A,
                                 const Eigen::MatrixXd G,
                                 const Eigen::VectorXd tau_v,
                                 const Eigen::MatrixXd Y,
                                 const Eigen::VectorXd alpha,
                                 bool final) {
  // Find dimension quantities
  int niter = nsamp + nburn;
  int V = tau_v.size();
  int Q = template_mean.cols();
  int ntime = Y.cols();
  // Initialize output values
  Eigen::VectorXd A_sum = Eigen::VectorXd::Zero(Q);
  Eigen::MatrixXd AtA_sum = Eigen::MatrixXd::Zero(Q,Q);
  Eigen::VectorXd yAS_sum = Eigen::VectorXd::Zero(V);
  Eigen::VectorXd AS_sq = Eigen::VectorXd::Zero(V);
  // Calculate static quantities
  Eigen::VectorXd tau_inv = tau_v.cwiseInverse();
  Eigen::MatrixXd G_tau_inv = tau_inv.asDiagonal();
  Eigen::MatrixXd G_inv = G.inverse();
  Eigen::VectorXd alphaGinv = G_inv * alpha;
  Eigen::MatrixXd YtG_tau_inv = Y.transpose() * G_tau_inv;
  // Initialize intermediate objects
  Eigen::MatrixXd sig_inv_A = Eigen::MatrixXd::Zero(Q,Q);
  Eigen::MatrixXd chol_sig_A(Q,Q), chol_sig_S(Q,Q);
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
  for(int i=1;i <= niter; i++) {
    // Update A
    // Rcout << "Update A" << std::endl;
    SGti = S.transpose() * G_tau_inv;
    sig_inv_A = SGti * S;
    sig_inv_A += G_inv;
    chol_sig_inv_A.compute(sig_inv_A);
    chol_sig_A = chol_sig_inv_A.matrixL().toDenseMatrix().inverse();
    YGS = YtG_tau_inv * S;
    for(int t = 1;t <= ntime; t++) {
      ygs_alphaGinv = YGS.row(t-1) + alphaGinv.transpose();
      mu_at = chol_sig_inv_A.solve(ygs_alphaGinv);
      Z = rnorm(Q);
      Eigen::Map<Eigen::VectorXd> ZZ = as<Eigen::Map<Eigen::VectorXd> >(Z);
      At = chol_sig_A * ZZ;
      At += mu_at;
      A.row(t-1) = At;
    }
    // Update S
    // Rcout << "Update A" << std::endl;
    AtA = A.transpose() * A;
    for(int v=1; v <= V; v++) {
      sig_inv_S = AtA;
      sig_inv_S /= tau_v(v-1);
      // sig_inv_S += template_var.row(v-1).asDiagonal();
      sig_inv_S += G_sv_inv.row(v-1).asDiagonal();
      chol_sig_inv_S.compute(sig_inv_S);
      chol_sig_S = chol_sig_inv_S.matrixL().toDenseMatrix().inverse();
      AtYvtempVarMean = A.transpose() * Y.row(v-1).transpose();
      AtYvtempVarMean /= tau_v(v-1);
      tVarMean = template_var.row(v-1).asDiagonal() * template_mean.row(v-1).transpose();
      AtYvtempVarMean += tVarMean;
      mu_sv = chol_sig_inv_S.solve(AtYvtempVarMean);
      Z = rnorm(Q);
      Eigen::Map<Eigen::VectorXd> ZZ = as<Eigen::Map<Eigen::VectorXd> >(Z);
      Sv = chol_sig_S * ZZ;
      Sv += mu_sv;
      S.row(v-1) = Sv;
    }
    // Store summaries
    if(i > nburn) {
      AtS = A * S.transpose();
      // AtA = A * A.transpose();
      YtAtS = Y.transpose().cwiseProduct(AtS);
      for(int t = 0; t < ntime; t++) {
        A_sum += A.row(t);
        yAS_sum += YtAtS.row(t);
        AS_sq += AtS.row(t).cwiseAbs2();
      }
      AtA_sum += AtA;
    }
  }
  A_sum /= nsamp;
  AtA_sum /= nsamp;
  yAS_sum /= nsamp;
  AS_sq /= nsamp;
  Rcpp::NumericVector A_sumX(wrap(A_sum));
  SEXP aat = Rcpp::wrap(AtA_sum);
  Rcpp::NumericMatrix AtA_sumX(aat);
  Rcpp::NumericVector yAS_sumX(wrap(yAS_sum));
  Rcpp::NumericVector AS_sqX(wrap(AS_sq));
  return Rcpp::List::create(Rcpp::Named("A_sum") = A_sumX,
                            Rcpp::Named("AtA_sum") = AtA_sumX,
                            Rcpp::Named("yAS_sum") = yAS_sumX,
                            Rcpp::Named("AS_sq_sum") = AS_sqX);
}
