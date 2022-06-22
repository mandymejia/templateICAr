#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Rcpp::List UpdateTheta_FCtemplateICAcpp(Eigen::MatrixXd template_mean,
//                                         Eigen::MatrixXd template_var,
//                                         Rcpp::List template_FC,
//                                         Eigen::VectorXd prior_params,
//                                         Eigen::MatrixXd BOLD,
//                                         Rcpp::List post_sums,
//                                         bool verbose){
//
//
//   int nvox = BOLD.rows();
//   int ntime = BOLD.cols();
//   int nICs = template_mean.cols();
//
//   // UPDATE TAU^2 (ERROR VARIANCE)
//
//   Eigen::VectorXd tau_sq_new(nvox);
//   double alpha_tau = prior_params(0);
//   double beta_tau = prior_params(1);
//   Eigen::VectorXd AS_sq = Eigen::VectorXd (post_sums["AS_sq"]);
//   Eigen::VectorXd yAS = Eigen::VectorXd (post_sums["yAS_sum"]);
//   Eigen::VectorXd A = Eigen::VectorXd (post_sums["A_sum"]);
//   Eigen::VectorXd AAt = Eigen::VectorXd (post_sums["AAt_sum"]);
//   Eigen::VectorXd alpha_new(nICs);
//   double nu0 = double (template_FC["nu"]);
//   Eigen::MatrixXd psi0 = Eigen::MatrixXd (template_FC["psi"]);
//
//   for(int v=0;v < nvox; v++){
//     tau_sq_new(v) = 1/(ntime + 2*alpha_tau + 2) * (AS_sq(v) - 2*yAS(v) + 2*beta_tau);
//   }
//
// // UPDATE ALPHA (TEMPORAL INTERCEPT)
//
//   alpha_new = A / ntime;
//
//
// // UPDATE G (TEMPORAL COVARIANCE, I.E. FUNCTIONAL CONNECTIVITY)
//
//       Eigen::MatrixXd alpha_alpha_t = alpha_new * alpha_new.transpose();
//       Eigen::MatrixXd tmp = 2* A * alpha_new.transpose() - AAt - ntime*alpha_alpha_t - psi0;
//       Eigen::MatrixXd G_new = (ntime + nu0 + nICs + 1) * tmp.inverse();
//
// // RETURN NEW PARAMETER ESTIMATES
//       // theta_new <- list(
//       //     tau_sq <- tau_sq_new,
//       //             alpha <- alpha_new,
//       //                    G <- G_new
//       // )
//       //   return(theta_new)
//
//       Rcpp::NumericVector tau_sq_newX(wrap(tau_sq_new));
//       Rcpp::NumericVector alpha_newX(wrap(alpha_new));
//       Rcpp::NumericVector G_newX(wrap(G_new));
//       return Rcpp::List::create(Rcpp::Named("tau_sq") = tau_sq_newX,
//                                 Rcpp::Named("alpha") = alpha_newX,
//                                 Rcpp::Named("G") = G_newX);
// }

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
  // Rcout << "Dimension check: niter = " << niter << ", V = " << V << ", Q = " << Q << ", ntime = " << ntime << std::endl;
  // Initialize output values
  Eigen::VectorXd A_sum = Eigen::VectorXd::Zero(Q);
  Eigen::VectorXd AAt_sum = Eigen::VectorXd::Zero(ntime);
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
  Eigen::MatrixXd AAt(ntime,ntime), AtS(ntime, V), YtAtS(ntime,V);
  Eigen::VectorXd mu_at(Q), mu_sv(Q), tVarRow(Q), At(Q), Sv(Q);
  Eigen::VectorXd ygs_alphaGinv(Q), AtYvtempVarMean(Q), tVarMean(Q);
  NumericVector Z(Q);
  Eigen::LLT<Eigen::MatrixXd> chol_sig_inv_A, chol_sig_inv_S;
  // Start the Gibbs sampler
  for(int i=1;i <= niter; i++) {
    // Rcout << "Iteration " << i << std::endl;
    // Update A
    SGti = S.transpose() * G_tau_inv;
    sig_inv_A = SGti * S;
    sig_inv_A += G_inv;
    // sig_A = sig_inv_A.inverse();
    // chol_sig_A.compute(sig_A);
    // Rcout << sig_A << std::endl;
    // Eigen::LLT<Eigen::MatrixXd> chol_sig_inv_A(sig_inv_A);
    chol_sig_inv_A.compute(sig_inv_A);
    chol_sig_A = chol_sig_inv_A.matrixL().toDenseMatrix().inverse();
    YGS = YtG_tau_inv * S;
    // ygs_alphaGinv = YGS.row(ntime - 1) + alphaGinv.transpose();
    // mu_at = sig_A * ygs_alphaGinv;
    // Rcout << mu_at.transpose() << std::endl;
    // Rcout << YGS.row(ntime - 1) << std::endl;
    // Rcout << "Updating A" << std::endl;
    for(int t = 1;t <= ntime; t++) {
      // Eigen::VectorXd ygs_alphaGinv = Eigen::VectorXd::Zero(Q);
      ygs_alphaGinv = YGS.row(t-1) + alphaGinv.transpose();
      // mu_a.row(t) = chol_sig_inv_A.solve(ygs_alphaGinv);
      mu_at = chol_sig_inv_A.solve(ygs_alphaGinv);
      // mu_at = sig_A * ygs_alphaGinv;
      // if(t == 1) Rcout << "mu_a1: " << mu_at.transpose() << std::endl;
      // mu_at = sig_inv_A.HouseholderQR().solve(ygs_alphaGinv);
      Z = rnorm(Q);
      Eigen::Map<Eigen::VectorXd> ZZ = as<Eigen::Map<Eigen::VectorXd> >(Z);
      // At = chol_sig_inv_A.solve(ZZ);
      // At = chol_sig_A.matrixL() * ZZ;
      // At = chol_sig_inv_A.matrixL().toDenseMatrix().inverse() * ZZ;
      At = chol_sig_A * ZZ;
      At += mu_at;
      // A.row(t-1) = chol_sig_inv_A.solve(ZZ);
      // A.row(t-1) += mu_at;
      A.row(t-1) = At;
    }
    // Rcout << ygs_alphaGinv.transpose() << std::endl;
    // Rcout << mu_at.transpose() << std::endl;
    // Update S
    AtA = A.transpose() * A;
    // Rcout << "Updating S" << std::endl;;
    for(int v=1; v <= V; v++) {
      // Rcout << "Vertex " << v << std::endl;
      sig_inv_S = AtA;
      sig_inv_S /= tau_v(v-1);
      sig_inv_S += template_var.row(v-1).asDiagonal();
      // sig_S = sig_inv_S.inverse();
      // chol_sig_S.compute(sig_S);
      // Rcout << "Found sig_inv_S...";
      chol_sig_inv_S.compute(sig_inv_S);
      chol_sig_S = chol_sig_inv_S.matrixL().toDenseMatrix().inverse();
      // Rcout << "Cholesky computed...";
      AtYvtempVarMean = A.transpose() * Y.row(v-1).transpose();
      AtYvtempVarMean /= tau_v(v-1);
      tVarMean = template_var.row(v-1).asDiagonal() * template_mean.row(v-1).transpose();
      AtYvtempVarMean += tVarMean;
      // Rcout << "RHS computed...";
      mu_sv = chol_sig_inv_S.solve(AtYvtempVarMean);
      // mu_sv = sig_S * AtYvtempVarMean;
      // if(v == 1) Rcout << "mu_s1: " << mu_sv.transpose() << std::endl;
      // Rcout << "mu computed...";
      Z = rnorm(Q);
      Eigen::Map<Eigen::VectorXd> ZZ = as<Eigen::Map<Eigen::VectorXd> >(Z);
      // S.row(v-1) = chol_sig_inv_S.solve(ZZ);
      // S.row(v-1) += mu_sv;
      Sv = chol_sig_S * ZZ;
      Sv += mu_sv;
      S.row(v-1) = Sv;
    }
    // Rcout << mu_sv.transpose() << std::endl;
    // Rcout << A_sum.transpose() << std::endl;
    // Store summaries
    if(i > nburn) {
      AtS = A * S.transpose();
      AAt = A * A.transpose();
      YtAtS = Y.transpose().cwiseProduct(AtS);
      for(int t = 0; t < ntime; t++) {
        A_sum += A.row(t);
        AAt_sum += AAt.row(t);
        yAS_sum += YtAtS.row(t);
        AS_sq += AtS.row(t).cwiseAbs2();
      }
    }
  }
  A_sum = A_sum / nsamp;
  AAt_sum = AAt_sum / nsamp;
  yAS_sum = yAS_sum / nsamp;
  AS_sq = AS_sq / nsamp;
  // Rcout << "A_sum:" << std::endl;
  // Rcout << A_sum.transpose() << std::endl;
  // Rcout << "AAt_sum:" << std::endl;
  // Rcout << AAt_sum.transpose() << std::endl;
  // Rcout << "yAS_sum:" << std::endl;
  // Rcout << yAS_sum.transpose() << std::endl;
  // Rcout << "AS_sq:" << std::endl;
  // Rcout << AS_sq.transpose() << std::endl;
  Rcpp::NumericVector A_sumX(wrap(A_sum));
  Rcpp::NumericVector AAt_sumX(wrap(AAt_sum));
  Rcpp::NumericVector yAS_sumX(wrap(yAS_sum));
  Rcpp::NumericVector AS_sqX(wrap(AS_sq));
  return Rcpp::List::create(Rcpp::Named("A_sum") = A_sumX,
                            Rcpp::Named("AAt_sum") = AAt_sumX,
                            Rcpp::Named("yAS_sum") = yAS_sumX,
                            Rcpp::Named("AS_sq") = AS_sqX);
}
