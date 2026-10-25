#include <RcppArmadillo.h>
#include <limits>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List zadr_irls_cpp(
    const arma::mat& y,
    const arma::mat& x,
    const Rcpp::List& Pidx,
    const arma::mat& logy,
    double tol = 1e-6,
    int maxit = 100
) {

  const int n = y.n_rows;
  const int D = y.n_cols;
  const int p = D - 1;
  const int K = x.n_cols;
  const int dim = K * p;

  // Initial values: exactly as in the original R function
  arma::mat beta(K, p, fill::zeros);

  double phi = 1.0;

  double loglik_old =
    -std::numeric_limits<double>::infinity();

    double loglik = NA_REAL;
    int iters = 0;

    for (int iter = 1; iter <= maxit; ++iter) {

      iters = iter;

      // --------------------------------------------------
      // 1. Linear predictors and probabilities
      // --------------------------------------------------

      arma::mat eta = x * beta;
      arma::mat exp_eta = arma::exp(eta);

      arma::vec sum_exp = arma::sum(exp_eta, 1);

      arma::mat mu(n, D);

      mu.col(0) = 1.0 / (1.0 + sum_exp);

      for (int j = 0; j < p; ++j) {
        mu.col(j + 1) =
          exp_eta.col(j) / (1.0 + sum_exp);
      }

      const double psi1_phi = R::trigamma(phi);

      // --------------------------------------------------
      // 2. Initialize likelihood and derivatives
      // --------------------------------------------------

      loglik = 0.0;

      arma::vec S_vec(dim, fill::zeros);
      arma::mat I_mat(dim, dim, fill::zeros);

      double S_phi = 0.0;
      double H_phi = 0.0;

      // --------------------------------------------------
      // 3. Observation loop
      // --------------------------------------------------

      for (int i = 0; i < n; ++i) {

        // Positive components; R indices are 1-based
        IntegerVector Pi_R = Pidx[i];

        const int mPi = Pi_R.size();

        arma::uvec Pi(mPi);

        for (int r = 0; r < mPi; ++r) {
          Pi[r] = Pi_R[r] - 1;
        }

        const arma::rowvec x_i = x.row(i);
        const arma::rowvec mu_i = mu.row(i);

        // ----------------------------------------------
        // Conditional probabilities
        // ----------------------------------------------

        double Mi = 0.0;

        for (int r = 0; r < mPi; ++r) {
          Mi += mu_i[Pi[r]];
        }

        arma::vec mu2(mPi);
        arma::vec alpha2(mPi);
        arma::vec y_log(mPi);

        for (int r = 0; r < mPi; ++r) {

          mu2[r] = mu_i[Pi[r]] / Mi;
          alpha2[r] = phi * mu2[r];
          y_log[r] = logy(i, Pi[r]);
        }

        // ----------------------------------------------
        // Log-likelihood
        // ----------------------------------------------

        double sum_lgamma = 0.0;
        double sum_alpha_y = 0.0;

        for (int r = 0; r < mPi; ++r) {

          sum_lgamma += R::lgammafn(alpha2[r]);

          sum_alpha_y +=
            (alpha2[r] - 1.0) * y_log[r];
        }

        loglik +=
          R::lgammafn(phi) - sum_lgamma + sum_alpha_y;

        // ----------------------------------------------
        // Derivatives with respect to phi
        // ----------------------------------------------

        double S_phi_i = R::digamma(phi);
        double H_phi_i = psi1_phi;

        for (int r = 0; r < mPi; ++r) {

          const double mr = mu2[r];

          S_phi_i -=
            mr * R::digamma(alpha2[r]);

          S_phi_i += mr * y_log[r];

          H_phi_i -=
            mr * mr * R::trigamma(alpha2[r]);
        }

        S_phi += S_phi_i;
        H_phi += H_phi_i;

        // ----------------------------------------------
        // J: derivative of mu with respect to eta
        // ----------------------------------------------

        arma::mat J(D, p, fill::zeros);

        for (int d = 0; d < D; ++d) {

          for (int k = 0; k < p; ++k) {

            const double indicator =
              (d == k + 1) ? 1.0 : 0.0;

            J(d, k) =
              mu_i[d] * (indicator - mu_i[k + 1]);
          }
        }

        // ----------------------------------------------
        // Jp and J2
        // ----------------------------------------------

        arma::mat Jp(mPi, p);

        for (int r = 0; r < mPi; ++r) {
          Jp.row(r) = J.row(Pi[r]);
        }

        const arma::rowvec Jsum = arma::sum(Jp, 0);

        arma::mat J2(mPi, p);

        for (int r = 0; r < mPi; ++r) {

          for (int k = 0; k < p; ++k) {

            J2(r, k) =
              (Jp(r, k) - mu2[r] * Jsum[k]) / Mi;
          }
        }

        // ----------------------------------------------
        // h_i and C
        // ----------------------------------------------

        arma::vec h_i(mPi);
        arma::vec psi1_alpha2(mPi);

        for (int r = 0; r < mPi; ++r) {

          h_i[r] =
            y_log[r] - R::digamma(alpha2[r]);

          psi1_alpha2[r] =
            R::trigamma(alpha2[r]);
        }

        // Preserve the original R construction:
        // C = diag(trigamma(alpha2)) - trigamma(phi)

        arma::mat C(mPi, mPi);
        C.fill(-psi1_phi);

        C.diag() += psi1_alpha2;

        // Preserve the original multiplication order
        arma::mat JtCJ = (J2.t() * C) * J2;

        arma::vec Jt_h = J2.t() * h_i;

        // ----------------------------------------------
        // Accumulate score
        // ----------------------------------------------

        for (int j = 0; j < p; ++j) {

          for (int k = 0; k < K; ++k) {

            S_vec[k + K * j] +=
              phi * x_i[k] * Jt_h[j];
          }
        }

        // ----------------------------------------------
        // Accumulate information matrix
        //
        // Equivalent to:
        // I_mat <- I_mat +
        //   phi^2 * kronecker(JtCJ, tcrossprod(x_i))
        //
        // Calculate each block directly, avoiding
        // explicit Kronecker products and txi storage.
        // ----------------------------------------------

        const double phi2 = phi * phi;

        for (int a = 0; a < p; ++a) {

          for (int b = 0; b < p; ++b) {

            const double coeff = phi2 * JtCJ(a, b);

            for (int r = 0; r < K; ++r) {

              const int row = r + K * a;

              for (int c = 0; c < K; ++c) {

                const int col = c + K * b;

                I_mat(row, col) +=
                  coeff * x_i[r] * x_i[c];
              }
            }
          }
        }

      } // end observation loop

      // --------------------------------------------------
      // 4. Convergence test
      //
      // IMPORTANT: test before updating beta and phi,
      // just as in the original R implementation.
      // --------------------------------------------------

      if (std::abs(loglik - loglik_old) < tol) {
        break;
      }

      // --------------------------------------------------
      // 5. Newton update for beta
      // --------------------------------------------------

      arma::vec beta_vec = arma::vectorise(beta);

      arma::vec delta = arma::solve(I_mat, S_vec);

      beta_vec += delta;

      beta = arma::reshape(beta_vec, K, p);

      // --------------------------------------------------
      // 6. Newton update for phi
      // --------------------------------------------------

      double phi_new = phi - S_phi / H_phi;

      if (phi_new <= 0.0) {
        phi_new = 1e-4;
      }

      phi = phi_new;

      // Same convergence bookkeeping as the original R
      loglik_old = loglik;
    }

    return Rcpp::List::create(
      Rcpp::Named("be") = beta,
      Rcpp::Named("phi") = phi,
      Rcpp::Named("loglik") = loglik,
      Rcpp::Named("iters") = iters
    );
}
