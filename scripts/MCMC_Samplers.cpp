// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double lkernel_beta_i_cpp(const arma::vec& beta2_i, 
                          const arma::vec& mu_beta, 
                          const arma::mat& iSigma_beta) {
  arma::vec diff = beta2_i - mu_beta;
  return -0.5 * arma::as_scalar(diff.t() * iSigma_beta * diff);
}

// [[Rcpp::export]]
double llong_i_cpp(const arma::vec& y_i,
                   const arma::vec& t_i,
                   const arma::vec& beta_i,
                   double sigma2,
                   const arma::vec& kappa_i) {
  arma::vec psi_i = (beta_i(0) + kappa_i * beta_i(2)) + beta_i(1) * t_i;

  arma::vec resid = y_i - psi_i;
  
  // Calculate log-likelihood
  double n = static_cast<double>(y_i.n_elem);
  double ll = -0.5 * n * std::log(2.0 * M_PI * sigma2) - 
    arma::dot(resid, resid) / (2.0 * sigma2);
  
  return ll;
}

// [[Rcpp::export]]
arma::vec log_min_cpp(const arma::vec& x) {
  arma::vec result = arma::log(arma::clamp(x, DBL_MIN, INFINITY));
  return result;
}

// [[Rcpp::export]]
double lsurv_i_cpp(double T_i,
                   int delta_i,
                   const arma::rowvec& X_i,
                   const arma::vec& alpha,
                   const arma::vec& beta_i,
                   double gamma,
                   const arma::vec& lambda,
                   const arma::vec& u,
                   const arma::field<arma::mat>& int_bounds_i,
                   int k_T_i) {
  
  int K = lambda.n_elem;
  double X_i_alpha = arma::as_scalar(X_i * alpha);
  double int_sum = 0.0;
  
  // Loop over the piecewise intervals
  for (int k = 0; k < K; k++) {
    if (u(k) >= T_i) break;
    
    double lambda_k = std::max(lambda(k), DBL_MIN);
    
    // Extract the k-th block matrix (may be 0x2 if empty)
    arma::mat blocks = int_bounds_i(k);
    
    if (blocks.n_rows == 0) continue;
    
    double b0 = beta_i(0);
    double b1 = beta_i(1);
    
    if (std::abs(gamma * b1) < 1e-6) {
      arma::vec UB = blocks.col(1);  
      arma::vec LB = blocks.col(0); 
      double sum_diff = arma::accu(UB - LB);
      int_sum += lambda_k * std::exp(gamma * b0) * sum_diff;
    } else {
      arma::vec LB = blocks.col(0);
      arma::vec UB = blocks.col(1);
      
      arma::vec piece = arma::exp(gamma * b1 * LB) % 
                        arma::expm1(gamma * b1 * (UB - LB)) / (gamma * b1);
      
      int_sum += lambda_k * std::exp(gamma * b0) * arma::accu(piece);
    }
  }
  
  // Note: k_T_i is 1-based in R, so subtract 1 for 0-based C++ indexing
  arma::vec lambda_vec = {lambda(k_T_i - 1)};
  double log_lambda = log_min_cpp(lambda_vec)(0);
  
  double log_density = delta_i * (log_lambda + 
                                 (gamma * (beta_i(0) + beta_i(1) * T_i) + 
                                 X_i_alpha)) - (std::exp(X_i_alpha) * int_sum);
  
  return log_density;
}

//Function for sampling from a standard normal distribution
// [[Rcpp::export]]
arma::vec rnormSNRcpp(int n) {
  arma::vec out(n);
  for (int i = 0; i < n; i++) out(i) = R::rnorm(0, 1);
  return out;
}

//Function for sampling from a multivariate normal distribution
// [[Rcpp::export]]
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma) {
  int ncols = sigma.n_cols;
  arma::vec Yvec(n * ncols);
  Yvec = rnormSNRcpp(n * ncols);
  arma::mat Y = arma::reshape(Yvec, n, ncols);
  return arma::trans(arma::repmat(mean, 1, n).t() + Y * arma::chol(sigma));
}

// [[Rcpp::export]]
arma::vec sample_beta_i_cpp(const arma::vec& beta_i_curr,
                        const arma::vec& y_i,
                        const arma::vec& t_i,
                        double T_i,
                        int delta_i,
                        const arma::rowvec& X_i,
                        const arma::vec& alpha,
                        double gamma,
                        const arma::vec& lambda,
                        const arma::vec& u,
                        const arma::vec& mu_beta,
                        const arma::vec& kappa_i,
                        const arma::field<arma::mat>& int_bounds_i,
                        int k_T_i,
                        const arma::mat& Sigma_beta,
                        double sigma2,
                        const arma::mat& Sigma_prop) {
  
  // Propose new beta_i from MVN centered at current beta_i
  arma::vec beta2_curr = beta_i_curr.subvec(0, 1);  // Restrict to first 2 elements
  arma::mat beta2_prop_mat = rmvnormRcpp(1, beta2_curr, Sigma_prop);
  
  // Debug Code
  // Rcpp::Rcout << "beta2_prop_mat dimensions: " << beta2_prop_mat.n_rows << " x " << beta2_prop_mat.n_cols << "\n";
  // Rcpp::Rcout << "beta2_prop_mat:\n" << beta2_prop_mat << "\n";
  
  arma::vec beta2_prop = beta2_prop_mat.col(0);
  
  // Debug Code
  // Rcpp::Rcout << "beta2_prop dimensions: " << beta2_prop.n_elem << "\n";
  // Rcpp::Rcout << "beta2_prop:\n" << beta2_prop << "\n";
  
  arma::vec beta_i_prop(3);
  beta_i_prop(0) = beta2_prop(0);
  beta_i_prop(1) = beta2_prop(1);
  beta_i_prop(2) = beta_i_curr(2);  // Keep spike element unchanged
  
  // Compute log-posterior at current beta_i
  arma::mat iSigma_beta = arma::inv(Sigma_beta);
  double lprior_curr = lkernel_beta_i_cpp(beta2_curr, mu_beta, iSigma_beta);
  double llong_curr = llong_i_cpp(y_i, t_i, beta_i_curr, sigma2, kappa_i);
  double lsurv_curr = lsurv_i_cpp(T_i, delta_i, X_i, alpha, beta_i_curr, gamma, 
                              lambda, u, int_bounds_i, k_T_i);
  double lpost_curr = lprior_curr + llong_curr + lsurv_curr;
  
  // Compute log-posterior at proposal beta_i*
  double lprior_prop = lkernel_beta_i_cpp(beta2_prop, mu_beta, iSigma_beta);
  double llong_prop = llong_i_cpp(y_i, t_i, beta_i_prop, sigma2, kappa_i);
  double lsurv_prop = lsurv_i_cpp(T_i, delta_i, X_i, alpha, beta_i_prop, gamma,
                              lambda, u, int_bounds_i, k_T_i);
  double lpost_prop = lprior_prop + llong_prop + lsurv_prop;
  
  // Accept/reject
  double log_ratio = lpost_prop - lpost_curr;
  double log_u = std::log(arma::randu());  // log(uniform(0,1))
  
  if (log_u < log_ratio) {
    return beta_i_prop;
  } else { 
    return beta_i_curr;
  }
}

// [[Rcpp::export]]
arma::vec sample_beta_spike_i_cpp(const arma::vec& y_i,
                                  const arma::vec& t_i,
                                  const arma::uvec& active_idx_i,
                                  const arma::vec& beta_i,
                                  double sigma_sq,
                                  double mu_spike,
                                  double nu2_spike,
                                  double max_t = 6.01) {
  
  arma::vec beta_i_new = beta_i;  // Create a copy to modify
  
  // Check if active_idx_i is empty
  if (active_idx_i.n_elem == 0) {
    beta_i_new(2) = 0.0;
    return beta_i_new;
  } 
  
  // Residuals excluding the spike term
  double b0 = beta_i(0);
  double b1 = beta_i(1);
  arma::vec r_i = y_i - (b0 + b1 * t_i);
  
  // Get number of active observations
  int n_act = active_idx_i.n_elem;
  
  // Posterior precision and variance
  double prec_post = (1.0 / nu2_spike) + n_act / sigma_sq;
  double var_post = 1.0 / prec_post;
  
  // Sum of residuals at active indices
  // Req. conversion from 1-based to 0-based indexing
  double sum_r = 0.0;
  for (arma::uword i = 0; i < active_idx_i.n_elem; i++) {
    sum_r += r_i(active_idx_i(i) - 1);  // Subtract 1 for 0-based indexing
  } 
  
  // Posterior mean
  double mean_post = var_post * ((mu_spike / nu2_spike) + (sum_r / sigma_sq));
  
  // Sample new spike coefficient
  beta_i_new(2) = mean_post + std::sqrt(var_post) * R::rnorm(0, 1);
  
  return beta_i_new;
} 

// [[Rcpp::export]]
double sample_mu_beta0_cpp(const arma::mat& beta,
                           const arma::mat& Sigma_beta,
                           double mu_beta1,
                           double tau0_sq,
                           double b0) {
  
  arma::vec beta0 = beta.col(0);
  arma::vec beta1 = beta.col(1);
  int N = beta.n_rows;
  
  // Extract elements from Sigma_beta
  double sigma_beta0_sq = Sigma_beta(0, 0);
  double sigma_beta1_sq = Sigma_beta(1, 1);
  double sigma_beta0 = std::sqrt(sigma_beta0_sq);
  double sigma_beta1 = std::sqrt(sigma_beta1_sq);
  double rho = Sigma_beta(0, 1) / (sigma_beta0 * sigma_beta1);
  
  // Compute posterior parameters
  double a = (N / (sigma_beta0_sq * (1.0 - rho * rho))) + (1.0 / tau0_sq);
  
  arma::vec beta_diff = beta1 - mu_beta1;
  arma::vec adjusted = beta0 - ((rho * sigma_beta0 / sigma_beta1) * beta_diff);
  double sum_adjusted = arma::accu(adjusted);
  
  double b = (1.0 / (sigma_beta0_sq * (1.0 - rho * rho))) * sum_adjusted + (b0 / tau0_sq);
  
  // Sample from posterior
  double mu_beta0 = (b / a) + std::sqrt(1.0 / a) * R::rnorm(0, 1);
  
  return mu_beta0;
}

// [[Rcpp::export]]
double sample_mu_beta1_cpp(const arma::mat& beta,
                           const arma::mat& Sigma_beta,
                           double mu_beta0,
                           double tau1_sq,
                           double b1) {
  
  arma::vec beta0 = beta.col(0);  
  arma::vec beta1 = beta.col(1);
  int N = beta.n_rows;
  
  // Extract elements from Sigma_beta
  double sigma_beta0_sq = Sigma_beta(0, 0);
  double sigma_beta1_sq = Sigma_beta(1, 1);
  double sigma_beta0 = std::sqrt(sigma_beta0_sq);
  double sigma_beta1 = std::sqrt(sigma_beta1_sq);
  double rho = Sigma_beta(0, 1) / (sigma_beta0 * sigma_beta1);
  
  // Compute posterior parameters
  double a = (N / (sigma_beta1_sq * (1.0 - rho * rho))) + (1.0 / tau1_sq);
  
  arma::vec beta_diff = beta0 - mu_beta0;
  arma::vec adjusted = beta1 - ((rho * sigma_beta1 / sigma_beta0) * beta_diff);
  double sum_adjusted = arma::accu(adjusted);
  
  double b = (1.0 / (sigma_beta1_sq * (1.0 - rho * rho))) * sum_adjusted + (b1 / tau1_sq);
  
  // Sample from posterior
  double mu_beta1 = (b / a) + std::sqrt(1.0 / a) * R::rnorm(0, 1);
  
  return mu_beta1;
}

// [[Rcpp::export]]
double rchisqRcpp(double df) {
  return R::rchisq(df);
}

// Sample from a Wishart distribution using the Bartlett decomposition
// [[Rcpp::export]]
arma::mat rwishRcpp(double nu0, arma::mat const& S0) {
  int p = S0.n_rows;
  arma::mat L = arma::chol(S0);
  arma::mat A(p, p, arma::fill::zeros);
  for (int i = 0; i < p; i++) A(i, i) = sqrt(rchisqRcpp(nu0 - i));
  if (p > 1) {
    arma::vec RandSN = rnormSNRcpp(p * (p - 1) / 2);
    int counter = 0;
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < i; j++){
        A(j, i) = RandSN(counter);
        counter++;
      }
    }
  }
  arma::mat AL = A * L;
  return arma::trans(AL) * AL;
}

// [[Rcpp::export]]
arma::mat sample_Sigma_beta_cpp(double nu0,
                                const arma::mat& S0,
                                const arma::mat& beta,
                                const arma::vec& mu_beta) {
  int N = beta.n_rows;

  arma::mat Smu(2, 2, arma::fill::zeros);

  for (int i = 0; i < N; i++) {
    arma::vec beta_i = arma::trans(beta.row(i).cols(0, 1));  // [beta0_i, beta1_i]
    arma::vec diff = beta_i - mu_beta;
    Smu += diff * diff.t();
  }

  return arma::inv(rwishRcpp(nu0 + N, arma::inv(S0 + Smu)));
}

// [[Rcpp::export]]
double sample_sigma_sq_cpp(double a, double b,
                           const arma::mat& beta,
                           const arma::field<arma::vec>& Y,
                           const arma::field<arma::vec>& t,
                           const arma::field<arma::vec>& kappa,
                           double max_t = 6.01) {

  int N = Y.n_elem;

  double total_n = 0.0;
  for (int i = 0; i < N; i++) total_n += static_cast<double>(Y(i).n_elem);

  double shape = a + (0.5 * total_n);

  double sum_resids = 0.0;

  for (int i = 0; i < N; i++) {
    const arma::vec& y_i = Y(i);
    const arma::vec& t_i = t(i);
    const arma::vec& kappa_i = kappa(i);

    double b0 = beta(i, 0);
    double b1 = beta(i, 1);
    double b2 = beta(i, 2);

    arma::vec psi_i = (b0 + kappa_i * b2) + b1 * t_i;
    arma::vec resid = y_i - psi_i;

    sum_resids += arma::dot(resid, resid);
  }

  double rate = (.5 * sum_resids) + b;

  // R's rgamma(n, shape, rate) uses rate parameterization; R::rgamma() takes (shape, scale)
  double sigma_sq = 1.0 / R::rgamma(shape, 1.0 / rate);

  return sigma_sq;
}

// [[Rcpp::export]]
double sample_lambda_k_cpp(int k,
                           const arma::vec& u,
                           const arma::vec& T,
                           const arma::ivec& delta,
                           const arma::mat& X,
                           const arma::vec& alpha,
                           const arma::mat& beta,
                           double gamma,
                           const Rcpp::List& int_bounds,
                           double a_k,
                           double b_k,
                           double max_t = 6.01) {

  // k is 1-based (as in R)
  double u_lower = u(k - 1);
  double u_upper = u(k);

  int N = X.n_rows;

  int n_k = 0;
  for (int i = 0; i < N; i++) {
    if (delta(i) == 1 && T(i) >= u_lower && T(i) < u_upper) {
      n_k++;
    }
  }

  double total_integral = 0.0;

  for (int i = 0; i < N; i++) {
    double T_i = T(i);
    if (T_i < u_lower) continue;

    arma::rowvec X_i = X.row(i);
    double X_i_alpha = arma::as_scalar(X_i * alpha);

    double b0 = beta(i, 0);
    double b1 = beta(i, 1);

    Rcpp::List patient_bounds = int_bounds[i];
    arma::mat blocks = Rcpp::as<arma::mat>(patient_bounds[k - 1]);

    if (blocks.n_rows == 0) continue;

    arma::vec LB = blocks.col(0);
    arma::vec UB = blocks.col(1);

    if (std::abs(gamma * b1) < 1e-6) {
      total_integral += std::exp(X_i_alpha + gamma * b0) * arma::accu(UB - LB);
    } else {
      arma::vec piece = (arma::exp(gamma * b1 * UB) - arma::exp(gamma * b1 * LB)) / (gamma * b1);
      total_integral += std::exp(X_i_alpha + gamma * b0) * arma::accu(piece);
    }
  }

  double shape_post = a_k + n_k;
  double rate_post = std::max(b_k + total_integral, DBL_MIN);

  double lambda_k = R::rgamma(shape_post, 1.0 / rate_post);

  return lambda_k;
}

// [[Rcpp::export]]
double lsurv_all_cpp(const arma::vec& T,
                     const arma::ivec& delta,
                     const arma::mat& X,
                     const arma::vec& alpha,
                     const arma::mat& beta,
                     double gamma,
                     const arma::vec& lambda,
                     const arma::vec& u,
                     const Rcpp::List& int_bounds,
                     const arma::ivec& k_T) {

  int N = T.n_elem;
  double total_log_density = 0.0;

  for (int i = 0; i < N; i++) {
    arma::field<arma::mat> int_bounds_i = Rcpp::as<arma::field<arma::mat>>(int_bounds[i]);

    double li = lsurv_i_cpp(
      T(i),
      delta(i),
      X.row(i),
      alpha,
      arma::vec(beta.row(i).t()),
      gamma,
      lambda,
      u,
      int_bounds_i,
      k_T(i)
    );

    if (!std::isfinite(li)) {
      Rcpp::Rcout << "Non-finite ll at i = " << (i + 1) << "\n";
    }

    total_log_density += li;
  }

  return total_log_density;
}

// [[Rcpp::export]]
arma::vec sample_alpha_cpp(arma::vec alpha_curr,
                           const arma::mat& beta,
                           const arma::vec& T,
                           const arma::ivec& delta,
                           const arma::mat& X,
                           double gamma,
                           const arma::vec& lambda,
                           const arma::vec& u,
                           const Rcpp::List& int_bounds,
                           const arma::ivec& k_T,
                           const arma::vec& mu_alpha,
                           const arma::mat& Sigma_alpha,
                           const arma::mat& Sigma_prop_alpha) {

  int p = alpha_curr.n_elem;
  arma::vec alpha_prop = alpha_curr;

  double lsurv_curr = lsurv_all_cpp(T, delta, X, alpha_curr, beta, gamma, lambda, u, int_bounds, k_T);

  for (int j = 0; j < p; j++) {
    alpha_prop(j) = R::rnorm(alpha_curr(j), Sigma_prop_alpha(j, j));

    double lprior_curr = R::dnorm(alpha_curr(j), mu_alpha(j), Sigma_alpha(j, j), true);
    double lpost_curr = lsurv_curr + lprior_curr;

    double lsurv_prop = lsurv_all_cpp(T, delta, X, alpha_prop, beta, gamma, lambda, u, int_bounds, k_T);
    double lprior_prop = R::dnorm(alpha_prop(j), mu_alpha(j), Sigma_alpha(j, j), true);
    double lpost_prop = lsurv_prop + lprior_prop;

    double log_ratio = lpost_prop - lpost_curr;
    if (std::log(R::runif(0, 1)) < log_ratio) {
      alpha_curr(j) = alpha_prop(j);
      lsurv_curr = lsurv_prop;  // propagate updated log-surv density to next component
    } else {
      alpha_prop(j) = alpha_curr(j);
    }
  }

  return alpha_curr;
}

// [[Rcpp::export]]
double sample_gamma_cpp(double gamma_curr,
                        const arma::mat& beta,
                        const arma::vec& T,
                        const arma::ivec& delta,
                        const arma::mat& X,
                        const arma::vec& alpha,
                        const arma::vec& lambda,
                        const arma::vec& u,
                        const Rcpp::List& int_bounds,
                        const arma::ivec& k_T,
                        double mu_gamma,
                        double sigma_sq_gamma,
                        double sigma_sq_prop) {

  double gamma_prop = R::rnorm(gamma_curr, std::sqrt(sigma_sq_prop));

  double lsurv_curr = lsurv_all_cpp(T, delta, X, alpha, beta, gamma_curr, lambda, u, int_bounds, k_T);
  double lprior_curr = R::dnorm(gamma_curr, mu_gamma, std::sqrt(sigma_sq_gamma), true);
  double lpost_curr = lsurv_curr + lprior_curr;

  double lsurv_prop = lsurv_all_cpp(T, delta, X, alpha, beta, gamma_prop, lambda, u, int_bounds, k_T);
  double lprior_prop = R::dnorm(gamma_prop, mu_gamma, std::sqrt(sigma_sq_gamma), true);
  double lpost_prop = lsurv_prop + lprior_prop;

  double log_ratio = lpost_prop - lpost_curr;

  if (std::log(R::runif(0, 1)) < log_ratio) {
    return gamma_prop;
  } else {
    return gamma_curr;
  }
}

// [[Rcpp::export]]
double pilot_adapt_cpp(double tuning_param, double accept_rate) {

  // adjust tuning parameter by scaling existing parameter based on current acceptance rate
  if (accept_rate >= 0.90) {
    tuning_param = tuning_param * 1.3;
  } else if ((accept_rate >= 0.75) && (accept_rate < 0.90)) {
    tuning_param = tuning_param * 1.2;
  } else if ((accept_rate >= 0.45) && (accept_rate < 0.75)) {
    tuning_param = tuning_param * 1.1;
  } else if ((accept_rate <= 0.25) && (accept_rate > 0.15)) {
    tuning_param = tuning_param * 0.9;
  } else if ((accept_rate <= 0.15) && (accept_rate > 0.10)) {
    tuning_param = tuning_param * 0.8;
  } else if (accept_rate <= 0.10) {
    tuning_param = tuning_param * 0.7;
  }

  return tuning_param;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
beta2_i <- c(1.2, -0.3)
mu_beta <- c(0, 0)
iSigma_beta <- matrix(c(2, 0.5, 0.5, 1), 2, 2)

lkernel_beta_i_cpp(beta2_i, mu_beta, iSigma_beta)

y_i <- c(5.2, 4.8, 4.5, 4.1, 3.9)
t_i <- c(0, 1, 2, 3, 4)
beta_i <- c(4.0, -0.3, 1.0)
sigma2 <- 0.5
kappa_i <- c(0, 0, 1, 0, 0)

llong_i_cpp(y_i, t_i, beta_i, sigma2, kappa_i)

x <- c(1e-300, 1e-200, 1e-100, 1e-10, 1.0)
log_min_cpp(x)

T_i <- 5.0
delta_i <- 1
X_i <- matrix(c(1, 0.5, -0.3), nrow = 1)  # 1x3 design matrix
alpha <- c(0.2, -0.1, 0.15)
beta_i <- c(4.0, -0.3, 1.0)
gamma <- 0.8
lambda <- c(0.1, 0.15, 0.2, 0.25)
u <- c(0, 2, 4, 6)
k_T_i <- 3

int_bounds_i <- list(
  matrix(c(0, 1, 1, 2), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
  matrix(c(2, 3, 3, 4), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
  matrix(c(4, 5), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
  matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("LB", "UB")))
)

lsurv_i_cpp(T_i, delta_i, X_i, alpha, beta_i, gamma, lambda, u, int_bounds_i, k_T_i)

beta_i_curr <- c(3.0, -0.5, -1.2)
y_i <- c(5.2, 4.8, 4.5, 4.1, 3.9)
t_i <- c(0, 1, 2, 3, 4)
T_i <- 5.0
delta_i <- 1L
X_i <- matrix(c(1, 0.5, -0.3), nrow = 1)
alpha <- c(0.2, -0.1, 0.15)
gamma <- 0.8
lambda <- c(0.1, 0.15, 0.2, 0.25)
u <- c(0, 2, 4, 6)
mu_beta <- c(3.0, -0.5)
kappa_i <- c(0, 0, 1, 1, 1)
int_bounds_i <- list(
  matrix(c(0, 1, 1, 2), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
  matrix(c(2, 3, 3, 4), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
  matrix(c(4, 5), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
  matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("LB", "UB")))
)
k_T_i <- 3L
Sigma_beta <- matrix(c(1.0, 0.3, 0.3, 0.8), nrow = 2)
sigma2 <- 0.5
Sigma_prop <- diag(c(0.1, 0.1))

sample_beta_i_cpp(beta_i_curr, y_i, t_i, T_i, delta_i, X_i, alpha, gamma,
                  lambda, u, mu_beta, kappa_i, int_bounds_i, k_T_i, Sigma_beta, sigma2, Sigma_prop)

y_i <- c(5.2, 4.8, 4.5, 3.2, 2.9)
t_i <- c(0, 1, 2, 3, 4)
active_idx_i <- c(3, 4, 5)  # 1-based indices in R (observations 3, 4, 5 are active)
beta_i <- c(5.0, -0.5, 1.0)  # [intercept, slope, spike]
sigma_sq <- 0.5
mu_spike <- -1.2
nu2_spike <- 1.0

sample_beta_spike_i_cpp(y_i, t_i, active_idx_i, beta_i, sigma_sq, mu_spike, nu2_spike)

nu0 <- 5
S0 <- matrix(c(1.0, 0.3, 0.3, 0.8), nrow = 2)

rwishRcpp(nu0, S0)

beta <- matrix(c(4.0, 3.8, 4.2, -0.3, -0.5, -0.2, 1.0, 0.5, -1.0), nrow = 3, ncol = 3)
mu_beta <- c(4.0, -0.3)

sample_Sigma_beta_cpp(nu0, S0, beta, mu_beta)

Y_field <- list(c(5.2, 4.8, 4.5), c(6.0, 5.5))
t_field <- list(c(0, 1, 2), c(0, 1))
kappa_field <- list(c(0, 0, 1), c(0, 1))
beta_sig <- matrix(c(4.0, 5.0, -0.3, -0.5, 1.0, 0.5), nrow = 2, ncol = 3)
a <- 0.001
b <- 0.001

sample_sigma_sq_cpp(a, b, beta_sig, Y_field, t_field, kappa_field)

u <- c(0, 2, 4, 6)
T_vec <- c(5.0, 3.0)
delta_vec <- c(1L, 0L)
X_mat <- matrix(c(1, 1, 0.5, -0.2, -0.3, 0.1), nrow = 2, ncol = 3)
alpha_vec <- c(0.2, -0.1, 0.15)
beta_lam <- matrix(c(4.0, 3.5, -0.3, -0.4, 1.0, 0.0), nrow = 2, ncol = 3)
gamma_val <- 0.8
a_k <- 0.001
b_k <- 0.001

int_bounds_list <- list(
  list(
    matrix(c(0, 1, 1, 2), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
    matrix(c(2, 3, 3, 4), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
    matrix(c(4, 5), ncol = 2, dimnames = list(NULL, c("LB", "UB")))
  ),
  list(
    matrix(c(0, 1), ncol = 2, dimnames = list(NULL, c("LB", "UB"))),
    matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("LB", "UB"))),
    matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("LB", "UB")))
  )
)

sample_lambda_k_cpp(1, u, T_vec, delta_vec, X_mat, alpha_vec, beta_lam, gamma_val,
                    int_bounds_list, a_k, b_k)

k_T_vec <- c(3L, 2L)

lsurv_all_cpp(T_vec, delta_vec, X_mat, alpha_vec, beta_lam, gamma_val, lambda,
             u, int_bounds_list, k_T_vec)

mu_alpha_vec <- c(0, 0, 0)
Sigma_alpha_mat <- diag(100, 3)
Sigma_prop_alpha_mat <- diag(0.1, 3)

sample_alpha_cpp(alpha_vec, beta_lam, T_vec, delta_vec, X_mat, gamma_val, lambda,
                 u, int_bounds_list, k_T_vec, mu_alpha_vec, Sigma_alpha_mat, Sigma_prop_alpha_mat)

mu_gamma_val <- 0
sigma_sq_gamma_val <- 1
sigma_sq_prop_val <- 0.001

sample_gamma_cpp(gamma_val, beta_lam, T_vec, delta_vec, X_mat, alpha_vec, lambda,
                 u, int_bounds_list, k_T_vec, mu_gamma_val, sigma_sq_gamma_val, sigma_sq_prop_val)

pilot_adapt_cpp(0.1, 0.95)
pilot_adapt_cpp(0.1, 0.80)
pilot_adapt_cpp(0.1, 0.50)
pilot_adapt_cpp(0.1, 0.20)
pilot_adapt_cpp(0.1, 0.12)
pilot_adapt_cpp(0.1, 0.05)

*/
