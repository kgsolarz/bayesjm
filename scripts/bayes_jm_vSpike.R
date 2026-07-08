library(readr)
library(coda)
library(MASS)
library(dplyr)
library(forcats)
library(survival)
library(invgamma)
library(ggplot2)
library(tidyr)


pain_data <- read.csv("/hpc/group/berchucklab/kgs35/data/pain_030226.csv")
active_trt <- read_csv("/hpc/group/berchucklab/kgs35/data/surg_030226.csv")
elix <- read_csv(("/hpc/group/berchucklab/kgs35/data/elix_030226.csv"))


### Data pre-processing
pain_data <- pain_data |>
  left_join(elix, by = "obs")

pain_data <- pain_data |>
  mutate(stage = if_else(is.na(stage), "unknown", as.character(stage)),
         insurance = case_when(
           is.na(insurance) ~ "other",
           insurance == "managed care" ~ "medicaid",
           TRUE ~ insurance
         ),
         race = if_else(race == "white", race, "non-white")) |>
  mutate(stage = as.factor(stage),
         race = as.factor(race),
         ethnicity = as.factor(ethnicity),
         insurance = as.factor(insurance),
  ) |>
  mutate(ethnicity = fct_relevel(ethnicity, "not hispanic", "hispanic"),
         race = fct_relevel(race, "white", "non-white")) |>
  filter(sex == "female") |>
  filter(!is.na(adi)) |>
  filter(adi != "GQ", adi != "PH") |>
  mutate(adi = as.numeric(adi))

patients <- unique(pain_data$obs) ## Cohort of 2484 patients
N = length(patients)

### Construct covariate matrix for all patients (demos only)
demo_data <- pain_data |>
  dplyr::select(-c(pain:endocrine), -pain_time) |>
  distinct() |>
  arrange(obs)

demo_data <- demo_data |>
  dplyr::select(obs, stage, age, race, insurance, adi, elix_mort)

demo_data$age <- scale(demo_data$age)
demo_data$adi <- scale(demo_data$adi)

key <- demo_data |>
  distinct(obs) |>
  arrange(obs) |>    
  mutate(obs_new = row_number())

# Relabel demo_data
demo_data <- demo_data |>
  left_join(key, by = "obs") |>
  mutate(obs = obs_new) |>
  dplyr::select(-obs_new)

# Relabel pain_data to match the same subset + new IDs
pain_data <- pain_data |>
  inner_join(key, by = "obs") |>
  mutate(obs = obs_new) |>
  dplyr::select(-obs_new)

demo_data <- model.matrix(~ ., data = demo_data)

p <- dim(demo_data)[2] - 2 ## Where X is a Nxp matrix, and X_i is 1xp

## Compute baseline-haz related data:
event_data <- pain_data |>
  dplyr::select(obs, event, time_final) |>
  distinct() |>
  filter(event == 1)

time_data <- pain_data |>
  dplyr::select(obs, event, time_final) |>
  distinct() |>
  arrange(obs)

summary = as.numeric(summary(event_data$time_final))

## Discuss w/Sam: cutoff points are not even - potentially an issue given censoring distribution
u <- c(0, summary[2], summary[3], summary[4], max(time_data$time_final) + .01) # add padding to upper bound on time to avoid downstream indexing issues
max_t <- max(time_data$time_final + .01)

## Initialize alpha using Cox PH coeff.
demo_df <- data.frame(demo_data)
demo_df <- demo_df |>
  left_join(time_data, by = "obs")

baseline_mod <- coxph(Surv(time_final, event) ~ stage1 + stage2 + stage3 + stage4 + stageunknown
                      + insurancemedicaid + insurancemedicare + insuranceother
                      + age + elix_mort + adi + racenon.white,
                      data = demo_df[, c(ncol(demo_df)-1, ncol(demo_df), 3:(ncol(demo_df)-2))])

## Addtl. Cox Regressions
# mod_stage <- coxph(Surv(time_final, event) ~ stage1 + stage2 + stage3 + stage4 + stageunknown, data = demo_df[, c(15, 16, 3:14)])
# mod_insurance <- coxph(Surv(time_final, event) ~ insurancemanaged.care + insurancemedicare + insuranceother, data = demo_df[, c(15, 16, 3:14)])
# mod_race <- coxph(Surv(time_final, event) ~ raceblack + raceother, data = demo_df[, c(15, 16, 3:14)])
# mod_age <- coxph(Surv(time_final, event) ~ age, data = demo_df[, c(15, 16, 3:14)])


alpha <- matrix(baseline_mod$coefficients , nrow = p, ncol = 1)

## Finalize design matrix, X
X <- demo_data[ , 3:(ncol(demo_df)-2)]

## Initialize lambda vector using empirical hazard rates in each time interval defined by cutoff points above [left closed, right open)
lambda <- c()

for (i in 1: (length(u) - 1)) {
  events_k <- event_data |>
    filter(time_final >= u[i] & time_final < u[i + 1])

  events <- length(unique(events_k$obs))

  risk_set <- time_data |>
    filter(time_final >= u[i])

  n_risk <- length(unique(risk_set$obs))

  lambda_k <- events / (n_risk * (u[i+1] - u[i]))

  lambda <- c(lambda, lambda_k)
}

## Initialize beta vector at the patient-level by performing OLS reg on individual pain data
beta <- matrix(NA, nrow = N, ncol = 3)

## Gen. list of lists for indexing
pain_data <- pain_data |>
  group_by(obs) |>
  arrange(obs, pain_time)

Y <- split(pain_data$pain, pain_data$obs)
t <- split(pain_data$pain_time, pain_data$obs)

agg_resids <- 0 # Initialize

for (i in 1:N) {
  df <- data.frame(pain = Y[[i]], time = t[[i]])

  fit <- lm(pain ~ time, data = df)
  beta[i, ] <- c(coef(fit), 0)  # Store intercept and slope for each beta_i

  agg_resids <- c(agg_resids, residuals(fit))
}

sigma_sq <- var(agg_resids)

mu_beta0 <- mean(beta[, 1])
mu_beta1 <- mean(beta[, 2])
mu_beta <- matrix(c(mu_beta0, mu_beta1), nrow = 1, ncol = 2)

Sigma_beta11 <- var(beta[, 1])
Sigma_beta22 <- var(beta[, 2])
Sigma_beta = matrix(c(Sigma_beta11, 0,  0, Sigma_beta22), nrow = 2, ncol = 2)

# # Generate surg. timeline
# active_trt <- active_trt |>
#   distinct()
# 
# sum_surg <- active_trt |>
#   group_by(obs) |>
#   summarize(count = n())
# 
# active_trt <- active_trt |>
#   mutate(end_trt = surg_time + (14/365))
# 
# active_trt <- active_trt |>
#   pivot_longer(cols = -obs,
#                names_to = "trt_t",
#                values_to = "w_vec") |>
#   dplyr::select(obs, w_vec, trt_t)
# 
# w <- split(active_trt$w_vec, active_trt$obs)

# 14 days in "years" 
window_len <- 14/365

surg_clean <- active_trt |>
  inner_join(key, by = "obs") |>
  mutate(obs = obs_new) |>
  dplyr::select(-obs_new)

# Start with unique surgeries per patient
surg_clean <- surg_clean |>
  distinct(obs, surg_time) |>
  arrange(obs, surg_time)

merged_windows <- surg_clean |>
  group_by(obs) |>
  arrange(surg_time, .by_group = TRUE) |>
  mutate(
    lag_surg = lag(surg_time),
    gap = surg_time - lag_surg,
    new_cluster = if_else(is.na(gap) | gap > window_len, 1L, 0L),
    cluster_id = cumsum(new_cluster)
  ) |>
  group_by(obs, cluster_id) |>
  summarize(
    window_start = first(surg_time),
    window_end   = max(surg_time) + window_len, 
    .groups = "drop"
  )

clean_trt <- merged_windows |>
    dplyr::select(-cluster_id) |>
    pivot_longer(cols = -obs,
                 names_to = "trt_t",
                 values_to = "w_vec") |>
    dplyr::select(obs, w_vec, trt_t)

clean_trt <- clean_trt |>
  filter(w_vec <= 6.00)

w <- split(clean_trt$w_vec, clean_trt$obs)

cuts <- vector("list", N)
Z <- integer(N)

z <- vector("list", N)
kappa <- vector("list", N)
active_idx <- vector("list", N)

int_bounds <- vector("list", N)

K <- length(lambda)
k_T <- integer(N)

sort_w <- function(w_i) {
  if (is.null(w_i) || length(w_i) == 0 || all(is.na(w_i))) return(numeric(0))
  sort(unique(w_i[is.finite(w_i)]))
}

for (i in 1:N) {
  
  # Patient-level vectors
  w_i <- sort_w(w[[as.character(i)]])
  cuts_i <- c(0, w_i, max_t)
  Z_i <- length(cuts_i) - 1
  
  cuts[[i]] <- cuts_i
  Z[i] <- Z_i
  
  # Define z/kappa at "observed" times (times corr. to recorded pain scores)
  t_i <- t[[i]]
  z_i <- findInterval(t_i, cuts_i, left.open = FALSE, rightmost.closed = FALSE, all.inside = TRUE)
  kappa_i <- as.numeric(z_i %% 2 == 0)
  
  z[[i]] <- z_i
  kappa[[i]] <- kappa_i
  active_idx[[i]] <- which(kappa_i == 1)
  
  # survival integral bounds (LB,UB) by baseline interval k, excl. active trt. windows (kappa==1)
  T_i <- time_data$time_final[time_data$obs == i][1]
  k_T[i] <- findInterval(T_i, u, all.inside = TRUE)
  bki <- vector("list", K)
  
  for (k in 1:K) {
    if (u[k] >= T_i) {
      bki[[k]] <- matrix(numeric(0), ncol = 2)
      next
    }
    
    ub_k <- min(T_i, u[k + 1])
    
    blocks_k <- vector("list", Z_i)
    m <- 0L
    
    for (j in 1:Z_i) {
      kappa_seg <- as.numeric(j %% 2 == 0)
      if (kappa_seg == 1) next  # excl. active trt. windows
      
      LB <- max(u[k], cuts_i[j], 0)
      UB <- min(ub_k, cuts_i[j + 1])
      if (UB <= LB) next
      
      m <- m + 1L
      blocks_k[[m]] <- c(LB, UB)
    }
    
    if (m == 0L) {
      bki[[k]] <- matrix(numeric(0), ncol = 2)
    } else {
      bki[[k]] <- do.call(rbind, blocks_k[seq_len(m)])
      colnames(bki[[k]]) <- c("LB", "UB")
    }
  }
  
  int_bounds[[i]] <- bki
}

### END OF DATA PRE-PROCESSING: BEGIN FUNCTION DEFINITION

### Define helper functions:

## Log of the kernel of the distribution of beta_i, for sampling from the posterior of beta_i
lkernel_beta_i <- function(beta2_i, mu_beta, iSigma_beta) { ## Requires inversion of the covariance matrix prior to inclusion as a function argument
  mu_beta <- as.numeric(mu_beta)
  diff <- beta2_i - mu_beta
  diff <- matrix(diff, ncol = 1)
  return(-0.5 * crossprod(diff, iSigma_beta %*% diff))
}

## Log of the longitudinal component of the joint likelihood function; contribution from the i^th patient
llong_i <- function(y_i, t_i, beta_i, sigma2, kappa_i) {
  ## y_i is vector of obs. pain scores, 
  ## t_i is vector of corresponding times for the i^th patient, 
  ## kappa_i is the vector of indicators (indicating active trt. windows, kappa == 1)
  
  psi_i <- (beta_i[1] + kappa_i * beta_i[3]) + beta_i[2] * t_i 
  resid <- y_i - psi_i
  ll <- (-0.5 * length(y_i) * log(2 * pi * sigma2)) - (sum(resid^2) / (2 * sigma2))
  return(ll)
}

log_min <- function(x) log(pmax(x, .Machine$double.xmin))

## Log of the survival component of the joint likelihood function
lsurv_i <- function(T_i, delta_i, X_i, alpha, beta_i, gamma, lambda, u, int_bounds_i, k_T_i) {
  ## T_i = min((T_i)^*, C_i), delta_i = binary indicator for event, X_i = design matrix for i^th patient
  K <- length(lambda)
  
  X_i_alpha <- as.numeric(X_i %*% alpha)
  
  int_sum <- 0 # Initialize integral at 0
  
  # Loop over the piecewise intervals & update the computed integral
  for (k in 1:K) {
    if (u[k] >= T_i) break
    lambda_k <- max(lambda[k], .Machine$double.xmin)
    
    blocks <- int_bounds_i[[k]]   # matrix with cols LB,UB or 0x2
    if (nrow(blocks) == 0) next
    
    b0 <- beta_i[1]
    b1 <- beta_i[2]
    
    if (abs(gamma * b1) < 1e-6) {
      int_sum <- int_sum + lambda_k * exp(gamma * b0) * sum(blocks[, "UB"] - blocks[, "LB"])
    } else {
      LB <- blocks[, "LB"]
      UB <- blocks[, "UB"]
      piece <- exp(gamma * b1 * LB) * expm1(gamma * b1 * (UB - LB)) / (gamma * b1)
      int_sum <- int_sum + lambda_k * exp(gamma * b0) * sum(piece)
    }
  }
  
  # Compute the full log-density for the survival component
  log_density <- delta_i * (log_min(lambda[k_T_i]) + (gamma * (beta_i[1] + beta_i[2] * T_i) + X_i_alpha)) - 
    (exp(X_i_alpha) * int_sum)
  
  return(log_density)
}

## Function to sample beta_i
sample_beta_i <- function(beta_i_curr, y_i, t_i, T_i, delta_i, X_i, 
                          alpha, gamma, lambda, u, mu_beta, kappa_i, int_bounds_i, k_T_i,
                          Sigma_beta, sigma2, Sigma_prop = diag(c(0.1, 0.1))) {
  
  ## Propose new beta_i, beta_i*, from a MVN distribution centered at current beta_i
  beta2_curr <- beta_i_curr[1:2]
  beta2_prop <- mvrnorm(1, mu = beta2_curr, Sigma = Sigma_prop)
  beta_i_prop  <- c(beta2_prop, beta_i_curr[3])
  
  ## Compute log-posterior at current beta_i
  iSigma_beta <- solve(Sigma_beta)
  lprior_curr <- lkernel_beta_i(beta2_curr, mu_beta, iSigma_beta)
  llong_curr <- llong_i(y_i, t_i, beta_i_curr, sigma2, kappa_i)
  # print(paste("log long density at current for patient =", llong_curr))
  lsurv_curr <- lsurv_i(T_i, delta_i, X_i, alpha, beta_i_curr, gamma, lambda, u, int_bounds_i, k_T_i)
  # print(paste("log surv density at current for patient =", lsurv_curr))
  lpost_curr <- lprior_curr + llong_curr + lsurv_curr
  
  ## Compute log-posterior at proposal beta_i*
  lprior_prop <- lkernel_beta_i(beta2_prop, mu_beta, iSigma_beta)
  llong_prop <- llong_i(y_i, t_i, beta_i_prop, sigma2, kappa_i)
  lsurv_prop <- lsurv_i(T_i, delta_i, X_i, alpha, beta_i_prop, gamma, lambda, u, int_bounds_i, k_T_i)
  # print(lsurv_prop)
  lpost_prop <- lprior_prop + llong_prop + lsurv_prop
  
  # Accept/reject 
  log_ratio <- lpost_prop - lpost_curr
  # print(log_ratio)
  if(log(runif(1)) < log_ratio) {
    return(beta_i_prop)
  } else {
    return(beta_i_curr)
  }
}

## Function to sample "spike parameter":

sample_beta_spike_i <- function(y_i, t_i, active_idx_i, beta_i, sigma_sq,
                           mu_spike, nu2_spike, max_t = 6.01) {
    
    if(is.null(active_idx_i) | length(active_idx_i) == 0 | all(is.na(active_idx_i))) {
      beta_i[3] <- 0
      return(beta_i)
    }
    
    # Residuals excluding the spike term
    b0 <- beta_i[1]
    b1 <- beta_i[2]
    r_i <- y_i - (b0 + b1 * t_i)
    
    # Productize over observations with kappa == 1
    idx <- active_idx_i
    n_act <- length(idx)
    
    prec_post <- (1 / nu2_spike) + n_act / sigma_sq
    var_post <- 1 / prec_post
    
    sum_r <- sum(r_i[idx])
    mean_post <- var_post * ((mu_spike / nu2_spike) + (sum_r / sigma_sq))
    
    beta_i[3] <- rnorm(1, mean_post, sqrt(var_post))
    return(beta_i)
}

## Function to sample mu_beta0
sample_mu_beta0 <- function(beta, Sigma_beta, mu_beta1, tau0_sq, b0) {
  beta0 <- beta[, 1]  # Extract beta_0 vector
  beta1 <- beta[, 2]  # Extract beta_1 vector
  N <- nrow(beta)
  
  # Extract elements from Sigma
  sigma_beta0_sq <- Sigma_beta[1, 1]
  sigma_beta1_sq <- Sigma_beta[2, 2]
  sigma_beta0 <- sqrt(sigma_beta0_sq)
  sigma_beta1 <- sqrt(sigma_beta1_sq)
  rho <- Sigma_beta[1, 2] / (sigma_beta0 * sigma_beta1)
  
  a <- (N / (sigma_beta0_sq * (1 - rho^2))) + (1 / tau0_sq)
  b <- (1 / (sigma_beta0_sq * (1 - rho^2)))*(sum(beta0 - ((rho*sigma_beta0 / sigma_beta1)*(beta1 - mu_beta1)))) + (b0 / tau0_sq)
  
  mu_beta0 <- rnorm(1, mean = b / a, sd = sqrt(1 / a))
  return(mu_beta0)
}

## Function to sample mu_beta1
sample_mu_beta1 <- function(beta, Sigma_beta, mu_beta0, tau1_sq, b1) {
  beta0 <- beta[, 1]  # Extract beta_0i vector
  beta1 <- beta[, 2]  # Extract beta_1i vector
  N <- nrow(beta)
  
  # Extract elements from Sigma
  sigma_beta0_sq <- Sigma_beta[1, 1]
  sigma_beta1_sq <- Sigma_beta[2, 2]
  sigma_beta0 <- sqrt(sigma_beta0_sq)
  sigma_beta1 <- sqrt(sigma_beta1_sq)
  rho <- Sigma_beta[1, 2] / (sigma_beta0 * sigma_beta1)
  
  a <- (N / (sigma_beta1_sq * (1 - rho^2))) + (1 / tau1_sq)
  b <- (1 / (sigma_beta1_sq * (1 - rho^2)))*(sum(beta1 - ((rho*sigma_beta1 / sigma_beta0)*(beta0 - mu_beta0)))) + (b1 / tau1_sq)
  
  mu_beta1 <- rnorm(1, mean = b / a, sd = sqrt(1 / a))
  
  return(mu_beta1)
}

## Function to sample from the Wishart distribution:
rwish <- function(n, nu0, S0) {
  ## n = number of random samples to draw
  sS0 <- chol(S0) ## Cholesky decomposition of S0; requires S0 is a symm. PD sq. matrix
  S <- array(dim = c(dim(S0), n))
  for(i in 1:n) {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    print(Z)
    S[,,i] <- t(Z)%*%Z
  }
  S[,,1:n]
}

## Function to sample Sigma_beta:
sample_Sigma_beta <- function(nu0, S0, beta, mu_beta) {
  N <- nrow(beta)
  
  Smu <- 0 # Initialize
  
  for (i in 1:N) {
    beta_i <- beta[i, 1:2 , drop = FALSE]
    
    diff <- beta_i - mu_beta
    
    comp <- t(diff)%*%diff
    
    Smu <- Smu + comp
  }
  
  Sigma_beta <- solve(rwish(1, nu0 + N, 
                            solve(S0 + Smu))) ## Depends on dims of BETA / mu_beta; matrix op must match desired dims of Sigma_beta
  return(Sigma_beta)
}

## Function to sample sigma^2
sample_sigma_sq <- function(a, b, beta, Y, t, kappa, max_t = 6.01) {
  ## Note: Y and t should be objects of type "list", where both are lists of lists 
  ## (ex: Y is a list of length N, where the i^th element of Y is a list of length n_i)
  N <- length(Y)
  shape <- a + (0.5 * sum(sapply(Y, length)))
  
  sum_resids <- 0
  
  for (i in 1:N) {
    y_i <- Y[[i]]
    t_i <- t[[i]]
    kappa_i <- kappa[[i]]
    
    beta_i <- beta[i, ]
    psi_i <- (beta_i[1] + kappa_i * beta_i[3]) + beta_i[2] * t_i
    resid <- y_i - psi_i
    
    sum_element <- sum((resid)^2)
    sum_resids <- sum_resids + sum_element
  }
  
  rate <- (.5 * sum_resids) + b
  sigma_sq <- invgamma::rinvgamma(1, shape, rate)
  
  return(sigma_sq)
}

## Function to sample lambda_k:
sample_lambda_k <- function(k, u, time_data, X, alpha, beta, gamma, int_bounds, max_t = 6.01, a_k, b_k) {
  
  ## Bounds for the k^th interval of time
  u_lower <- u[k]
  u_upper <- u[k + 1]
  
  events_k <- time_data$obs[time_data$event == 1 & time_data$time_final >= u_lower & time_data$time_final < u_upper]
  n_k <- length(unique(events_k))
  
  ## Compute integral term for all patients in risk set at u_lower; patients who have exp. event / censoring prior to u_lower do not contribute (indicator func.)
  total_integral <- 0 # Initialize 
  for (i in 1:nrow(X)) {
    T_i <- time_data$time_final[time_data$obs == i][1]
    
    if (T_i < u_lower) next
    
    X_i <- matrix(X[i, ], nrow = 1)
    X_i_alpha <- as.numeric(X_i %*% alpha)
    
    b0 <- beta[i, 1]
    b1 <- beta[i, 2]
    
    blocks <- int_bounds[[i]][[k]]
    if (nrow(blocks) == 0) next
    
    if (abs(gamma * b1) < 1e-6) {
      total_integral <- total_integral +
        exp(X_i_alpha + gamma * b0) * sum(blocks[, "UB"] - blocks[, "LB"])
    } else {
      LB <- blocks[, "LB"]
      UB <- blocks[, "UB"]
      piece <- (exp(gamma * b1 * UB) - exp(gamma * b1 * LB)) / (gamma * b1)
      total_integral <- total_integral + exp(X_i_alpha + gamma * b0) * sum(piece)
    }

  }
  
  shape_post <- a_k + n_k
  rate_post <- max(b_k + total_integral, .Machine$double.xmin)
  
  lambda_k <- rgamma(1, shape = shape_post, rate = rate_post)
  return(lambda_k)
}

## Log of the survival component of the joint likelihood function for ALL individuals
lsurv_all <- function(T, delta, X, alpha, beta, gamma, lambda, u, int_bounds, k_T) {
  N <- length(T)
  total_log_density <- 0 # Initialize
  
  # density_vec <- numeric(N)
  
  for (i in 1:N) {
    li <- lsurv_i(
      T_i = T[i],
      delta_i = delta[i],
      X_i = X[i, , drop = FALSE],
      alpha = alpha,
      beta_i = beta[i, ],
      gamma = gamma,
      lambda = lambda,
      u = u,
      int_bounds_i = int_bounds[[i]],
      k_T_i = k_T[i]
    )
    
    if (!is.finite(li)) {
      cat("Non-finite ll at i =", i)
    }
    
    total_log_density <- total_log_density + li
    
    # density_vec[i] = total_log_density
  }
  
  return(total_log_density)
}

## Function to sample alpha, px1 matrix of coeff. for demo covs.:
sample_alpha <- function(alpha_curr, beta, T, delta, X, gamma, lambda, u, int_bounds, k_T,
                         mu_alpha, Sigma_alpha, Sigma_prop_alpha) {
  
  alpha_prop <- alpha_curr
  # Generate proposal centered at current value
  
  lsurv_curr <- lsurv_all(T, delta, X, alpha_curr, beta, gamma, lambda, u, int_bounds, k_T)
  
  for (j in 1:p) {
    alpha_prop[j] <- rnorm(1, mean = alpha_curr[j], sd = Sigma_prop_alpha[j, j]) 
    
    # Compute log-posterior at current value
    lsurv_curr <- lsurv_curr
    lprior_curr <- dnorm(alpha_curr[j], mean = mu_alpha[j], sd = Sigma_alpha[j, j], log = TRUE)
    lpost_curr <- lsurv_curr + lprior_curr
    
    # Compute log-posterior at proposal value
    lsurv_prop <- lsurv_all(T, delta, X, alpha_prop, beta, gamma, lambda, u, int_bounds, k_T)
    lprior_prop <- dnorm(alpha_prop[j], mean = mu_alpha[j], sd = Sigma_alpha[j, j], log = TRUE)
    lpost_prop = lsurv_prop + lprior_prop
    
    # Accept/reject
    log_ratio <- lpost_prop - lpost_curr
    if (log(runif(1)) < log_ratio) {
      alpha_curr[j] = alpha_prop[j]
    } else {
      alpha_prop[j] = alpha_curr[j]
    }
  }
  return(alpha_curr)
}

## Function to sample gamma, a scalar:
sample_gamma <- function(gamma_curr, beta, T, delta, X, alpha, lambda, u, int_bounds, k_T,
                         mu_gamma, sigma_sq_gamma, sigma_sq_prop) {
  
  # Generate proposal centered at current value
  gamma_prop = rnorm(1, mean = gamma_curr, sd = sqrt(sigma_sq_prop))

  # Compute log-posterior at current value
  lsurv_curr <- lsurv_all(T, delta, X, alpha, beta, gamma_curr, lambda, u, int_bounds, k_T)
  # print(paste("lsurv_all at current", lsurv_curr))
  lprior_curr <- dnorm(gamma_curr, mean = mu_gamma, sd = sqrt(sigma_sq_gamma), log = TRUE)
  # print(paste("lprior at current", lprior_curr))
  lpost_curr <- lsurv_curr + lprior_curr
  
  # Compute log-posterior at proposal value
  lsurv_prop <- lsurv_all(T, delta, X, alpha, beta, gamma_prop, lambda, u, int_bounds, k_T)
  # print(paste("lsurv_all at proposal", lsurv_prop))
  lprior_prop <- dnorm(gamma_prop, mean = mu_gamma, sd = sqrt(sigma_sq_gamma), log = TRUE)
  # print(paste("lprior at proposal", lprior_prop))
  lpost_prop <- lsurv_prop + lprior_prop
  
  # Accept/reject
  log_ratio <- lpost_prop - lpost_curr
  
  if(log(runif(1)) < log_ratio) {
    return(gamma_prop)
  } else {
    return(gamma_curr)
  }
}

## Function to adaptively tune parameters
pilot_adapt <- function(tuning_param, accept_rate) {
  
  ## adjust tuning parameter by scaling existing parameter based on current acceptance rate
  if (accept_rate >= 0.90) {
    tuning_param <- tuning_param * 1.3
  } else if ((accept_rate >= 0.75 ) & (accept_rate < 0.90 )) {
    tuning_param <- tuning_param * 1.2
  } else if ((accept_rate >= 0.45 ) & (accept_rate < 0.75 )) {
    tuning_param <- tuning_param * 1.1
  } else if ((accept_rate <= 0.25 ) & (accept_rate > 0.15 )) {
    tuning_param <- tuning_param * 0.9
  } else if ((accept_rate <= 0.15 ) & (accept_rate > 0.10 )) {
    tuning_param <- tuning_param * 0.8
  } else if (accept_rate <= 0.10) {
    tuning_param <- tuning_param * 0.7
  }
  
  return(tuning_param)
}

### Define hyperparameters

## beta_spike hyperparameters
mu_spike  <- 0
nu2_spike <- 10

## mu_beta hyperparameters; flat distribution
b0 <- 0 
b1 <- 0
tau0_sq <- 100
tau1_sq <- 100

## lambda_k hyperparameters
a_k <- .001
b_k <- .001

## gamma hyperparameters
mu_gamma <- 0
sigma_sq_gamma <- 1

## alpha hyperparameters
mu_alpha <- matrix(0, nrow = p, ncol = 1)
Sigma_alpha <- diag(100, p)

## Sigma_beta hyperparameters
nu0 <- 4 # must be larger than p, where p = 2 for beta
S0 <- diag(2) # pxp identity matrix

## sigma_sq hyperparameters
a <- .001
b <- .001

## Check initialization in environment:
# beta, lambda, mu_beta, Sigma_beta, alpha
# missing: gamma

## Remaining parameters to initialize
gamma <- 0.001 # Cannot be 0 otherwise div0 errors will occur in computation of log surv density; choose small number to assume no assoc.

## Initial proposal tuning parameters
Sigma_prop_alpha <- diag(0.1, p)
sigma_sq_prop <- .001

K <- length(lambda)

S <- 255000
burn_in <- 5000
n_keep <- 25000
adapt_interval <- 1000

thin <- (S - burn_in) / n_keep

S_keep <- n_keep  # number of stored posterior draws

### MCMC SAMPLING
set.seed(847)

## Track acceptance counts for relevant parameters (adaptive tuning)
accept_alpha <- rep(0, p)
adapt_counter_alpha <- rep(0, p)

accept_gamma <- 0
adapt_counter_gamma <- 0

Sigma_prop_beta <- diag(c(0.1, 0.1))
accept_beta <- 0
adapt_counter_beta <- 0

## Set up parameter storage w/appropriate dims
BETA <- array(NA, dim = c(N, 3, S))
MU_BETA <- matrix(NA, nrow = S, ncol = 2)
SIGMA_BETA <- array(NA, dim = c(2, 2, S))
SIGMA_SQ  <- numeric(S)
ALPHA <- matrix(NA, nrow = S, ncol = p)
GAMMA <- numeric(S)
LAMBDA <- matrix(NA, nrow = S, ncol = K)

output_dir <- "/hpc/group/berchucklab/kgs35/outputs"

param_csv <- file.path(output_dir, "mcmc_params_ALL.csv")

m <- 0

for(s in 1:S) {
  
  #1. Sample beta_i for each patient:
  for (i in 1:N) {
    beta_prev_i <- beta[i, ]
    
    beta[i, ] <- sample_beta_i(
      beta_i_curr = beta[i, ],
      y_i = Y[[i]],
      t_i = t[[i]],
      T_i = time_data$time_final[time_data$obs == i],
      delta_i = time_data$event[time_data$obs == i],
      X_i = X[i, , drop = FALSE], # dims are 1 x p
      alpha = alpha,
      gamma = gamma,
      lambda = lambda,
      u = u,
      mu_beta = mu_beta,
      kappa_i = kappa[[i]],
      int_bounds_i = int_bounds[[i]],
      k_T_i = k_T[i],
      Sigma_beta = Sigma_beta,
      sigma2 = sigma_sq,
      Sigma_prop = Sigma_prop_beta
    )
    
    if (s <= burn_in) {
      accepted_i <- any(beta[i, 1:2] != beta_prev_i[1:2])
      accept_beta <- accept_beta + as.integer(accepted_i)
      adapt_counter_beta <- adapt_counter_beta + 1L
    }
  }
  
  ## Adapt tuning parameter for alpha if within burn-in
  if (s <= burn_in && s %% adapt_interval == 0) {
    
    accept_rate_beta <- accept_beta / adapt_counter_beta
    
    mult <- pilot_adapt(1, accept_rate_beta)
    Sigma_prop_beta <- Sigma_prop_beta * mult
    
    # reset
    accept_beta <- 0L
    adapt_counter_beta <- 0L
  }
  
  #2. Sample beta_spike_i for each patient
  for (i in 1:N) {
    beta[i, ] <- sample_beta_spike_i(
      y_i = Y[[i]],
      t_i = t[[i]],
      active_idx_i = active_idx[[i]],            
      beta_i = beta[i, ],
      sigma_sq = sigma_sq,
      mu_spike = mu_spike,
      nu2_spike = nu2_spike
    )
  }
  
  #3. Sample mu_beta:
  mu_beta0 <- sample_mu_beta0(
    beta = beta, 
    Sigma_beta = Sigma_beta, 
    mu_beta1 = mu_beta[1, 2], 
    tau0_sq = tau0_sq, 
    b0 = b0)
  
  mu_beta1 <- sample_mu_beta1(
    beta = beta, 
    Sigma_beta = Sigma_beta, 
    mu_beta0 = mu_beta[1, 1], 
    tau1_sq = tau1_sq, 
    b1 = b1)
  
  mu_beta <- matrix(c(mu_beta0, mu_beta1), nrow = 1, ncol = 2)
  
  #4. Sample Sigma_beta:
  Sigma_beta <- sample_Sigma_beta(nu0, S0, beta, mu_beta)
  
  #5. Sample sigma^2 
  sigma_sq <- sample_sigma_sq(
    a = a, 
    b = b,
    beta = beta, 
    Y = Y, 
    t = t, 
    kappa = kappa)
  
  #6. Sample lambda:
  for (k in 1:K) {
    lambda[k] <- sample_lambda_k(
      k = k, u = u, time_data = time_data, X = X,
      alpha = alpha, beta = beta, gamma = gamma,
      int_bounds = int_bounds, max_t = 6.01, a_k = a_k, b_k = b_k
    )
  }
  
  #7. Sample alpha:
  alpha_prev <- alpha
  alpha <- sample_alpha(
    alpha_curr = alpha,
    beta = beta,
    T = time_data$time_final,
    delta = time_data$event,
    X = X,
    gamma = gamma,
    lambda = lambda,
    u = u,
    int_bounds = int_bounds,
    k_T = k_T,
    mu_alpha = mu_alpha,
    Sigma_alpha = Sigma_alpha,
    Sigma_prop_alpha = Sigma_prop_alpha
  )
  
  ## Adapt tuning parameter for alpha if within burn-in
  if (s <= burn_in) {
    # Check for acceptances 
    for (j in 1:p) {
      if (alpha[j] != alpha_prev[j]) {
        accept_alpha[j] <- accept_alpha[j] + 1
      }
      adapt_counter_alpha[j] <- adapt_counter_alpha[j] + 1
    }
    
    # Every adapt_interval iterations, adapt tuning
    if (s %% adapt_interval == 0) {
      accept_rate_alpha <- accept_alpha / adapt_counter_alpha
      for (j in 1:p) {
        Sigma_prop_alpha[j, j] <- pilot_adapt(Sigma_prop_alpha[j, j], accept_rate_alpha[j])
      }
      # Reset acceptance counters for next phase
      accept_alpha <- rep(0, p)
      adapt_counter_alpha <- rep(0, p)
    }
  }
  
  #8. Sample gamma:
  gamma_prev <- gamma
  gamma <- sample_gamma(
    gamma_curr = gamma,
    beta = beta,
    T = time_data$time_final,
    delta = time_data$event,
    X = X,
    alpha = alpha,
    lambda = lambda,
    u = u,
    int_bounds = int_bounds,
    k_T = k_T,
    mu_gamma = mu_gamma,
    sigma_sq_gamma = sigma_sq_gamma,
    sigma_sq_prop = sigma_sq_prop
  )
  
  ## Adapt tuning parameter for gamma if within burn-in
  if (s <= burn_in) {
    # Check for acceptances 
    if (gamma != gamma_prev) {
      accept_gamma <- accept_gamma + 1
    }
    adapt_counter_gamma <- adapt_counter_gamma + 1
    
    # Every adapt_interval iterations, adapt tuning
    if (s %% adapt_interval == 0) {
      accept_rate_gamma <- accept_gamma / adapt_counter_gamma
      sigma_sq_prop <- pilot_adapt(sigma_sq_prop, accept_rate_gamma)
      
      # Reset acceptance counters for next phase
      accept_gamma <- 0
      adapt_counter_gamma <- 0
    }
  }
  
  #9. Store all unknowns: 
  
  keep <- (s > burn_in) && ((s - burn_in) %% thin == 0L)
  
  if (keep) {
    m <- m + 1L
    
    # 9. Store all unknowns (to thinned idx, m)
    BETA[, , m] <- beta
    MU_BETA[m, ] <- mu_beta
    SIGMA_BETA[, , m] <- Sigma_beta
    SIGMA_SQ[m] <- sigma_sq
    ALPHA[m, ] <- t(alpha)
    GAMMA[m] <- gamma
    LAMBDA[m, ] <- lambda
    
    # Optional: write CSV every so often (based on m, not s)
    if (m %% 50L == 0L) {
      cat("Iter:", s, " Stored draw:", m, "\n")
      
      alpha_vec  <- as.numeric(alpha)
      names(alpha_vec) <- paste0("alpha_", seq_len(p))
      
      lambda_vec <- as.numeric(lambda)
      names(lambda_vec) <- paste0("lambda_", seq_len(K))
      
      row_df <- data.frame(
        iter = s,          
        draw = m, 
        sigma_sq = sigma_sq,
        gamma = gamma,
        mu_beta0 = mu_beta[1, 1],
        mu_beta1 = mu_beta[1, 2],
        t(alpha_vec),
        t(lambda_vec),
        check.names = FALSE
      )
      
      write.table(
        row_df,
        file = param_csv,
        sep = ",",
        row.names = FALSE,
        col.names = !file.exists(param_csv),
        append = file.exists(param_csv)
      )
    }
    
    if (m >= S_keep) break
  }
  
  if (s %% 5000L == 0L) {
    cat("Iteration:", s, " (stored:", m, ")\n")
  }
}

## Save MCMC outputs
save(BETA, MU_BETA, SIGMA_BETA, SIGMA_SQ, ALPHA, GAMMA, LAMBDA,
     file = file.path(output_dir, "mcmc_outputs_ALL_030726.RData"))
