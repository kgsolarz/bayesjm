library(Rcpp)
library(RcppArmadillo)
library(readr)
library(coda)
library(MASS)
library(dplyr)
library(forcats)
library(survival)
library(invgamma)
library(ggplot2)
library(tidyr)

## Compile the RcppArmadillo samplers -- adjust path to match where this
## file lives on the cluster relative to MCMC_Samplers.cpp
Rcpp::sourceCpp("/hpc/group/berchucklab/kgs35/scripts/MCMC_Samplers.cpp")


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

### END OF DATA PRE-PROCESSING

### Pre-flatten the objects that get reused every MCMC iteration into the
### exact shapes the RcppArmadillo samplers expect, so we're not re-deriving
### them (or re-filtering time_data) inside the S x N hot loop.
##
## Note: nothing here needs reshaping beyond type coercion -- RcppArmadillo's
## arma::vec/arma::mat converters accept R matrices (e.g. alpha's p x 1 shape,
## mu_beta's 1x2 shape) transparently, so alpha/mu_alpha/mu_beta are passed
## to the _cpp functions completely unchanged below.

T_vec <- time_data$time_final           # already ordered 1:N to match beta/Y/t
delta_vec <- as.integer(time_data$event)
X_mat <- as.matrix(X)
k_T_int <- as.integer(k_T)

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

### MCMC SAMPLING (RcppArmadillo-driven)
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

param_csv <- file.path(output_dir, "mcmc_params_ALL_cpp.csv")

m <- 0

for(s in 1:S) {

  #1. Sample beta_i for each patient:
  for (i in 1:N) {
    beta_prev_i <- beta[i, ]

    beta[i, ] <- sample_beta_i_cpp(
      beta_i_curr = beta[i, ],
      y_i = Y[[i]],
      t_i = t[[i]],
      T_i = T_vec[i],
      delta_i = delta_vec[i],
      X_i = X_mat[i, , drop = FALSE], # dims are 1 x p
      alpha = alpha,
      gamma = gamma,
      lambda = lambda,
      u = u,
      mu_beta = mu_beta,
      kappa_i = kappa[[i]],
      int_bounds_i = int_bounds[[i]],
      k_T_i = k_T_int[i],
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

  ## Adapt tuning parameter for beta if within burn-in
  if (s <= burn_in && s %% adapt_interval == 0) {

    accept_rate_beta <- accept_beta / adapt_counter_beta

    mult <- pilot_adapt_cpp(1, accept_rate_beta)
    Sigma_prop_beta <- Sigma_prop_beta * mult

    # reset
    accept_beta <- 0L
    adapt_counter_beta <- 0L
  }

  #2. Sample beta_spike_i for each patient
  for (i in 1:N) {
    beta[i, ] <- sample_beta_spike_i_cpp(
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
  mu_beta0 <- sample_mu_beta0_cpp(
    beta = beta,
    Sigma_beta = Sigma_beta,
    mu_beta1 = mu_beta[1, 2],
    tau0_sq = tau0_sq,
    b0 = b0)

  mu_beta1 <- sample_mu_beta1_cpp(
    beta = beta,
    Sigma_beta = Sigma_beta,
    mu_beta0 = mu_beta[1, 1],
    tau1_sq = tau1_sq,
    b1 = b1)

  mu_beta <- matrix(c(mu_beta0, mu_beta1), nrow = 1, ncol = 2)

  #4. Sample Sigma_beta:
  Sigma_beta <- sample_Sigma_beta_cpp(nu0, S0, beta, mu_beta)

  #5. Sample sigma^2
  sigma_sq <- sample_sigma_sq_cpp(
    a = a,
    b = b,
    beta = beta,
    Y = Y,
    t = t,
    kappa = kappa)

  #6. Sample lambda:
  for (k in 1:K) {
    lambda[k] <- sample_lambda_k_cpp(
      k = k, u = u, T = T_vec, delta = delta_vec, X = X_mat,
      alpha = alpha, beta = beta, gamma = gamma,
      int_bounds = int_bounds, a_k = a_k, b_k = b_k
    )
  }

  #7. Sample alpha:
  alpha_prev <- alpha
  alpha <- sample_alpha_cpp(
    alpha_curr = alpha,
    beta = beta,
    T = T_vec,
    delta = delta_vec,
    X = X_mat,
    gamma = gamma,
    lambda = lambda,
    u = u,
    int_bounds = int_bounds,
    k_T = k_T_int,
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
        Sigma_prop_alpha[j, j] <- pilot_adapt_cpp(Sigma_prop_alpha[j, j], accept_rate_alpha[j])
      }
      # Reset acceptance counters for next phase
      accept_alpha <- rep(0, p)
      adapt_counter_alpha <- rep(0, p)
    }
  }

  #8. Sample gamma:
  gamma_prev <- gamma
  gamma <- sample_gamma_cpp(
    gamma_curr = gamma,
    beta = beta,
    T = T_vec,
    delta = delta_vec,
    X = X_mat,
    alpha = alpha,
    lambda = lambda,
    u = u,
    int_bounds = int_bounds,
    k_T = k_T_int,
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
      sigma_sq_prop <- pilot_adapt_cpp(sigma_sq_prop, accept_rate_gamma)

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
     file = file.path(output_dir, "mcmc_outputs_ALL_cpp_030726.RData"))
