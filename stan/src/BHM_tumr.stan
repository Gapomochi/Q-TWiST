// tumor_growth_curve.stan

// ------------------------------------------------------------
// DATA
// ------------------------------------------------------------
data {
  int<lower=1> N;                      // total observations
  int<lower=1> N_subj;                 // number of subjects
  int<lower=1> K;                      // number of treatment groups

  array[N] int<lower=1, upper=N_subj> id;        // subject index
  array[N_subj] int<lower=1, upper=K> trt_subj;  // subject's treatment group

  vector[N] y;                         // tumor volume (or log volume)
  vector[N] t;                         // time
}

// ------------------------------------------------------------
// PARAMETERS
// ------------------------------------------------------------
parameters {
  // treatment-level means (population)
  vector[K] Int;        // mean intercept for each treatment
  vector[K] Slope;      // mean slope for each treatment

  // level-1 noise
  real<lower=0> sdLevel1;

  // random effects (intercept, slope) - non-centered
  vector<lower=0>[2] tau;                 // SDs of RE (intercept, slope)
  cholesky_factor_corr[2] L_corr;         // correlation (Cholesky)
  array[N_subj] vector[2] z;              // i.i.d. N(0,1)
}

// ------------------------------------------------------------
// TRANSFORMED PARAMETERS
// ------------------------------------------------------------
transformed parameters {
  array[N_subj] vector[2] beta;           // subject-specific (int, slope)
  matrix[2,2] L = diag_pre_multiply(tau, L_corr);

  for (i in 1:N_subj) {
    vector[2] mu_i;
    int g = trt_subj[i];

    mu_i[1] = Int[g];
    mu_i[2] = Slope[g];

    beta[i] = mu_i + L * z[i];
  }
}

// ------------------------------------------------------------
// MODEL
// ------------------------------------------------------------
model {
  // Priors (you can tune later)
  Int   ~ normal(0, 10);
  Slope ~ normal(0, 10);

  sdLevel1 ~ normal(0, 5);          // half-normal via <lower=0>

  tau    ~ normal(0, 5);            // half-normal
  L_corr ~ lkj_corr_cholesky(2);    // correlation prior
  for (i in 1:N_subj) z[i] ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    y[n] ~ normal(beta[id[n]][1] + beta[id[n]][2] * t[n], sdLevel1);
  }
}

// ------------------------------------------------------------
// GENERATED QUANTITIES
// ------------------------------------------------------------
generated quantities {
  matrix[K,K] IntDiff;
  matrix[K,K] SlopeDiff;

  // pairwise differences: (row - col) = (i - j)
  for (i in 1:K) {
    for (j in 1:K) {
      IntDiff[i,j]   = Int[i]   - Int[j];
      SlopeDiff[i,j] = Slope[i] - Slope[j];
    }
  }
}
