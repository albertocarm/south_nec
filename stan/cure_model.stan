// Promotion-time cure model, Weibull promotion time, time-dependent chemotherapy
functions {
  real rcs2(real x, vector k) {
    real kd2 = square(k[3] - k[1]);
    return (pow(fmax(x - k[1], 0), 3)
            - pow(fmax(x - k[2], 0), 3) * (k[3] - k[1]) / (k[3] - k[2])
            + pow(fmax(x - k[3], 0), 3) * (k[2] - k[1]) / (k[3] - k[2])) / kd2;
  }
}
data {
  int<lower=1> N;
  int<lower=0> K;
  int<lower=1> Kc;
  int<lower=0> L;
  matrix[N, K] X;                    // clonogenic covariates
  matrix[N, Kc] W;                   // chemotherapy, clonogenic
  matrix[N, L] Z;                    // kinetic covariates
  matrix[N, 1] V;                    // chemotherapy, kinetic
  vector<lower=0>[N] t;
  array[N] int<lower=0, upper=1> d;
  array[N] int<lower=0, upper=2> path;  // 0 unexposed, 1 exposed from surgery, 2 exposed from s
  vector<lower=0>[N] s;
  vector<lower=0>[K] sd_beta;
  vector<lower=0>[Kc] sd_delta;
  vector<lower=0>[L] sd_gamma;
  real<lower=0> sd_eta;
  real<lower=0> sd_int;
  int<lower=0> Nk;                   // missing Ki-67
  array[Nk] int<lower=1, upper=N> ik;
  vector[Nk] k_lo;
  vector[Nk] k_hi;
  vector[3] k_knots;
  vector[4] k_std;
  array[5] int<lower=0> ck;          // Ki-67 columns: X linear, X spline, Z linear, Z spline, W
  int<lower=0> Na;                   // missing age
  array[Na] int<lower=1, upper=N> ia;
  vector[Na] a_lo;
  vector[Na] a_hi;
  vector[3] a_knots;
  vector[4] a_std;
  array[2] int<lower=0> ca;          // age columns: X linear, X spline
}
parameters {
  real beta0;
  vector[K] beta;
  vector[Kc] delta;
  real gamma0;
  vector[L] gamma;
  vector[1] eta;
  real<lower=0> alpha;
  vector<lower=k_lo, upper=k_hi>[Nk] ki_miss;
  vector<lower=a_lo, upper=a_hi>[Na] age_miss;
}
transformed parameters {
  matrix[N, K] Xf;
  matrix[N, Kc] Wf;
  matrix[N, L] Zf;
  vector[N] lt0;                     // log theta before chemotherapy
  vector[N] lt1;                     // log theta after chemotherapy
  vector[N] ll0;                     // log lambda before chemotherapy
  vector[N] ll1;                     // log lambda after chemotherapy
  Xf = X;
  Wf = W;
  Zf = Z;
  for (j in 1:Nk) {
    int i = ik[j];
    real za = (ki_miss[j] - k_std[1]) / k_std[2];
    real zb = (rcs2(ki_miss[j], k_knots) - k_std[3]) / k_std[4];
    if (ck[1] > 0) Xf[i, ck[1]] = za;
    if (ck[2] > 0) Xf[i, ck[2]] = zb;
    if (ck[3] > 0) Zf[i, ck[3]] = za;
    if (ck[4] > 0) Zf[i, ck[4]] = zb;
    if (ck[5] > 0) Wf[i, ck[5]] = za;
  }
  for (j in 1:Na) {
    int i = ia[j];
    real za = (age_miss[j] - a_std[1]) / a_std[2];
    real zb = (rcs2(age_miss[j], a_knots) - a_std[3]) / a_std[4];
    if (ca[1] > 0) Xf[i, ca[1]] = za;
    if (ca[2] > 0) Xf[i, ca[2]] = zb;
  }
  lt0 = beta0 + Xf * beta;
  lt1 = lt0 + Wf * delta;
  ll0 = gamma0 + Zf * gamma;
  ll1 = ll0 + V * eta;
}
model {
  beta0 ~ normal(0, sd_int);
  gamma0 ~ normal(0, sd_int);
  beta ~ normal(0, sd_beta);
  delta ~ normal(0, sd_delta);
  gamma ~ normal(0, sd_gamma);
  eta ~ normal(0, sd_eta);
  alpha ~ gamma(2, 1);
  ki_miss ~ normal(k_std[1], k_std[2]);
  age_miss ~ normal(a_std[1], a_std[2]);
  for (i in 1:N) {
    real l0 = exp(ll0[i]);
    real l1 = exp(ll1[i]);
    real th0 = exp(lt0[i]);
    real th1 = exp(lt1[i]);
    real H;
    real lh;
    if (path[i] == 0) {
      real L0 = pow(t[i] / l0, alpha);
      H = th0 * (1 - exp(-L0));
      lh = lt0[i] + log(alpha) - log(l0) + (alpha - 1) * log(t[i] / l0) - L0;
    } else if (path[i] == 1) {
      real L1 = pow(t[i] / l1, alpha);
      H = th1 * (1 - exp(-L1));
      lh = lt1[i] + log(alpha) - log(l1) + (alpha - 1) * log(t[i] / l1) - L1;
    } else {
      real L0 = pow(s[i] / l0, alpha);
      real dL1 = pow(t[i] / l1, alpha) - pow(s[i] / l1, alpha);
      H = th0 * (1 - exp(-L0)) + th1 * exp(-L0) * (1 - exp(-dL1));
      lh = lt1[i] - L0 - dL1 + log(alpha) - log(l1) + (alpha - 1) * log(t[i] / l1);
    }
    target += -H;
    if (d[i] == 1) target += lh;
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] cure;
  for (i in 1:N) {
    real l0 = exp(ll0[i]);
    real l1 = exp(ll1[i]);
    real th0 = exp(lt0[i]);
    real th1 = exp(lt1[i]);
    real H;
    real Hinf;
    real lh;
    if (path[i] == 0) {
      real L0 = pow(t[i] / l0, alpha);
      H = th0 * (1 - exp(-L0));
      Hinf = th0;
      lh = lt0[i] + log(alpha) - log(l0) + (alpha - 1) * log(t[i] / l0) - L0;
    } else if (path[i] == 1) {
      real L1 = pow(t[i] / l1, alpha);
      H = th1 * (1 - exp(-L1));
      Hinf = th1;
      lh = lt1[i] + log(alpha) - log(l1) + (alpha - 1) * log(t[i] / l1) - L1;
    } else {
      real L0 = pow(s[i] / l0, alpha);
      real dL1 = pow(t[i] / l1, alpha) - pow(s[i] / l1, alpha);
      H = th0 * (1 - exp(-L0)) + th1 * exp(-L0) * (1 - exp(-dL1));
      Hinf = th0 * (1 - exp(-L0)) + th1 * exp(-L0);
      lh = lt1[i] - L0 - dL1 + log(alpha) - log(l1) + (alpha - 1) * log(t[i] / l1);
    }
    cure[i] = exp(-Hinf);
    log_lik[i] = (d[i] == 1 ? lh : 0) - H;
  }
}
