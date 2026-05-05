riley_continuous_shrinkage <- function(n, p, r2_adj) {
  if (any(n <= p + 1)) stop("n must be greater than p + 1.")
  if (p < 2) stop("p must be at least 2.")
  if (r2_adj <= 0 || r2_adj >= 1) stop("r2_adj must be between 0 and 1.")

  r2_app <- ((n - p - 1) * r2_adj + p) / (n - 1)
  1 + (p - 2) / (n * log(1 - r2_app))
}

riley_continuous_n_shrinkage <- function(p, r2_adj, target_sc = 0.9,
                                         lower = p + 2, upper = 1e7) {
  if (p < 2) stop("p must be at least 2.")
  if (r2_adj <= 0 || r2_adj >= 1) stop("r2_adj must be between 0 and 1.")
  if (target_sc <= 0 || target_sc >= 1) stop("target_sc must be between 0 and 1.")

  f <- function(n) riley_continuous_shrinkage(n, p, r2_adj) - target_sc

  lo <- max(ceiling(lower), p + 2)
  hi <- ceiling(upper)
  if (f(lo) >= 0) return(lo)
  if (f(hi) < 0) stop("Upper bound too small to reach the target shrinkage.")

  ceiling(uniroot(f, lower = lo, upper = hi)$root)
}

riley_continuous_n_r2diff <- function(p, r2_adj, delta = 0.05) {
  if (p < 1) stop("p must be positive.")
  if (r2_adj < 0 || r2_adj >= 1) stop("r2_adj must be in [0, 1).")
  if (delta <= 0 || delta >= 1) stop("delta must be between 0 and 1.")
  ceiling(1 + p * (1 - r2_adj) / delta)
}

riley_continuous_mmoe_sigma <- function(n, p, alpha = 0.05) {
  if (n <= p + 1) stop("n must be greater than p + 1.")
  df <- n - p - 1
  q_hi <- stats::qchisq(1 - alpha / 2, df)
  q_lo <- stats::qchisq(alpha / 2, df)
  sqrt(max(df / q_hi, df / q_lo))
}

riley_continuous_n_sigma <- function(p, target_mmoe = 1.1, alpha = 0.05,
                                     lower = p + 2, upper = 1e7) {
  if (p < 0) stop("p must be non-negative.")
  if (target_mmoe <= 1) stop("target_mmoe must be greater than 1.")

  f <- function(n) riley_continuous_mmoe_sigma(n, p, alpha) - target_mmoe

  lo <- ceiling(lower)
  hi <- ceiling(upper)
  if (f(lo) <= 0) return(lo)
  if (f(hi) > 0) stop("Upper bound too small to reach the target MMOE.")

  ceiling(uniroot(f, lower = lo, upper = hi)$root)
}

riley_continuous_intercept_ci <- function(n, p, alpha_hat, sigma2_null,
                                          r2_adj, alpha = 0.05) {
  if (n <= p + 1) stop("n must be greater than p + 1.")
  if (sigma2_null <= 0) stop("sigma2_null must be positive.")
  if (r2_adj < 0 || r2_adj >= 1) stop("r2_adj must be in [0, 1).")

  df <- n - p - 1
  se <- sqrt(sigma2_null * (1 - r2_adj) / n)
  crit <- stats::qt(1 - alpha / 2, df)
  c(lower = alpha_hat - crit * se, upper = alpha_hat + crit * se)
}

riley_continuous_min_sample_size <- function(p, r2_adj,
                                             target_sc = 0.9,
                                             delta_r2 = 0.05,
                                             target_mmoe = 1.1,
                                             alpha = 0.05) {
  n_i <- riley_continuous_n_shrinkage(p, r2_adj, target_sc)
  n_ii <- riley_continuous_n_r2diff(p, r2_adj, delta_r2)
  n_iii <- riley_continuous_n_sigma(p, target_mmoe, alpha)

  list(
    criterion_i = n_i,
    criterion_ii = n_ii,
    criterion_iii = n_iii,
    minimum = max(n_i, n_ii, n_iii),
    subjects_per_parameter = max(n_i, n_ii, n_iii) / p
  )
}

# Example from the paper:
# riley_continuous_min_sample_size(p = 25, r2_adj = 0.2)
