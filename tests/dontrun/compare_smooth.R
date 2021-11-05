# ---- negative log EL functions ----

# with original check function
neg_logel_qr_orig <- function(beta, y, X, tau,
                              max_iter = 100, rel_tol = 1e-7) {
  G <- qr_evalG_R(y, X, tau, beta)
  genel <- flexEL::GenEL$new(n_obs = nrow(X), n_eqs = ncol(X))
  genel$max_iter <- max_iter
  genel$rel_tol <- rel_tol
  -genel$logel(G)
}

# with flexEL's smoothed check function
neg_logel_qr_smooth <- function(beta, y, X, tau, alpha,
                                max_iter = 100, rel_tol = 1e-7) {
  G <- qr_evalG_smooth_R(y, X, tau, beta, alpha)
  genel <- flexEL::GenEL$new(n_obs = nrow(X), n_eqs = ncol(X))
  genel$max_iter <- max_iter
  genel$rel_tol <- rel_tol
  -genel$logel(G)
}

# with Zheng's smoothed check function
neg_logel_qr_sfun <- function(beta, y, X, tau, alpha,
                              max_iter = 100, rel_tol = 1e-7) {
  G <- qr_evalG_sfun(y, X, beta, tau, alpha)
  genel <- flexEL::GenEL$new(n_obs = nrow(X), n_eqs = ncol(X))
  genel$max_iter <- max_iter
  genel$rel_tol <- rel_tol
  -genel$logel(G)
}

# with shifted Zheng's smoothed check function
neg_logel_qr_zfun <- function(beta, y, X, tau, alpha,
                              max_iter = 100, rel_tol = 1e-7) {
  G <- qr_evalG_zfun(y, X, beta, tau, alpha)
  genel <- flexEL::GenEL$new(n_obs = nrow(X), n_eqs = ncol(X))
  genel$max_iter <- max_iter
  genel$rel_tol <- rel_tol
  -genel$logel(G)
}

