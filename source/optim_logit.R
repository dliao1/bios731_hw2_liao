neg_log_likelihood <- function(beta, x, y) {
  pi <- exp(x %*% beta) / (1 + exp(x %*% beta))
  return(-sum(y * log(pi) + (1 - y) * log(1 - pi)))  # Negative log-likelihood
}

run_optim <- function(beta, x, y, alpha) {
  tic()
  optim_fit <- optim(par = beta,
                     fn = neg_log_likelihood,
                     x = x,
                     y = y,
                     method = "BFGS",
                     hessian = TRUE)

  time_elapsed <- toc(quiet = TRUE)

  solution = optim_fit$par
  se_beta <- sqrt(diag(solve(optim_fit$hessian)))
  ci_l <- optim_fit$par - qnorm(1 - alpha/2) * se_beta
  ci_u <- optim_fit$par + qnorm(1 - alpha/2) * se_beta

  beta0_ci_l <- ci_l[[1]]
  beta0_ci_u <- ci_u[[1]]
  beta1_ci_l <- ci_l[[2]]
  beta1_ci_u <- ci_u[[2]]


  return(list(beta0_sol = solution[[1]],
              beta1_sol = solution[[2]],
              se_beta0 = se_beta[[1]],
              se_beta1 = se_beta[[2]],
              beta0_ci_l = beta0_ci_l,
              beta0_ci_u = beta0_ci_u,
              beta1_ci_l = beta1_ci_l,
              beta1_ci_u = beta1_ci_u,
              time_elapsed = time_elapsed$toc - time_elapsed$tic))

}
