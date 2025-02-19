run_glm <- function(beta, x, y, alpha) {
  tic()

  glm_fit <- glm(y ~ x, family = binomial)
  time_elapsed <- toc(quiet = TRUE)

  solution <- coef(glm_fit)
  beta0_sol = solution[[1]]
  beta1_sol = solution[[3]]

  se_beta <- sqrt(diag(vcov(glm_fit)))
  ci_l = solution - qnorm(1 - alpha/2) * se_beta
  ci_u = solution + qnorm(1 - alpha/2) * se_beta

  beta0_ci_l <- ci_l[[1]]
  beta0_ci_u <- ci_u[[1]]
  beta1_ci_l <- ci_l[[3]]
  beta1_ci_u <- ci_u[[3]]

  return(list(
    beta0_sol = beta0_sol,
    beta1_sol = beta1_sol,
    se_beta0 = se_beta[[1]],
    se_beta1 = se_beta[[3]],
    beta0_ci_l = beta0_ci_l,
    beta0_ci_u = beta0_ci_u,
    beta1_ci_l = beta1_ci_l,
    beta1_ci_u = beta1_ci_u,
    time_elapsed = time_elapsed$toc - time_elapsed$tic
  ))
}
