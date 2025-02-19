run_glm <- function(beta, x, y, alpha) {
  tic()

  glm_fit <- glm(y ~ x, family = binomial)
  time_elapsed <- toc(quiet = TRUE)

  glm_beta <- coef(glm_fit)
  glm_se <- sqrt(diag(vcov(glm_fit)))
  glm_ci <- confint(glm_fit)

  return(list(
    sol = glm_beta,
    se = glm_se,
    ci = glm_ci
  ))
}
