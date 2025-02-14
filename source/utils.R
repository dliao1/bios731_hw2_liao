gen_logit_data <- function(n, beta0, beta1) {
  x <- rnorm(n)
  x <- cbind(1, x)
  beta_true <- matrix(c(beta0, beta1), ncol = 1)

  pi <- exp(x %*% beta_true) / (1 + exp(x %*% beta_true))

  y <- rbinom(n, size = 1, prob = pi)

  list(
    x = x,
    y = y
  )
}
