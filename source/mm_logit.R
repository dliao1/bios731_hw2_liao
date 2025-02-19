mm_logit <- function(theta, x, y, alpha, max_iter = 200) {
  iter = 1
  tol = Inf

  objective = rep(NA, length = max_iter)
  p <- ncol(x)
  tic()

  for(iter in 1:max_iter){

    pi <- exp(x %*% theta) / (1 + exp(x %*% theta))

    theta_new <- theta

    for (j in 1:p) {

      root_func <- function(theta_j) {
        sum1 <- sum(y * x[, j])

        sum2 <- sum((pi * x[, j] * exp(-p * x[, j] * theta[j])) * exp(p * x[, j] * theta_j))

        return(sum1 - sum2)
      }

      # Solve for theta_j using uniroot!
      solution <- uniroot(root_func, interval = c(-10, 10), tol = 1e-6)

      theta_new[j] <- solution$root
    }

    if (iter > 1){
      tol = abs(sum(theta_new - theta))
    }

    if (tol < 1e-6) {
      break
    }

    theta <- theta_new

  }

  time_elapsed <- toc(quiet = TRUE)

  # Gets hessian for inference, similar to newton
  hessian_ls = as.list(rep(NA, length.out = length(nrow(x))))
  for(i in 1:nrow(x)){
    hessian_ls[[i]] <- as.numeric(pi[i, ] * (1 - pi[i, ])) * tcrossprod(x[i,], x[i,])
  }

  hessian <- -1 * Reduce("+", hessian_ls) # Sum across all observations
  se_beta <- sqrt(diag(solve(-hessian)))

  solution = theta
  ci_l = solution - qnorm(1 - alpha/2) * se_beta
  ci_u = solution + qnorm(1 - alpha/2) * se_beta
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
              niter = iter,
              converged = (iter < max_iter),
              time_elapsed = time_elapsed$toc - time_elapsed$tic
  ))

}

