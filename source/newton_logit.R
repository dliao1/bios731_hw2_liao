newton_logit <- function(beta, x, y, alpha, tol = 1e-8, max_iter = 200) {

  beta_cur = beta
  beta_history = gradient_vec = matrix(NA, nrow = max_iter,
                                       ncol = length(beta))
  se_beta <- NA

  tic()
  for (iter in 1:max_iter) {

    # store results
    beta_history[iter,] = beta_cur

    # Compute the gradient and hessian
    new_pi <- exp(x %*% beta_cur) / (1 + exp(x %*% beta_cur))

    gradient <- t(x) %*% (y - new_pi)
    hessian_ls = as.list(rep(NA, length.out = length(y)))

    for(i in 1:length(y)){
      hessian_ls[[i]] <- as.numeric(new_pi[i, ] * (1 - new_pi[i, ])) * tcrossprod(x[i,], x[i,])
    }

    hessian <- -1 * Reduce("+", hessian_ls)

    gradient_vec[iter,] = gradient

    if(abs(sum(beta_history[iter] - beta_history[iter - 1])) < tol & iter > 1){
      se_beta <- sqrt(diag(solve(-hessian)))
      break
    }

    # Update the solution
    beta_cur = beta_cur - solve(hessian) %*% gradient
  }

  time_elapsed <- toc(quiet = TRUE)

  solution = beta_cur
  ci_l = beta_cur - qnorm(1 - alpha/2) * se_beta
  ci_u = beta_cur + qnorm(1 - alpha/2) * se_beta
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
              beta_history = beta_history,
              gradient = gradient_vec,
              converged = (iter < max_iter),
              niter = iter,
              time_elapsed = time_elapsed$toc - time_elapsed$tic))
}
