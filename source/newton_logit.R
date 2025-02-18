newton_logit <- function(beta, x, y, tol = 1e-8, max_iter = 100) {

  beta_cur = beta
  beta_history = gradient_vec = matrix(NA, nrow = max_iter,
                                       ncol = length(beta))
  se_beta <- NA
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

    # Change stopping criterion ?? either converges super fast or not at all..
    if(sqrt(sum(gradient^2)) < tol){
      message("Converged in", iter, "iterations.\n")
      se_beta <- sqrt(diag(solve(-hessian)))
      print(se_beta)
      break
    }

    # Update the solution
    beta_cur = beta_cur - solve(hessian) %*% gradient
  }

  return(list(solution = beta_cur,
              beta_history = beta_history,
              gradient = gradient_vec,
              converged = (iter < max_iter),
              niter = iter,
              se_beta = se_beta))
}
