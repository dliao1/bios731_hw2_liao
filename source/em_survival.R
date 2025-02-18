em_censored_data <- function(y, delta, guess, tol = 1e-12, max_iter = 100) {

  lambda <- guess
  n <- length(y)
  tol_criteria <- Inf

  # define vectors to store elements of interest
  observed_ll = rep(NA, length = max_iter)
  lambda_hist = rep(NA, length = max_iter)

  for (iter in 1:max_iter) {

    # E-step: Expectation of censored survival times
    expected_t <- y + (1 - delta) * lambda  # E[t_i | censored]

    # M-step: Update lambda
    lambda_new <- sum(expected_t) / n

    # Use Q function to monitor convergence
    observed_ll[iter] <- -n * log(lambda_new) - sum(expected_t) / lambda_new

    # Stores lambda history
    lambda_hist[iter] <- lambda_new

    if(iter > 1){
      tol_criteria = observed_ll[iter] - observed_ll[iter - 1]

      # Uses squared difference between current and previous iteration of beta
      tol_criteria = (lambda_hist[iter] - lambda_hist[iter - 1])^2
    }

    if (tol_criteria < tol) {
      #cat("Converged in", iter, "iterations.\n")
      break
    }

    lambda <- lambda_new
  }

  return(list(solution = lambda,
              lambda_history = lambda_hist,
              converged = (iter < max_iter),
              niter = iter))
}

get_em_var <- function(y, delta, lambda_estim, guess, B) {
  n <- length(y)
  bootstrap_estimates <- rep(NA, B)


  for (b in 1:B) {
    # Resamples data
    sample_indices <- sample(1:n, size = n, replace = TRUE)
    y_boot <- y[sample_indices]
    delta_boot <- delta[sample_indices]

    # Calculates new parameter estimate for sample data
    lambda_boot <- em_censored_data(y = y_boot,
                                    delta = delta_boot,
                                    guess = guess)


    bootstrap_estimates[b] <- lambda_boot$solution

  }

  # Gets variance of the bootstrapped MLEs
  lambda_var <- var(bootstrap_estimates)

  return(lambda_var)
}
