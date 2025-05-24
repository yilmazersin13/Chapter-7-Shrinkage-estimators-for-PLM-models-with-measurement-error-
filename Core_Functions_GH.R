# Semiparametric Models with EIV under Right-Censoring: Simulation Study
# Core functions for comparing naive vs deconvolution approaches

# Data generation with measurement error and censoring
generate_data <- function(n, p, sigma_u = 0.1, sigma_eps = 0.5, censoring_rate = 0.2, high_dim = FALSE) {
  T_true <- runif(n, 0, 1)
  U <- rnorm(n, 0, sigma_u)
  W <- T_true + U  # observed covariate with measurement error
  
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- scale(X)
  
  # First 5 coeffs are 1, rest are 0
  beta_true <- c(rep(1, 5), rep(0, p-5))
  
  # True nonparametric function
  m_true <- sin(2 * pi * T_true)
  
  epsilon <- rnorm(n, 0, sigma_eps)
  Z <- X %*% beta_true + m_true + epsilon
  
  # Generate censoring - bit hacky but works
  generate_censoring <- function(Z, censoring_rate) {
    n <- length(Z)
    C <- rep(Inf, n)
    
    if (censoring_rate > 0) {
      n_censor <- floor(n * censoring_rate)
      censor_indices <- sample(1:n, n_censor)
      
      for (i in censor_indices) {
        if (!is.na(Z[i]) && is.finite(Z[i])) {
          min_val <- min(Z[Z > -Inf & Z < Inf], na.rm = TRUE) * 0.5
          if (min_val < Z[i]) {
            C[i] <- runif(1, min = min_val, max = Z[i])
          } else {
            C[i] <- Z[i] * 0.9
          }
        }
      }
    }
    return(C)
  }
  
  C <- generate_censoring(Z, censoring_rate)
  Y <- pmin(Z, C)
  delta <- as.numeric(Z <= C)
  
  cat("Actual censoring rate:", mean(delta == 0), "\n")
  
  return(list(
    Y = Y, X = X, W = W, T_true = T_true,
    delta = delta, beta_true = beta_true,
    m_true = m_true, Z = Z
  ))
}

# Buckley-James for handling censored responses
buckley_james <- function(Y, X, W, delta, max_iter = 10, tol = 1e-5) {
  n <- length(Y)
  Y_star <- Y
  
  for (iter in 1:max_iter) {
    fit <- lm(Y_star ~ X)
    resid <- Y - predict(fit)
    
    km_fit <- survival::survfit(survival::Surv(resid, delta) ~ 1)
    Y_old <- Y_star
    
    # Update censored observations
    for (i in which(delta == 0)) {
      larger_resid_idx <- which(resid >= resid[i])
      if (length(larger_resid_idx) > 0) {
        weights <- km_fit$surv[match(floor(rank(resid[larger_resid_idx])), floor(rank(km_fit$time)))]
        if (any(is.na(weights))) weights[is.na(weights)] <- min(weights, na.rm = TRUE)
        if (length(weights) > 0 && !all(is.na(weights))) {
          weights <- weights / sum(weights)
          Y_star[i] <- sum(Y[larger_resid_idx] * weights, na.rm = TRUE)
        }
      }
    }
    
    if (sqrt(sum((Y_star - Y_old)^2)) < tol) break
  }
  
  return(Y_star)
}

# Deconvolution kernel - the tricky part for handling measurement error
deconvolution_kernel <- function(u, h, sigma_u) {
  n <- length(u)
  K_U <- numeric(n)
  
  for (i in 1:n) {
    t_range <- seq(-5/h, 5/h, length.out = 100)
    integrand_vals <- numeric(length(t_range))
    
    for (j in seq_along(t_range)) {
      t <- t_range[j]
      phi_K <- exp(-t^2/2)  # Fourier transform of Gaussian kernel
      phi_U <- exp(-(sigma_u^2 * (t/h)^2)/2)  # Characteristic function of measurement error
      
      if (abs(phi_U) < 1e-10) {
        integrand_vals[j] <- 0
      } else {
        integrand_vals[j] <- cos(t * u[i]) * phi_K / phi_U
      }
    }
    
    dt <- diff(t_range)[1]
    K_U[i] <- sum(integrand_vals) * dt / (2 * pi)
  }
  
  # Normalize and ensure positivity
  K_U[K_U < 0] <- 0
  if (sum(K_U) > 0) {
    K_U <- K_U / sum(K_U)
  } else {
    # Fallback to naive if deconvolution fails
    K_U <- dnorm(u, 0, h)
    K_U <- K_U / sum(K_U)
  }
  
  return(K_U)
}

# Kernel smoothing with option for naive vs deconvolution
kernel_smoothing <- function(W, Y, h, sigma_u, naive = FALSE) {
  n <- length(W)
  grid <- seq(min(W), max(W), length.out = 100)
  m_hat <- numeric(length(grid))
  
  if (naive) {
    for (j in seq_along(grid)) {
      weights <- dnorm((W - grid[j])/h) / h
      weights <- weights / sum(weights)
      m_hat[j] <- sum(weights * Y)
    }
  } else {
    for (j in seq_along(grid)) {
      K_U <- deconvolution_kernel((W - grid[j])/h, h, sigma_u)
      if (sum(K_U) > 0) {
        weights <- K_U / sum(K_U)
        m_hat[j] <- sum(weights * Y)
      }
    }
  }
  
  # Return prediction function
  predict_func <- function(W_new) {
    m_pred <- numeric(length(W_new))
    for (i in seq_along(W_new)) {
      m_pred[i] <- approx(grid, m_hat, W_new[i], rule = 2)$y
    }
    return(m_pred)
  }
  
  return(list(grid = grid, m_hat = m_hat, predict = predict_func))
}

# Local polynomial alternative to kernel smoothing
local_poly_smooth <- function(W, Y, h, sigma_u, degree = 1, naive = FALSE) {
  n <- length(W)
  grid <- seq(min(W), max(W), length.out = 100)
  m_hat <- numeric(length(grid))
  
  for (j in seq_along(grid)) {
    X_poly <- matrix(NA, n, degree + 1)
    for (d in 0:degree) {
      X_poly[, d + 1] <- (W - grid[j])^d
    }
    
    if (naive) {
      weights <- dnorm((W - grid[j])/h) / h
    } else {
      weights <- deconvolution_kernel((W - grid[j])/h, h, sigma_u)
    }
    
    weights <- diag(weights)
    
    beta_hat <- tryCatch({
      solve(t(X_poly) %*% weights %*% X_poly) %*% t(X_poly) %*% weights %*% Y
    }, error = function(e) {
      rep(NA, degree + 1)
    })
    
    m_hat[j] <- beta_hat[1]
  }
  
  predict_func <- function(W_new) {
    m_pred <- numeric(length(W_new))
    for (i in seq_along(W_new)) {
      m_pred[i] <- approx(grid, m_hat, W_new[i], rule = 2)$y
    }
    return(m_pred)
  }
  
  return(list(grid = grid, m_hat = m_hat, predict = predict_func))
}

# Parametric estimation: Full Model (FM) and Submodel (SM)
estimate_parametric <- function(Y, X, W, Y_star, m_hat_func, method = "FM", 
                                high_dim = FALSE, lambda_seq = NULL) {
  n <- length(Y)
  p <- ncol(X)
  
  partial_resid <- Y_star - m_hat_func(W)
  
  if (method == "FM") {
    if (high_dim) {
      cv_ridge <- glmnet::cv.glmnet(X, partial_resid, alpha = 0)
      lambda_opt <- cv_ridge$lambda.min
      fit <- glmnet::glmnet(X, partial_resid, alpha = 0, lambda = lambda_opt)
      beta_hat <- as.vector(fit$beta)
    } else {
      fit <- lm(partial_resid ~ X - 1)
      beta_hat <- coef(fit)
    }
  } else if (method == "SM") {
    if (high_dim) {
      if (is.null(lambda_seq)) {
        lambda_seq <- 10^seq(-3, 0, length.out = 50)
      }
      cv_lasso <- glmnet::cv.glmnet(X, partial_resid, alpha = 1, lambda = lambda_seq)
      lambda_opt <- cv_lasso$lambda.min
      fit <- glmnet::glmnet(X, partial_resid, alpha = 1, lambda = lambda_opt)
      beta_hat <- as.vector(fit$beta)
    } else {
      # Stepwise AIC selection
      full_model <- lm(partial_resid ~ X - 1)
      step_model <- step(full_model, direction = "both", trace = 0, k = 2)
      beta_hat <- rep(0, p)
      selected_vars <- attr(terms(step_model), "term.labels")
      selected_idx <- as.numeric(gsub("X", "", selected_vars))
      beta_hat[selected_idx] <- coef(step_model)[selected_vars]
    }
  }
  
  return(beta_hat)
}

# Pretest, Shrinkage, and Positive Shrinkage estimators
compute_shrinkage_estimators <- function(beta_FM, beta_SM, Y, X, W, Y_star, m_hat_func) {
  n <- length(Y)
  p <- length(beta_FM)
  p2 <- sum(beta_SM == 0)
  
  resid_FM <- Y_star - X %*% beta_FM - m_hat_func(W)
  resid_SM <- Y_star - X %*% beta_SM - m_hat_func(W)
  
  # Wald-type test statistic
  T_stat <- n * (sum(resid_SM^2) - sum(resid_FM^2)) / sum(resid_FM^2)
  critical_value <- qchisq(0.95, df = p2)
  
  # Pretest: choose FM or SM based on test
  beta_PT <- ifelse(T_stat > critical_value, beta_FM, beta_SM)
  
  # Shrinkage factor
  delta <- 1 - (p2 - 2) / T_stat
  if (is.na(delta) || !is.finite(delta)) delta <- 0
  
  # Shrinkage estimator
  beta_S <- delta * beta_FM + (1 - delta) * beta_SM
  
  # Positive-part shrinkage (James-Stein style)
  delta_plus <- max(0, delta)
  beta_PS <- delta_plus * beta_FM + (1 - delta_plus) * beta_SM
  
  return(list(
    beta_PT = beta_PT,
    beta_S = beta_S,
    beta_PS = beta_PS,
    T_stat = T_stat,
    delta = delta
  ))
}

# Performance evaluation
evaluate_performance <- function(beta_hat, beta_true, m_hat, m_true, T_true, grid) {
  mse_beta <- mean((beta_hat - beta_true)^2)
  m_hat_at_T <- approx(grid, m_hat, T_true, rule = 2)$y
  mise <- mean((m_hat_at_T - m_true)^2)
  
  return(list(mse_beta = mse_beta, mise = mise))
}

# Main simulation runner
run_simulation <- function(n_sim = 100, n = 100, p = 10, h_range = c(0.05, 0.2), 
                           sigma_u = 0.1, sigma_eps = 0.5, censoring_rate = 0.2,
                           high_dim = FALSE, smooth_method = "kernel", degree = 1) {
  
  # Storage structure
  results <- list(
    naive = list(
      FM = list(mse_beta = numeric(n_sim), mise = numeric(n_sim)),
      SM = list(mse_beta = numeric(n_sim), mise = numeric(n_sim)),
      PT = list(mse_beta = numeric(n_sim), mise = numeric(n_sim)),
      S = list(mse_beta = numeric(n_sim), mise = numeric(n_sim)),
      PS = list(mse_beta = numeric(n_sim), mise = numeric(n_sim))
    ),
    conv = list(
      FM = list(mse_beta = numeric(n_sim), mise = numeric(n_sim)),
      SM = list(mse_beta = numeric(n_sim), mise = numeric(n_sim)),
      PT = list(mse_beta = numeric(n_sim), mise = numeric(n_sim)),
      S = list(mse_beta = numeric(n_sim), mise = numeric(n_sim)),
      PS = list(mse_beta = numeric(n_sim), mise = numeric(n_sim))
    )
  )
  
  lambda_seq <- 10^seq(-3, -0.5, length.out = 10)
  
  for (sim in 1:n_sim) {
    print(sim)
    data <- generate_data(n, p, sigma_u, sigma_eps, censoring_rate, high_dim)
    Y_star <- buckley_james(data$Y, data$X, data$W, data$delta)
    
    # Bandwidth selection via CV
    h_values <- seq(h_range[1], h_range[2], length.out = 5)
    cv_errors <- numeric(length(h_values))
    
    for (i in seq_along(h_values)) {
      h <- h_values[i]
      folds <- sample(1:5, n, replace = TRUE)
      cv_error <- 0
      
      for (fold in 1:5) {
        train_idx <- which(folds != fold)
        test_idx <- which(folds == fold)
        
        if (smooth_method == "kernel") {
          m_hat_naive <- kernel_smoothing(data$W[train_idx], Y_star[train_idx], h, sigma_u, naive = TRUE)
          m_hat_conv <- kernel_smoothing(data$W[train_idx], Y_star[train_idx], h, sigma_u, naive = FALSE)
        } else {
          m_hat_naive <- local_poly_smooth(data$W[train_idx], Y_star[train_idx], h, sigma_u, degree, naive = TRUE)
          m_hat_conv <- local_poly_smooth(data$W[train_idx], Y_star[train_idx], h, sigma_u, degree, naive = FALSE)
        }
        
        pred_naive <- m_hat_naive$predict(data$W[test_idx])
        pred_conv <- m_hat_conv$predict(data$W[test_idx])
        cv_error <- cv_error + mean((Y_star[test_idx] - pred_naive)^2 + (Y_star[test_idx] - pred_conv)^2)
      }
      
      cv_errors[i] <- cv_error
    }
    
    h_opt <- h_values[which.min(cv_errors)]
    
    # Fit with optimal bandwidth
    if (smooth_method == "kernel") {
      m_hat_naive <- kernel_smoothing(data$W, Y_star, h_opt, sigma_u, naive = TRUE)
      m_hat_conv <- kernel_smoothing(data$W, Y_star, h_opt, sigma_u, naive = FALSE)
    } else {
      m_hat_naive <- local_poly_smooth(data$W, Y_star, h_opt, sigma_u, degree, naive = TRUE)
      m_hat_conv <- local_poly_smooth(data$W, Y_star, h_opt, sigma_u, degree, naive = FALSE)
    }
    
    # Parametric estimation
    beta_FM_naive <- estimate_parametric(data$Y, data$X, data$W, Y_star, m_hat_naive$predict, "FM", high_dim)
    beta_SM_naive <- estimate_parametric(data$Y, data$X, data$W, Y_star, m_hat_naive$predict, "SM", high_dim, lambda_seq)
    shrink_naive <- compute_shrinkage_estimators(beta_FM_naive, beta_SM_naive, data$Y, data$X, data$W, Y_star, m_hat_naive$predict)
    
    beta_FM_conv <- estimate_parametric(data$Y, data$X, data$W, Y_star, m_hat_conv$predict, "FM", high_dim)
    beta_SM_conv <- estimate_parametric(data$Y, data$X, data$W, Y_star, m_hat_conv$predict, "SM", high_dim, lambda_seq)
    shrink_conv <- compute_shrinkage_estimators(beta_FM_conv, beta_SM_conv, data$Y, data$X, data$W, Y_star, m_hat_conv$predict)
    
    # Store results
    results$naive$FM$mse_beta[sim] <- mean((beta_FM_naive - data$beta_true)^2)
    results$naive$SM$mse_beta[sim] <- mean((beta_SM_naive - data$beta_true)^2)
    results$naive$PT$mse_beta[sim] <- mean((shrink_naive$beta_PT - data$beta_true)^2)
    results$naive$S$mse_beta[sim] <- mean((shrink_naive$beta_S - data$beta_true)^2)
    results$naive$PS$mse_beta[sim] <- mean((shrink_naive$beta_PS - data$beta_true)^2)
    
    m_hat_naive_at_T <- m_hat_naive$predict(data$T_true)
    results$naive$FM$mise[sim] <- mean((m_hat_naive_at_T - data$m_true)^2)
    results$naive$SM$mise[sim] <- results$naive$FM$mise[sim]
    results$naive$PT$mise[sim] <- results$naive$FM$mise[sim]
    results$naive$S$mise[sim] <- results$naive$FM$mise[sim]
    results$naive$PS$mise[sim] <- results$naive$FM$mise[sim]
    
    results$conv$FM$mse_beta[sim] <- mean((beta_FM_conv - data$beta_true)^2)
    results$conv$SM$mse_beta[sim] <- mean((beta_SM_conv - data$beta_true)^2)
    results$conv$PT$mse_beta[sim] <- mean((shrink_conv$beta_PT - data$beta_true)^2)
    results$conv$S$mse_beta[sim] <- mean((shrink_conv$beta_S - data$beta_true)^2)
    results$conv$PS$mse_beta[sim] <- mean((shrink_conv$beta_PS - data$beta_true)^2)
    
    m_hat_conv_at_T <- m_hat_conv$predict(data$T_true)
    results$conv$FM$mise[sim] <- mean((m_hat_conv_at_T - data$m_true)^2)
    results$conv$SM$mise[sim] <- results$conv$FM$mise[sim]
    results$conv$PT$mise[sim] <- results$conv$FM$mise[sim]
    results$conv$S$mise[sim] <- results$conv$FM$mise[sim]
    results$conv$PS$mise[sim] <- results$conv$FM$mise[sim]
    
    if (sim %% 10 == 0) {
      cat("Completed", sim, "of", n_sim, "simulations\n")
    }
  }
  
  # Summarize
  summary <- list(
    naive = list(
      mse_beta = sapply(results$naive, function(x) mean(x$mse_beta)),
      mise = sapply(results$naive, function(x) mean(x$mise))
    ),
    conv = list(
      mse_beta = sapply(results$conv, function(x) mean(x$mse_beta)),
      mise = sapply(results$conv, function(x) mean(x$mise))
    )
  )
  
  return(list(results = results, summary = summary))
}

# Example usage and results formatting
format_results <- function(results, title) {
  cat("\n", title, "\n")
  cat("Parametric component (MSE):\n")
  
  naive_mse <- results$summary$naive$mse_beta
  conv_mse <- results$summary$conv$mse_beta
  
  result_table <- data.frame(
    Estimator = c("FM", "SM", "PT", "S", "PS"),
    Naive = naive_mse,
    Deconvolution = conv_mse,
    Improvement = (naive_mse - conv_mse) / naive_mse * 100
  )
  
  result_table <- result_table[order(result_table$Deconvolution), ]
  print(result_table)
  
  cat("\nNonparametric component (MISE):\n")
  cat("Naive:", results$summary$naive$mise[1], "\n")
  cat("Deconvolution:", results$summary$conv$mise[1], "\n")
  cat("Improvement:", (results$summary$naive$mise[1] - results$summary$conv$mise[1]) / 
        results$summary$naive$mise[1] * 100, "%\n")
}

# Run simulations if sourced
if (FALSE) {  # Set to TRUE to run
  library(survival)
  library(glmnet)
  
  set.seed(123)
  
  # Low-dim case
  low_dim_results <- run_simulation(
    n_sim = 100, n = 50, p = 10,
    sigma_u = 0.1, censoring_rate = 0.2,
    high_dim = FALSE, smooth_method = "kernel"
  )
  
  # High-dim case  
  high_dim_results <- run_simulation(
    n_sim = 100, n = 50, p = 60,
    sigma_u = 0.1, censoring_rate = 0.2,
    high_dim = TRUE, smooth_method = "kernel"
  )
  
  format_results(low_dim_results, "Low-dimensional Results")
  format_results(high_dim_results, "High-dimensional Results")
}