# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


#' @title Bayesian Variable importance
#'
#' @description This script contains functions for Bayesian variable importance analysis for GLMM's,
#' utilizing INLA for finding the posteriors and functions for sampling and plotting
#'
#' @author August Arnstad <augustarnstad@gmail.com>
#' @seealso \url{http://r-pkgs.had.co.nz/} for more details on package creation.
#'
#' @importFrom stats scale

# -----------------------


#' Standardize Data
#'
#' This function standardizes the data by subtracting the mean and dividing by the standard deviation.
#'
#' @param data A data frame with numeric columns to be standardized.
#'
#' @return A data frame with standardized values.
#'
#' @examples
#' data <- data.frame(a = rnorm(100), b = rnorm(100))
#' standardized_data <- standardize_data(data)
#'
standardize_data <- function(data) {
  # Exclude non-numeric columns (if any)
  numeric_data <- data[, sapply(data, is.numeric)]

  # Compute means and standard deviations
  means <- colMeans(numeric_data, na.rm = TRUE)
  sds <- apply(numeric_data, 2, sd, na.rm = TRUE)

  # Standardize each feature
  standardized_data <- sweep(numeric_data, 2, means, "-")
  standardized_data <- sweep(standardized_data, 2, sds, "/")

  # Return the standardized data
  return(standardized_data)
}


#' Z-transform for Data
#'
#' This function applies an SVD decomposition to the fixed effects in data.
#'
#' @param X A matrix or data frame to be transformed.
#'
#' @return A list containing the transformed data, R matrix, and lambda matrix.
#'
#' @examples
#' data <- data.frame(a = rnorm(100), b = rnorm(100))
#' transformed_data <- SVD_decomp(data)
#'
SVD_decomp <- function(X){

  original_colnames <- colnames(X)
  X <- scale(X)

  # Calculate eigenvalues and eigenvectors
  e <- eigen(t(X) %*% X)
  Q <- e$vectors
  D <- diag(sqrt(e$values), nrow = nrow(Q))
  Dinv <- diag(1/sqrt(e$values), nrow = nrow(Q))

  # Calculate R_xx^(-0.5)
  R <- sqrt(nrow(X) - 1) * Q %*% Dinv %*% t(Q)

  # Calculate the transformed numerical fixed effects
  Z <- X %*% R
  colnames(Z) <- original_colnames <- colnames(X)


  lambda <- 1 / sqrt(nrow(X) - 1) * Q %*% D %*% t(Q)
  return(list(Z = Z, R = R, lambda = lambda))
}


#' Preprocess Data Based on Formula
#'
#' This function preprocesses the data based on a given formula. It standardizes the data,
#' extracts formula components, checks variable existence, identifies random effects,
#' and constructs the design matrix.
#'
#' @param data A data frame containing the data to be processed.
#' @param formula A formula specifying the model.
#'
#' @return A list containing the response vector, the design matrix, and the random effects.
#'
#' @examples
#' data <- data.frame(Y = rnorm(100), V2 = rnorm(100), V3 = rnorm(100), gamma = rnorm(100))
#' result <- preprocess_data_formula(data, Y ~ V2 + V3 + (1|gamma))
#'
preprocess_data_formula <- function(formula, data) {
  #Standardize data
  data <- standardize_data(data)
  # Extract formula components
  response <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]

  # Check if all variables in formula exist in data
  all_vars <- c(response, predictors)
  if (!all(all_vars %in% colnames(data))) {
    stop("Some variables in the formula are not found in the data.")
  }

  # Identify random effects using terms object
  term_objects <- terms(formula)
  random_effects <- gsub(".*\\|", "", attr(term_objects, "term.labels")[grep("\\|", attr(term_objects, "term.labels"))])
  random_effects <- trimws(random_effects)


  # Remove random effects terms from predictors
  predictors <- setdiff(predictors, random_effects)

  Y <- data[[response]]

  # Build the design matrix for fixed effects
  X <- model.matrix(~ . - 1, data = data[, predictors, drop = FALSE])


  # Extract random effects if they exist
  if (length(random_effects) > 0) {
    gamma <- lapply(random_effects, function(re) data[[re]])
    names(gamma) <- random_effects
  } else {
    gamma <- NULL
  }

  # Return the processed data
  list(Y = Y, X = X, random = gamma, predictor_names = predictors)
}


#' Run INLA Model
#'
#' This function processes the data based on the given formula and then applies the INLA model.
#'
#' @param formula A formula specifying the model.
#' @param data A preprocessed list containing the response vector, design matrix, and random effects.
#'
#' @return The results of the INLA model.
#'
run_INLA <- function(formula, data) {

  # Extract Z from the input list
  Z <- data$Z

  # Define the hyperparameter settings directly for the formula
  hyper_formula <- function(effect_name) {
    paste0("f(", effect_name, ", model = 'iid', hyper = list(prec = list(initial = log(0.8), prior = 'pc.prec', param = c(1,0.02))))")
  }

  # As it stands now, I place the same priors on all random effects, as I found it difficult to be able to give the user this functionality.
  # This must be fixed, kinda ugly now.
  # Also a bit problematic because here we take in the Z values, not the values transformed back. What should I do? I can transform back in
  # this step, but I fear that it is a bit too early. I can also manually change back later. For now, it is not done, as for this dataset
  # the matrix is very close to the identity.

  # Check if random effects exist in the data list
  if (!is.null(data$random)) {
    # Extract random effect names
    random_names <- names(data$random)

    # Construct the random effects parts of the formula
    random_effects_formulas <- sapply(random_names, hyper_formula)

    # Construct INLA formula with random effects
    inla_formula <- as.formula(paste("Y ~ Z +", paste(random_effects_formulas, collapse = " + ")))

  } else {
    # Construct INLA formula without random effects
    inla_formula <- as.formula("Y ~ Z")
  }

  inla_data <- list(Y = data$Y, Z = Z)
  # Add random effects to data list if they exist
  if (!is.null(data$random)) {
    inla_data <- c(inla_data, data$random)
  }

  # Run the INLA model
  model_Z <- inla(formula = inla_formula,
                  family = "gaussian",
                  data = inla_data,
                  control.compute = list(dic = F,
                                         return.marginals = TRUE,
                                         config = TRUE,
                                         return.marginals.predictor=TRUE),
                  control.family = list(hyper = list(theta = list(initial = log(0.5),
                                                                  prior = "pc.prec",
                                                                  param = c(sqrt(2),0.05)))),
                  control.predictor = list(compute = TRUE)
  )
  return(model_Z)
}


#' Bayesian Importance Analysis
#'
#' A wrapper function that handles data standardization, preprocessing, INLA modeling,
#' and potentially sampling from the posterior distribution.
#'
#' @param data A data frame to be analyzed.
#' @param formula A formula specifying the model.
#' @param ... Additional arguments.
#'
#' @return Returns the results of the INLA model.
#'
#' @examples
#' data <- data.frame(Y = rnorm(100), V1 = rnorm(100), V2 = rnorm(100))
#' results <- run_bayesian_imp(Y ~ V1 + V2, data, )
#'
run_bayesian_imp <- function(formula, data, ...) {
  require(INLA)

  # Standardize the data
  standardized_data <- standardize_data(data)

  # Preprocess the data based on the given formula
  preprocessed_data <- preprocess_data_formula(formula, standardized_data)

  # Apply SVD decomposition to the fixed effects design matrix
  svd_results <- SVD_decomp(preprocessed_data$X)

  # Add the results from SVD decomposition to the preprocessed data
  preprocessed_data$Z <- svd_results$Z
  preprocessed_data$R <- svd_results$R
  preprocessed_data$lambda <- svd_results$lambda

  # bayesian_imp_env$context <- list()
  # bayesian_imp_env$context$lambda_matrix <- svd_results$lambda
  # bayesian_imp_env$context$data_raw <- data
  # bayesian_imp_env$context$formula <- formula

  #Store the number of classes for all random effects.
  term_objects <- terms(formula)
  random_effects <- gsub(".*\\|", "", attr(term_objects, "term.labels")[grep("\\|", attr(term_objects, "term.labels"))])
  random_effects <- trimws(random_effects)

  # total_classes <- 0
  #
  # if (length(random_effects) != 0){
  #   for (random_effect in random_effects){
  #     num_unique <- length(unique(data[[random_effect]]))
  #     bayesian_imp_env$context$num_classes[[random_effect]] <- num_unique
  #
  #     # Update the total count
  #     total_classes <- total_classes + num_unique
  #   }
  # }
  # bayesian_imp_env$context$total_num_classes <- total_classes


  # Run the INLA model
  model_results <- run_INLA(formula, preprocessed_data)

  return(model_results)

  # Check if the user wants the model as output
  # if (!return_samples) {
  #   return(model_results)
  # }
  # else{
  #   # Sample from the INLA model
  #   model_samples <- sample_inla_model(model_results, n_samp)
  #
  #   # Return the samples
  #   return(list(samples = model_samples, model = model_results, lambda = bayesian_imp_env$context$lambda_matrix))
  # }
}


#' Plot Posterior Distributions
#'
#' A function to visualize the posterior distributions of fixed and random effects from the Bayesian importance analysis.
#'
#' NOTE: This plots the approximated marginals of the model parameters. If there is correlation among the predictors, consider using sample_posteriors instead.
#'
#' @param model_results The results from the INLA model, as obtained from \code{run_bayesian_imp}.
#' @param ... Additional arguments to customize the plot.
#'
#' @return A ggplot object visualizing the posterior distributions/importances. Initialize a variable to hold the plot and call the variable.
#'
#' @details
#' This function uses ggplot2 to visualize the posterior distributions of the fixed effects and random effects
#' obtained from the Bayesian importance analysis. It provides insights into the distribution of each parameter
#' and assists in understanding their significance.
#'
#' @examples
#' data <- data.frame(Y = rnorm(100), V1 = rnorm(100), V2 = rnorm(100))
#' model_results <- run_bayesian_imp(data, Y ~ V1 + V2)
#' plot_model <- plot_posteriors(model_results)
#' plot_model

plot_posteriors <- function(model, importance=FALSE, modelname="model") {
  # Get the marginals
  variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))
  fixed_marginals_list <- lapply(model$marginals.fixed, function(x) x)

  # Extract names
  random_effect_names <- names(model$marginals.hyperpar)
  fixed_effect_names <- names(model$marginals.fixed)
  random_effect_names <- gsub("Precision for ", "Variance of ", random_effect_names)

  # Get posterior means for random effects and fixed effects
  random_means <- 1/sapply(model$marginals.hyperpar, function(x) inla.zmarginal(x, silent=TRUE)$mean)
  if (!importance){
    fixed_means <- sapply(model$marginals.fixed, function(x) inla.zmarginal(x, silent=TRUE)$mean)
  }else{
    fixed_means <- sapply(model$marginals.fixed, function(x) inla.zmarginal(x, silent=TRUE)$mean)^2
  }
  #fixed_means <- sapply(model$marginals.fixed, function(x) inla.zmarginal(x)$mean)


  # Rename with mean values for the legend
  random_effect_names <- paste(random_effect_names, " (Posterior Mean:", round(random_means, 3), ")")
  fixed_effect_names <- paste(fixed_effect_names, " (Posterior Mean:", round(fixed_means, 3), ")")

  # random_effect_names <- c(expression(sigma[epsilon]^2 ~ " (Posterior Mean:" ~ round(random_means, 3)[1] ~ ")"),
  #                          expression(sigma[alpha]^2 ~ " (Posterior Mean:" ~ round(random_means, 3)[2] ~ ")"))
  #
  # fixed_effect_names <- c(expression(beta[1]^2 ~ " (Posterior Mean:" ~ round(fixed_means, 3)[2] ~ ")"),
  #                         expression(beta[2]^2 ~ " (Posterior Mean:" ~ round(fixed_means, 3)[3] ~ ")"),
  #                         expression(beta[3]^2 ~ " (Posterior Mean:" ~ round(fixed_means, 3)[4] ~ ")"))


  # Create data frames
  df_list <- lapply(1:length(variance_marginals_list), function(i) {
    data.frame(
      x = variance_marginals_list[[i]][, 1],
      y = variance_marginals_list[[i]][, 2],
      Effect = random_effect_names[i]
    )
  })

  if (!importance){
    df_list_fixed <- lapply(1:length(fixed_marginals_list), function(i) {
      data.frame(
        x = fixed_marginals_list[[i]][, 1],
        y = fixed_marginals_list[[i]][, 2],
        Effect = fixed_effect_names[i]
      )
    })
  }else{
    df_list_fixed <- lapply(1:length(fixed_marginals_list), function(i) {
      data.frame(
        x = fixed_marginals_list[[i]][, 1]^2,
        y = fixed_marginals_list[[i]][, 2],
        Effect = fixed_effect_names[i]
      )
    })
  }

  # Combine data frames
  df_combined <- do.call(rbind, c(df_list, df_list_fixed))

  if (!importance){
    plot_title = paste("Posterior distributions")
  }else{
    plot_title = paste("Posterior proportion of variance")
  }

  # Plot using ggplot
  p <- ggplot(df_combined, aes(x = x, y = y, color = Effect)) +
    geom_line() +
    labs(title = plot_title, x = "Values", y = "Density") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.title = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          strip.text = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    scale_color_manual(values = rainbow(length(df_list) + length(df_list_fixed)))
    #geom_vline(aes(xintercept = 1/6), linetype = "dashed", color = "black") +    #This is true importance for beta_1, gamma, eta, epsilon
    #geom_vline(aes(xintercept = 2/6), linetype = "dashed", color = "black")      #This is true importance for beta_2
  return(list(posterior_marginals = df_combined, posterior_plot = p))
}

#' Sample Posterior Distributions for Bayesian Models
#'
#' This function samples posterior distributions from a Bayesian model, calculates
#' the importance of each variable, and provides variance marginals for the random effects.
#' It uses a provided Bayesian model formula and data to generate these samples.
#'
#' @param formula An object of class `formula` (or one that can be coerced to that class):
#'                a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param n_samp The number of samples to be drawn from the posterior distribution.
#' @param n The number of observations in the data.
#' @param n_classes The number of total classes for all random effects in the model.
#'
#' @return A list containing three elements:
#'         \itemize{
#'           \item \code{beta}: A matrix of sampled beta coefficients for the fixed effects.
#'           \item \code{importance}: A matrix of importance values for each variable.
#'           \item \code{marginals}: A list of data frames, each containing the x and y
#'                                  values for the variance marginals of the random effects.
#'         }
#'
#' @examples
#' # Define the model formula
#' formula <- Y ~ V2 + V3 + V4 + (1 | gamma)
#'
#' # Specify the data frame and other parameters
#' data_bayes <- my_data_frame # Replace with your data frame
#' n_samp <- 100
#' n <- nrow(data_bayes)
#' n_classes <- 3 # Set as appropriate
#'
#' # Run the function
#' results <- sample_posteriors(formula, data_bayes, n_samp, n, n_classes)
#'
#' #' @export
sample_posteriors <- function(formula, data, n_samp, n, n_classes){

  model <- run_bayesian_imp(formula, data=data)

  formula_str = deparse(formula)
  split_formula <- strsplit(formula_str, "~")[[1]][2]
  fixed_effects_str <- ifelse(grepl("\\(", split_formula), strsplit(split_formula, "\\(")[[1]][1], split_formula)
  fixed_effects_str <- trimws(fixed_effects_str)
  fixed_effects_str <- gsub("\\+\\s*$", "", fixed_effects_str)
  fixed_effect_vars <- unlist(strsplit(fixed_effects_str, "\\s*\\+\\s*"))
  fixed_effects <- trimws(fixed_effect_vars)

  # Extracting random effects
  random_effects_str <- regmatches(split_formula, gregexpr("\\(.*?\\)", split_formula))
  random_effects <- lapply(random_effects_str, function(x) gsub("[\\(\\)]", "", x))
  random_effects <- unlist(random_effects)
  random_effects <- strsplit(random_effects, "\\|")

  # Extracting random effect variable names, excluding '1'
  random_effect_vars <- unique(unlist(lapply(random_effects, function(x) {
    vars <- trimws(x)
    # Taking the second part after '|', assuming the first part is '1' or grouping factor
    if (length(vars) > 1) {
      return(vars[2])
    } else {
      return(vars[1])
    }
  })))
  random_effect_vars

  random_effects_data <- data[, random_effect_vars, drop = FALSE]

  ncols_ran <- length(random_effects_data[1, ])
  nclasses_vec <- matrix(NA, nrow=ncols_ran, ncol=1)
  for (i in 1:ncols_ran){
    nclasses_vec[i] <- length(unique(random_effects_data[, i]))
  }

  n_classes=sum(nclasses_vec)

  X = data[fixed_effects]

  beta_mat <- matrix(NA, nrow=n_samp, ncol=ncol(X))
  importance_mat <- matrix(NA, nrow=n_samp, ncol=ncol(X))
  R2_mat <- matrix(NA, nrow=n_samp, ncol=1)
  R2_cond_mat <- matrix(NA, nrow=n_samp, ncol=1)

  SVD = BayesianImportance::SVD_decomp(X)

  lambda = SVD$lambda

  samps_Z <- inla.posterior.sample(model, n = n_samp)

  num_fixed <- ncol(X)

  for (i in 1:n_samp){
    beta <- samps_Z[[i]]$latent[(n + n_classes + 2):(n + n_classes + num_fixed + 1)]
    #sigma_sq = var(samps_Z[[i]]$latent[n:(n+n_classes)]) + as.numeric(1/samps_Z[[i]]$hyperpar['Precision for the Gaussian observations']) #sum(1/samps_Z[[i]]$hyperpar)
    beta_mat[i, ] <- beta
    importance_mat[i, ] <- lambda^2 %*% beta^2

    var_sum <- 0
    start_idx <- n + 1

    # Loop through each class and sum their variances
    for (j in 1:length(nclasses_vec)){
      end_idx <- start_idx + nclasses_vec[j] - 1
      var_sum <- var_sum + var(samps_Z[[i]]$latent[start_idx:end_idx])
      start_idx <- end_idx + 1
    }

    # Add the precision for Gaussian observations
    sigma_sq <- var_sum + as.numeric(1/samps_Z[[i]]$hyperpar['Precision for the Gaussian observations'])
    #print(sigma_sq)

    R2_mat[i, ] <- sum(importance_mat[i, ])/(sum(importance_mat[i, ]) + sigma_sq)
    R2_cond_mat[i, ] <- (sum(importance_mat[i, ])+var_sum)/(sum(importance_mat[i, ]) + sigma_sq)
  }

  variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))
  random_effect_names <- names(model$marginals.hyperpar)


  df_list <- lapply(1:length(variance_marginals_list), function(i) {
    data.frame(
      x = variance_marginals_list[[i]][, 1],
      y = variance_marginals_list[[i]][, 2],
      Effect = random_effect_names[i]
    )
  })

  return(list(beta = beta_mat, importance=importance_mat, marginals = df_list, r2=R2_mat, r2_cond=R2_cond_mat))
}
#' sample_posteriors <- function(formula, data, n_samp, n, n_classes){
#'
#'   model <- run_bayesian_imp(formula, data=data)
#'
#'   formula_str = deparse(formula)
#'   split_formula <- strsplit(formula_str, "~")[[1]][2]
#'   fixed_effects_str <- ifelse(grepl("\\(", split_formula), strsplit(split_formula, "\\(")[[1]][1], split_formula)
#'   fixed_effects_str <- trimws(fixed_effects_str)
#'   fixed_effects_str <- gsub("\\+\\s*$", "", fixed_effects_str)
#'   fixed_effect_vars <- unlist(strsplit(fixed_effects_str, "\\s*\\+\\s*"))
#'   fixed_effects <- trimws(fixed_effect_vars)
#'
#'   # Extracting random effects
#'   random_effects_str <- regmatches(split_formula, gregexpr("\\(.*?\\)", split_formula))
#'   random_effects <- lapply(random_effects_str, function(x) gsub("[\\(\\)]", "", x))
#'   random_effects <- unlist(random_effects)
#'   random_effects <- strsplit(random_effects, "\\|")
#'
#'   # Extracting random effect variable names, excluding '1'
#'   random_effect_vars <- unique(unlist(lapply(random_effects, function(x) {
#'     vars <- trimws(x)
#'     # Taking the second part after '|', assuming the first part is '1' or grouping factor
#'     if (length(vars) > 1) {
#'       return(vars[2])
#'     } else {
#'       return(vars[1])
#'     }
#'   })))
#'   random_effect_vars
#'
#'   random_effects_data <- data[, random_effect_vars, drop = FALSE]
#'
#'   ncols_ran <- length(random_effects_data[1, ])
#'   nclasses_vec <- matrix(NA, nrow=ncols_ran, ncol=1)
#'   for (i in 1:ncols_ran){
#'     nclasses_vec[i] <- length(unique(random_effects_data[, i]))
#'   }
#'
#'   n_classes=sum(nclasses_vec)
#'
#'   X = data[fixed_effects]
#'
#'   beta_mat <- matrix(NA, nrow=n_samp, ncol=ncol(X))
#'   importance_mat <- matrix(NA, nrow=n_samp, ncol=ncol(X))
#'   R2_mat <- matrix(NA, nrow=n_samp, ncol=1)
#'
#'   SVD = BayesianImportance::SVD_decomp(X)
#'
#'   lambda = SVD$lambda
#'
#'   samps_Z <- inla.posterior.sample(model, n = n_samp)
#'
#'   num_fixed <- ncol(X)
#'
#'   start_idx <- n + 1
#'
#'   for (i in 1:n_samp){
#'     beta <- samps_Z[[i]]$latent[(n + n_classes + 2):(n + n_classes + num_fixed + 1)]
#'     #sigma_sq = var(samps_Z[[i]]$latent[n:(n+n_classes)]) + as.numeric(1/samps_Z[[i]]$hyperpar['Precision for the Gaussian observations']) #sum(1/samps_Z[[i]]$hyperpar)
#'     beta_mat[i, ] <- beta
#'     importance_mat[i, ] <- lambda^2 %*% beta^2
#'
#'     var_sum <- 0
#'
#'     # Loop through each class and sum their variances
#'     for (j in 1:length(nclasses)){
#'       end_idx <- start_idx + nclasses[j] - 1
#'       var_sum <- var_sum + var(samps_Z[[i]]$latent[start_idx:end_idx])
#'       start_idx <- end_idx + 1
#'     }
#'
#'     # Add the precision for Gaussian observations
#'     sigma_sq <- var_sum + as.numeric(1/samps_Z[[i]]$hyperpar['Precision for the Gaussian observations'])
#'
#'
#'     R2_mat[i, ] <- sum(importance_mat[i, ])/(sum(importance_mat[i, ]) + sigma_sq)
#'   }
#'
#'   variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))
#'   random_effect_names <- names(model$marginals.hyperpar)
#'
#'
#'   df_list <- lapply(1:length(variance_marginals_list), function(i) {
#'     data.frame(
#'       x = variance_marginals_list[[i]][, 1],
#'       y = variance_marginals_list[[i]][, 2],
#'       Effect = random_effect_names[i]
#'     )
#'   })
#'
#'   return(list(beta = beta_mat, importance=importance_mat, marginals = df_list, r2=R2_mat))
#' }
#'

#' Sample Response Matrix from INLA Model
#'
#' This function samples a response matrix from the given INLA model.
#'
#' @param model An object of class `inla`, typically the result of fitting a model with the `inla` function.
#' @param n The number of predictors.
#' @param s The number of samples to draw for each predictor.
#'
#' @return A matrix where each row corresponds to a predictor and each column corresponds to a sampled value.
#' @export
sample_response_matrix <- function(model, n, s) {
  # Pre-allocate a matrix to store the results
  result_matrix <- matrix(0, nrow=n, ncol=s)

  for (i in 1:n) {
    max_digits <- nchar(n)
    zero_padding <- max(0, max_digits - nchar(i))
    predictor_name <- paste0("fitted.Predictor.", paste(rep("0", zero_padding), collapse = ""), i)

    # Now draw s samples at once
    samples <- inla.rmarginal(s, marginal = model$marginals.fitted.values[[predictor_name]])
    result_matrix[i, ] <- samples
  }

  return(result_matrix)
}


#' Sample Response Matrix from INLA Model
#'
#' This function samples a response matrix from the given INLA model.
#'
#' @param model An object of class `inla`, typically the result of fitting a model with the `inla` function.
#' @param n The number of predictors.
#' @param s The number of samples to draw for each predictor.
#'
#' @return A matrix where each row corresponds to a predictor and each column corresponds to a sampled value.
#' @export
gelman_r2_metrics <- function(model, s = 1000, plot = FALSE, modelname="model") {

  # Number of predictors
  n <- length(model$marginals.fitted.values)

  y_samples <- sample_response_matrix(model, n, s)

  y_var <- apply(y_samples, MARGIN=2, FUN=var)

  epsilon_var = 1/inla.rmarginal(s, marginal = model$marginals.hyperpar$"Precision for the Gaussian observations")

  gelman_r2_cond = y_var/(y_var + epsilon_var)


  # Plotting
  if (plot) {
    library(ggplot2)

    p <- ggplot() +
      #geom_histogram(aes(x = gelman_r2_cond), bins = 30, fill = "skyblue") +
      geom_density(aes(x = gelman_r2_cond), fill = "skyblue", alpha = 0.7) +
      labs(title = paste("Distribution of Gelman Conditional R² of", modelname), x = "Gelman R²", y = "Frequency") +
      theme_minimal()

    return(list(conditional_gelman_r2 = gelman_r2_cond, plot = p))
  }

  # Return the R² values
  return(gelman_r2_cond)
}








#' #' Sample from INLA Model
#' #'
#' #' This function draws samples from the posterior distribution of an INLA model.
#' #'
#' #' @param model_Z The INLA model.
#' #' @param n_samp Number of samples to draw.
#' #'
#' #' @return The drawn samples.
#' #'
#' sample_inla_model <- function(model_Z, n_samp) {
#'   # Sample from the posterior distribution of the INLA model
#'   samps_Z <- inla.posterior.sample(model_Z, n = n_samp)
#'
#'   return(samps_Z)
#' }
#'
#' #' Summarize Importance of Predictors
#' #'
#' #' This function provides a summary of the importance of predictors in the model.
#' #' It presents the estimates, importance, and standard errors for fixed effects,
#' #' as well as the variance and standard errors for random effects.
#' #'
#' #' @param model The results of the INLA model.
#' #'
#' #' @return Prints the summary table.
#' #'
#' summary_importances <- function(model){
#'
#'   lambda_matrix <- bayesian_imp_env$context$lambda_matrix
#'
#'   # Extract fixed effects coefficients, squared coefficients, and standard errors
#'   beta <- model$summary.fixed[, "mean"]
#'   importances <- c(beta[1]^2, lambda_matrix^2 %*% beta[2:length(beta)]^2)
#'   beta_se <- model$summary.fixed[, "sd"]
#'
#'   # Extract covariate names
#'   covariate_names <- rownames(model$summary.fixed)
#'
#'   # Extract hyperparameters (precisions for random effects) and compute variances and standard errors
#'   hyperparameters <- model$summary.hyperpar[, "mean"]
#'   variance_randoms <- 1/hyperparameters
#'   random_se <- model$summary.hyperpar[, "sd"]
#'
#'   # Construct the fixed effect portion
#'   fixed_effects <- data.frame(
#'     Predictor = covariate_names,
#'     Estimate = beta,
#'     Importance = importances,
#'     Std_Error = beta_se
#'   )
#'   fixed_effects_output <- capture.output(print(fixed_effects, digits = 4, row.names = FALSE))
#'
#'   # Construct the random effect portion (if available)
#'   if(length(hyperparameters) > 1) {
#'     random_effects <- data.frame(
#'       Random_Effect = c("Gaussian Observations", names(model$summary.random)),
#'       Precision = hyperparameters,
#'       Variance = variance_randoms,
#'       Std_Error = random_se
#'     )
#'     random_effects_output <- capture.output(print(random_effects, digits = 4, row.names = FALSE))
#'     total_sum_explained = sum(importances) + sum(variance_randoms[2:length(hyperparameters)])
#'   } else {
#'     random_effects <- data.frame(
#'       Random_Effect = "Gaussian Observations",
#'       Precision = hyperparameters,
#'       Variance = variance_randoms,
#'       Std_Error = random_se
#'     )
#'     random_effects_output <- capture.output(print(random_effects, digits = 4, row.names = FALSE))
#'     total_sum_explained = sum(importances)
#'   }
#'
#'   # Calculate the sum of variances for fixed and random effects
#'   #total_sum = sum(squared_beta)
#'
#'   # Construct the bottom portion with total variance and other potential relevant quantities
#'   bottom_output <- sprintf("Total variance captured by the model: %.4f", total_sum_explained)
#'
#'   # Combine the outputs
#'   final_output <- c(
#'     "Fixed Effects:",
#'     fixed_effects_output,
#'     "\nRandom Effects:",
#'     random_effects_output,
#'     "\n",
#'     bottom_output
#'   )
#'
#'   cat(final_output, sep = "\n")
#' }
#'
#'
#' #' Extract Information from Samples
#' #'
#' #' Extracts relevant information from the posterior samples, such as beta coefficients,
#' #' random effect variances, and Gelman's R^2.
#' #'
#' #' @param samps_Z Posterior samples.
#' #' @param model The INLA model.
#' #'
#' #' @return A list containing beta coefficients, random effect variances, Gelman's R^2,
#' #'         and other related information.
#' #'
#' extract_samples_info <- function(samps_Z, model) {
#'
#'   # I wish to lower the number of input parameters in this function, as it is surely redundant to take in this much.
#'   # The n_class term is a bit hard to handle to be honest. Not sure what the best way to store it from the beginning is.
#'   # This is a very brute force method as is
#'
#'   data <- bayesian_imp_env$context$data_raw
#'   formula <- bayesian_imp_env$context$formula
#'   n_class <- bayesian_imp_env$context$total_num_classes
#'
#'   # Parse the formula to extract fixed effects and random effects
#'   terms <- all.vars(formula)
#'   fixed_effects <- terms[!grepl("\\|", terms)] # exclude terms with '|'
#'   random_effects <- unique(gsub(".*\\((.*)\\|.*", "\\1", terms[grep("\\|", terms)])) # extract the random effect names
#'
#'   n_fixed_effects <- length(rownames(model$summary.fixed))
#'
#'   n_random_effects <- length(names(model$summary.random))
#'   n <- nrow(data)
#'
#'   # Beta coefficients
#'   beta <- matrix(NA, nrow=length(samps_Z), ncol=n_fixed_effects)
#'   colnames(beta) <- rownames(model$summary.fixed)
#'
#'   if (n_random_effects != 0){
#'     # Random effect variances
#'     random_var <- matrix(NA, nrow=length(samps_Z), ncol=n_random_effects)
#'     colnames(random_var) <- names(model$summary.random)
#'
#'     # Proportion of total variance due to each random effect
#'     random_prop <- matrix(NA, nrow=length(samps_Z), ncol=n_random_effects)
#'     colnames(random_prop) <- paste0(names(model$summary.random), "_prop")
#'   }
#'   else{
#'     # Random effect variances
#'     random_var <- NULL
#'     colnames(random_var) <- NULL
#'
#'     # Proportion of total variance due to each random effect
#'     random_prop <- NULL
#'     colnames(random_prop) <- NULL
#'   }
#'
#'   # Samples of the response variable
#'   Y_samps <- matrix(NA, nrow=length(samps_Z), ncol=n)
#'
#'   # Gelman R2 and sigma squared
#'   gelman_R2 <- numeric(length(samps_Z))
#'   sigma_sq <- numeric(length(samps_Z))
#'
#'   for (i in 1:length(samps_Z)) {
#'
#'     # Extracting beta coefficients
#'     for (j in 1:n_fixed_effects) {
#'       beta[i, j] <- samps_Z[[i]]$latent[(n + n_class + j)]
#'     }
#'
#'     # Transform back from SVD decomposition
#'
#'     # Extracting random effect variances and sigma squared
#'     sigma_sq[i] <- as.numeric(1/samps_Z[[i]]$hyperpar[[1]])
#'     if (n_random_effects != 0){
#'       for (j in 1:n_random_effects) {
#'         random_var[i, j] <- as.numeric(1/samps_Z[[i]]$hyperpar[[j + 1]])
#'       }
#'     }
#'     else{
#'       random_var <- NULL
#'     }
#'
#'     # Calculate Y_samps using the formula
#'     linear_predictor <- rep(beta[i, 1], n)
#'     for (j in 2:n_fixed_effects) {
#'       linear_predictor <- linear_predictor + beta[i, j] * data[[fixed_effects[j]]]
#'     }
#'
#'     total_Wu <- 0
#'     if (!is.null(random_var)){
#'       for (j in 1:n_random_effects) {
#'         Wu_j <- rnorm(n, mean=0, sd=sqrt(random_var[i, j]))
#'         total_Wu <- total_Wu + Wu_j
#'       }
#'     }
#'
#'     Y_samps[i, ] = linear_predictor + total_Wu
#'
#'     if (!is.null(random_var)) {
#'       total_random_var <- sum(random_var[i, ])
#'       total_variance <- var(Y_samps[i, ]) + sigma_sq[i] + total_random_var
#'
#'       for (j in 1:n_random_effects) {
#'         random_prop[i, j] <- random_var[i, j] / total_variance
#'       }
#'       gelman_R2[i] <- (total_random_var + var(Y_samps[i, ])) / total_variance
#'     } else {
#'       gelman_R2[i] <- var(Y_samps[i, ]) / (var(Y_samps[i, ]) + sigma_sq[i])
#'     }
#'   }
#'
#'   # Return all extracted values as a list
#'   return(list(beta = beta, random_var = random_var, random_prop = random_prop,
#'               Y_samps = Y_samps, gelman_R2 = gelman_R2, sigma_sq = sigma_sq, nclass = n_class))
#' }
#'
#'
#' #' Plot Summary Information of Samples
#' #'
#' #' This function creates three plots:
#' #' 1. A histogram of variable importance for each predictor.
#' #' 2. A histogram with density approximation for R2.
#' #' 3. A histogram with density approximation for Gelman R2.
#' #'
#' #' @param samples_matrix A list containing the posterior samples.
#' #'        It should include `beta`, `random_var`, and `gelman_R2` matrices.
#' #'
#' #' @return A list containing three ggplot2 plots:
#' #' `VariableImportancePlot`, `RDistributionPlot`, and `GelmanR2DistributionPlot`.
#' #'
#' #' @examples
#' #' # Assuming `samples` is a list containing the necessary matrices
#' #' plots <- plot_summary_info(samples)
#' #' print(plots$VariableImportancePlot)
#' #' print(plots$RDistributionPlot)
#' #' print(plots$GelmanR2DistributionPlot)
#' #'
#' #' @importFrom ggplot2 ggplot geom_histogram geom_density aes theme_minimal labs theme
#' #' @importFrom RColorBrewer scale_fill_brewer
#' #'
#' plot_summary_info <- function(samples_matrix) {
#'   require(ggplot2)
#'   require(RColorBrewer)
#'
#'   lambda_matrix <- bayesian_imp_env$context$lambda_matrix
#'
#'   all_beta_values_matrix <- matrix(0, nrow(samples_matrix$beta), ncol(samples_matrix$beta) - 1)
#'   # Go through each row of samples_matrix$beta and multiply with lambda_matrix
#'   for (i in 1:nrow(samples_matrix$beta)) {
#'     beta_row <- as.matrix(samples_matrix$beta[i, -1]^2)
#'     all_beta_values_matrix[i, ] <- as.vector(lambda_matrix^2 %*% beta_row)
#'   }
#'
#'   all_beta_values <- as.vector(all_beta_values_matrix)
#'
#'   all_beta_labels <- rep(paste0("Beta ", 2:(ncol(samples_matrix$beta))), each=nrow(samples_matrix$beta))
#'
#'   # Check if random effects are present
#'   if (!is.null(samples_matrix$random_var)) {
#'     all_random_values <- as.vector(samples_matrix$random_var)
#'     all_random_labels <- rep(paste0("sigma ", 1:ncol(samples_matrix$random_var)), each=nrow(samples_matrix$random_var))
#'   } else {
#'     all_random_values <- numeric(0)
#'     all_random_labels <- character(0)
#'   }
#'
#'   all_values <- c(all_beta_values, all_random_values)
#'   all_labels <- c(all_beta_labels, all_random_labels)
#'
#'   # Plotting Beta Distribution
#'   all_df <- data.frame(Group = all_labels, Value = all_values)
#'   all_df <- all_df[is.finite(all_df$Value), ]
#'   #all_df <- all_df[-is.nan(all_df$Value), ]
#'
#'
#'   all_dist <- ggplot(all_df, aes(x = Value, fill = Group)) +
#'     geom_histogram(aes(y = after_stat(density)), binwidth = 0.01, alpha = 0.7, position = "identity") +
#'     geom_density(aes(y = after_stat(density)), color = "black", size = 0.25, alpha=0.5) +
#'     scale_fill_brewer(palette = "Set1") +
#'     theme_minimal() +
#'     labs(title = "Histogram of variable importance for each predictor",
#'          x = "Value", y = "Density") +
#'     theme(legend.title = element_blank(),
#'           legend.position = c(0.85, 0.85))
#'
#'   # Plot 2: Histogram with Density Approximation for R2
#'   n_fixed_effects <- ncol(samples_matrix$beta) - 1
#'
#'   all_betas_squared <- rowSums(samples_matrix$beta[, -1]^2)
#'   if (!is.null(samples_matrix$random_var)) {
#'     all_random_vars <- rowSums(samples_matrix$random_var)
#'     x_label_terms <- c(paste0("Lambda^2 beta_", 1:n_fixed_effects, "^2"),
#'                        paste0("sigma[", names(samples_matrix$random_var), "]"))
#'   } else {
#'     all_random_vars <- numeric(length(all_betas_squared))
#'     x_label_terms <- paste0("Lambda^2 beta_", 1:n_fixed_effects, "^2")
#'   }
#'
#'   R_sq_df <- data.frame(Group = rep("R-squared", length(all_betas_squared)),
#'                         Value = all_betas_squared + all_random_vars)
#'
#'   x_label <- paste("Value (", paste(x_label_terms, collapse = " + "), ")")
#'
#'   R_dist <- ggplot(R_sq_df, aes(x = Value, fill = Group)) +
#'     geom_histogram(aes(y = after_stat(density)), binwidth = 0.01, alpha = 0.7, position = "identity") +
#'     geom_density(aes(y = after_stat(density)), color = "black", size = 0.25, alpha=0.5) +
#'     theme_minimal() +
#'     labs(title = "Histogram with Density Approximation for R2",
#'          x = x_label, y = "Density") +
#'     theme(legend.title = element_blank(),
#'           legend.position = c(0.85, 0.85))
#'
#'   # Plot 3: Histogram with Density Approximation for Gelman R2
#'   gelman_R2_df <- data.frame(R2 = samples_matrix$gelman_R2)
#'   gelman_R2_dist <- ggplot(gelman_R2_df, aes(x = R2)) +
#'     geom_histogram(aes(y = after_stat(density), fill = "Gelman R2"), binwidth = 0.004, alpha = 0.7, position = "identity") +
#'     geom_density(aes(y = after_stat(density)), color = "black", size = 0.25, alpha = 0.5) +
#'     theme_minimal() +
#'     labs(title = "Histogram with Density Approximation for Gelman R2",
#'          x = "Value (Gelman R2)", y = "Density")
#'
#'   # Return the plots
#'   list(VariableImportancePlot = all_dist, RDistributionPlot = R_dist, GelmanR2DistributionPlot = gelman_R2_dist)
#' }
#'











