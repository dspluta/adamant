
#' Adaptive Mantel test for linear association.
#'
#' \code{adamant} Returns a dataframe containing results of the adaptive
#' Mantel test for linear association between two sets of variables.
#'
#' @param X A matrix of covariates.
#' @param Y A vector of response values.
#' @return An object of class `adamant` that contains the test results.
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' Y <- rnorm(20)
#' lambdas_X <- c(0, 1, 1000, Inf)
#' adamant(X, Y, lambdas_X)
#' }
#' @import dplyr
#' @import purrr
#' @import magrittr
#' @export
adamant <- function(X, Y, lambdas_X = c(0), n_perms = 100, verbose = TRUE,
                         P_val_only = FALSE) {
  # ---------------------------------------------------------------------------
  # START Test Calculations
  # ---------------------------------------------------------------------------
  start_time <- Sys.time()

  # Set up vars for computation, following notation in Pluta et al. 2017 ------
  n_lambda <- length(lambdas_X)
  n_perms_p1 <- n_perms + 1
  n <- nrow(X)
  X <- scale(X)
  Y <- scale(Y)
  X_svd <- svd(X)
  U <- X_svd$u
  d <- X_svd$d
  Z <- crossprod(U, Y)

  # Generate needed number of permutations for Z ------------------------------
  #   First element of Z_perm is Z from original data
  #   Elements 2:(n_perms + 1) are Z calculated from permuting rows of U
  Z_perm <- prepend(map(1:n_perms,
                        function(b) {crossprod(U[sample(1:nrow(U)), ] , Y)}),
                    list(Z))

  # Calculate weights for each lambda in lambdas_X ----------------------------
  #   Set weights separately for lambda = Inf since it's special

  weights <- map(lambdas_X, ~d^2/(.x + d^2))
  if (Inf %in% lambdas_X)
    weights[[which(lambdas_X == Inf)]] <- d^2

  # Define `main` -------------------------------------------------------------
  #   `main` is used to calculate the adaptive P_val
  #   `main` is the returned data frame when P_val_only = FALSE
  #   when returned, `main` includes P_val_lambda^(b) for
  #   all lambda and b = 1:(n_perms + 1)
  main <- tibble(lambda = rep(lambdas_X, each = n_perms_p1),
                 b = rep(1:(n_perms_p1), n_lambda), MT = NA)
  for (i in 1:(n_lambda)) {
    main$MT[((i - 1) * n_perms_p1 + 1):(i * (n_perms_p1))] <-
      unlist(map(Z_perm, ~crossprod(weights[[i]]^2, .x^2)))
  }

  # Compute P_val_lambda for all lambda and all perms -------------------------
  main <- main %>%
    group_by(lambda) %>%
    mutate(P_val_lambda = 1 - percent_rank(MT))

  # Calculate the adaptive Mantel P_val ---------------------------------------
  P_val <- main %>%
    group_by(b) %>%
    summarize(P_val_b = min(P_val_lambda)) %>%
    transmute(b = b, P_val = percent_rank(P_val_b)) %>%
    filter(b == 1) %$%
    P_val

  # Identify the lambda that yields the P_val ---------------------------------
  lambda <- main %>% filter(b == 1) %>%
    group_by(b) %>%
    slice(which.min(P_val_lambda)) %$%
    lambda

  end_time <- Sys.time()
  #----------------------------------------------------------------------------
  # END  Test Calculations
  #----------------------------------------------------------------------------

  # Print output if verbose = TRUE --------------------------------------------
  time_taken <- (end_time - start_time)
  if (verbose) {
    print(glue::glue(
      '-------------------------------\n',
      'Adaptive Mantel Output\n',
      '-------------------------------\n',
      'P_val         = {P_val}\n',
      'n             = {n}\n',
      'p             = {p}\n',
      'rank(X^TX)    = {r}\n',
      'kappa         = {kappa}\n',
      'Best Lambda   = {lambda}\n',
      'n_perms       = {n_perms}\n',
      'time          = {time_taken} {units}\n',
      '-------------------------------\n\n',
      p = ncol(X), r = sum(round(d, digits = 12) != 0),
      kappa = round(max(d) / min(d[d > 0]), digits = 3),
      units = attr(time_taken, "units"),
      time_taken = round(time_taken, digits = 3)))
  }

  # Return P_val if P_val_only = TRUE, else return main -----------------------
  if (P_val_only) {
    return(P_val)
  } else {
    return(list(results = main, P_val = P_val, lambda = lambda,
                n_perms = n_perms, lambdas_X = lambdas_X))
  }
}
