#' Create the SKEPTIC estimates for Spearman's rho and Kendall's tau.
#' 
#' @param x a numeric vector or matrix.
#' @param y a numeric vector or matrix with compatible dimensions to \code{x}.
#' @param method a string indicating the type of correlation matrix.
#' @param ... additional arguments passed to \code{\link{stats::cor}}.
#' @return A (possibly indefinite) correlation matrix.
#' @export
npn_cor <- function(x, y = NULL, method = c("spearman", "kendall"), ...) {
  method <- method[1]
  cor_est <- stats::cor(x, y, method = method, ...)
  cor_adj <- switch(method,
                    kendall = 2 * sin(cor_est * pi / 6),
                    spearman = sin(cor_est * pi / 2))
  if (is.null(y)) {
    diag(cor_adj) <- 1
  }
  return(cor_adj)
}

#' Run iterative threshodling for Semi-Parametric Affinity Based Variable Selection (STAVE)
#' 
#' @param X_affinity a symmetric matrix of affinities between predictors.
#' @param y_affinity a vector of affinities between each predictor and the response.
#' @param lambda the tuning paramemter. Either a non-negative real number for soft or hard thresholding or a non-negative integer for quantile thresholding.
#' @param y_sc a positive real number representing the scaling value for \code{y}. Default is 1.
#' @param beta_init an initial value for the beta vector. Default is the zero vector.
#' @param eta parameter for ridge type thresholding.
#' @param accelerate use over-relaxation to accelerate the optimization.
#' @param line_search use line search.
#' @param step_size the step-size. The default is the inverse spectral norm of \code{X_affinity} when \code{ls_type} is "forward", or 1 when \code{ls_type} is "backward".
#' @param ls_beta by how much should the step-size increase or decrease in each line search iteration?
#' @param ls_eps you can miss the line search check by this much.
#' @param ls_max_iter the maximum number of line search iterations.
#' @param thresh_type the type of thresholding.
#' @param tol_type the type of tolerance checking.
#' @param max_iter the maximum number of iterations.
#' @param eps the maximum value of the tolerance before declaring convergence.
#' @param ls_type either backtracking (backward) or forward tracking (forward) line search.
#' @return An objecty of type \code{STAVE} with
#' \item{\code{beta_opt}}{optimal beta.}
#' \item{\code{iter}}{number of iterations.}
#' \item{\code{step_size}}{}
#' \item{\code{eta}}{}
#' \item{\code{loss_final}}{final loss value.}
#' \item{\code{loss_vals}}{vector of loss values at each iteration.}
#' \item{\code{converged}}{a boolean.}
#' \item{\code{tot_time}}{total time for iterative thresholding.}
#' \item{\code{nz_patt}}{the indices of estimated non-zero predictors.}
#' \item{\code{lambda}}{}
#' \item{\code{J}}{the estimated cardinality.}
#' \item{\code{IF}}{the inflation factor (for tuning).}
#' \item{\code{DF}}{the degrees of freedom (for tuning).}
#' @export
STAVE <- function(X_affinity, y_affinity, lambda, y_sc = NULL, beta_init = NULL, eta = 1e-4, accelerate = F, line_search = T, 
                        step_size = NULL, ls_beta = 1.5, ls_eps = 1e-10, ls_max_iter = 100, thresh_type = c('quantile_ridge', 'quantile', 'soft', 'hard', 'hard_ridge'), 
                        tol_type = c('beta', 'loss'), max_iter = 2000, eps = 1e-5, ls_type = c('forward', 'backward')) {
  if (nrow(X_affinity) != ncol(X_affinity)) {
    stop('X_affinity must be a square matrix')
  }
  if (!is.null(ncol(y_affinity))) {
    y_affinity <- as.vector(y_affinity) ##make sure y_affinity is a vector
  }
  thresh_type <- thresh_type[1]
  tol_type <- tol_type[1]
  ls_type <- ls_type[1]
  if (is.null(y_sc)) {
    y_sc <- 1
  }
  p <- nrow(X_affinity)
  if (is.null(beta_init)) {
    beta_init <- rep(0, p)
  }
  if (ls_type == 'backward' && ls_beta > 1) {
    ls_beta <- 1 / ls_beta
  }
  if (is.null(step_size) && ls_type == 'backward') {
    step_size <- 1
  }
  if (is.null(step_size) && ls_type == 'forward') {
    step_size <- 1 / base::norm(X_affinity, '2')
  }
  if (thresh_type %in% c('quantile_ridge', 'quantile')) { ## checking the form of lambda - make it an integer if quantile
    lambda <- max(min(floor(lambda), p), 0)
  }
  if (!(thresh_type %in% c('quantile_ridge', 'hard_ridge'))) {
    eta <- 0
  }
  start_tm <- Sys.time()
  loop_call <- TISP_inner(X_affinity, y_affinity, beta_init, y_sc = as.numeric(y_sc), lambda = lambda, q = lambda, eta = eta, step_size = step_size, accelerate = accelerate,
                            line_search = line_search, ls_beta = ls_beta, ls_eps = ls_eps, ls_maxit = ls_max_iter, thresh_type = thresh_type, tol_type = tol_type, maxit = max_iter,
                            eps = eps, ls_type = ls_type)
  end_tm <- Sys.time()
  tot_time <- end_tm - start_tm
  beta_opt <- loop_call$beta_opt
  nz_patt <- which(abs(beta_opt) > .Machine$double.eps)
  J <- length(nz_patt)
  IF <- J * (1 + log(p) - log(J))
  DF <- J
  # if (thresh_type %in% c('quantile_ridge', 'hard_ridge') && !is.null(X)) {
  #   X_sub <- X[, nz_patt]
  #   X_sub_affinity <- crossprod(X_sub)
  #   DF <- sum(diag(solve(X_sub_affinity + loop_call$eta * diag(J), X_sub_affinity)))
  # } else {
  #   DF <- J
  # }
  ret <- c(loop_call, list(tot_time = as.numeric(tot_time, unit = 'secs'), nz_patt = nz_patt, lambda = lambda, J = J, IF = IF, DF = DF))
  class(ret) <- 'STAVE'
  return(ret)
}

#' Print method for STAVE
#' @param x a \code{STAVE} object.
#' @param ... not used.
#' @export
print.STAVE <- function(x, ...) {
  cat('Estimated cardinality:', x$J,
      '\nIterations:', x$iter,
      '\nConverged:', x$converged,
      '\nTolerance:', x$tol,
      '\nTotal time (s):', x$tot_time)
}
