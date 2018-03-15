
STAVE_fval <- function(X_affinity, y_affinity, lambda, beta_init, ...) {
  res <- STAVE(X_affinity, y_affinity, lambda, beta_init = beta_init, ...)
  return(res$loss_final)
}

beta_init_sampler <- function(X_affinity, y_affinity, num_sub) {
  p <- nrow(X_affinity)
  samp_idx <- sample(1:p, num_sub, prob = abs(y_affinity) / sum(abs(y_affinity)))
  X_aff_sub <- X_affinity[samp_idx, samp_idx]
  y_aff_sub <- y_affinity[samp_idx]
  beta_sub <- solve(X_aff_sub, y_aff_sub)
  beta_samp <- rep(0, p)
  beta_samp[samp_idx] <- beta_sub
  return(beta_samp)
}


#' @export
multi_sampler <- function(X_affinity, y_affinity, lambda, num_sub = lambda, 
                          num_init = 50, iter_init = 5, num_run = 2, seed = 123, init_type = c('sampler', 'rnorm'), parallel = F, ...) {
  init_type <- init_type[1]
  p <- nrow(X_affinity)
  loss_inits <- rep(0, num_init)
  beta_inits <- matrix(0, nrow = p, ncol = num_init)
  beta_fix_init <- rep(0, p)
  dot_inp <- list(...)
  init_inp <- dot_inp
  if (!is.null(dot_inp$max_iter)) {
    init_inp$max_iter <- NULL
  }
  if (is.null(dot_inp$step_size)) {
    init_inp$step_size <- 1 / norm(X_affinity, '2')
  }
  if (parallel) {
    `%fun%` <- doRNG::`%dorng%`
  } else {
    `%fun%` <- foreach::`%do%`
  }
  i <- NULL
  set.seed(seed)
  if (init_type == 'sampler') {
    init_res <- foreach::foreach (i = 1:num_init) %fun% {
      beta_try <- try(beta_init_sampler(X_affinity, y_affinity, num_sub))
      if (class(beta_try) == "try-error" || any(any(is.na(beta_try)) || any(is.infinite(beta_try)))) {
        beta_try <- beta_fix_init
      }
      beta_inits <- beta_try
      loss_inits <- do.call('STAVE_fval', c(list(X_affinity, y_affinity, lambda, beta_try, max_iter = iter_init), init_inp))
      return(list(beta_inits = beta_inits, loss_inits = loss_inits))
    }
  } else if (init_type == 'rnorm') {
    init_res <- foreach::foreach (i = 1:num_init) %fun% {
      beta_inits <- stats::rnorm(p)
      # p <- nrow(X_affinity)
      # samp_idx <- sample(1:p, num_sub, prob = abs(y_affinity) / sum(abs(y_affinity)))
      # beta_inits <- beta_inits[samp_idx]
      loss_inits <- do.call('STAVE_fval', c(list(X_affinity, y_affinity, lambda, beta_inits, max_iter = iter_init), init_inp))
      return(list(beta_inits = beta_inits, loss_inits = loss_inits))
    }
  }
  loss_inits <- purrr::map_dbl(init_res, 'loss_inits')
  beta_inits <- purrr::map(init_res, 'beta_inits')
  ord_idx <- order(loss_inits)
  run_idx <- ord_idx[1:num_run]
  run_beta <- beta_inits[run_idx]
  set.seed(seed)
  run_mods <- foreach::foreach (i = 1:num_run) %fun% {
    res <- STAVE(X_affinity, y_affinity, lambda, beta_init = run_beta[[i]], ...)
    return(res)
  }
  run_conv <- purrr::map_lgl(run_mods, 'converged')
  if (any(run_conv) && any(!run_conv)) {
    run_mods <- run_mods[-which(!run_conv)]
  }
  run_loss <- purrr::map_dbl(run_mods, 'loss_final')
  final_mod <- run_mods[[which.min(run_loss)]]
  return(final_mod)
}
