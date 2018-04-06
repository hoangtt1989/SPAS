Semi-Parametric Affinity Based Variable Selection
================
Hoang Tran
3/16/2018

Installation
------------

``` r
devtools::install_github('hoangtt1989/SPAS')
```

Simulated Data Example
----------------------

``` r
library(SPAS)
library(MASS)
set.seed(123)
n <- 200
p <- 500
J <- 5
cor_str <- .1

cor_mat <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      cor_mat[i, j] <- cor_str^(abs(i - j))
    }
  }
X <- mvrnorm(n, rep(0, p), cor_mat)
beta <- rnorm(p) * 5
beta[(J+1):p] <- 0
y <- X %*% beta + rnorm(n)

X_trans <- sinh(X)
y_trans <- sinh(y)

X_affinity <- cor(X_trans, method = 'spearman')
y_affinity <- cor(X_trans, y_trans, method = 'spearman')

SPAS_fit <- SPAS(X_affinity, y_affinity, 10)
SPAS_fit
```

    ## Estimated cardinality: 10 
    ## Iterations: 122 
    ## Converged: TRUE 
    ## Tolerance: 9.023684e-06 
    ## Total time (s): 0.03906298

What percent of the true predictors were selected?

``` r
true_nz <- which(abs(beta - 0) > 1e-10)
length(intersect(true_nz, SPAS_fit$nz_patt)) / J * 100
```

    ## [1] 80
