#' Convert z statistics to correlations.
#' @export
zstat_to_cor <- function(zstat, n) {
  tmp_val <- zstat / sqrt(n - 1)
  sign(tmp_val) * abs(tmp_val) / sqrt(1 + tmp_val^2)
}
