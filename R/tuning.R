PIC_inner <- function(loss, DF, IF, mn, type = 'fractional', const1 = 1.5, const2 = .8) {
  complexity <- (DF * const1 + IF * const2) / mn
  if (complexity >= 1 || any(is.na(complexity) || is.nan(complexity) || is.infinite(complexity))) {
    val <- Inf
  } else {
    val <- switch(type,
                  additive = log(loss) + complexity,
                  fractional = loss / (1 - complexity),
                  GCV = loss / (1 - complexity)^2,
                  plug = loss + complexity * loss)
  }
  return(val)
}
#' Generic method for PIC
#' @export
PIC <- function(object, n, type = 'fractional', const1 = 1.5, const2 = .8, ...) {
  UseMethod('PIC', object)
}
#' STAVE method for PIC
#' @export
PIC.STAVE <- function(object, n, type = 'fractional', const1 = 1.5, const2 = .8, ...) {
  loss <- object$loss_final
  DF <- object$DF
  IF <- object$IF
  PIC_inner(loss, DF, IF, n, type, const1, const2)
}
