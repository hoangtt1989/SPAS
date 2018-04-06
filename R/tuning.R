#' Internal PIC function - not meant to be directly accessed by users.
#' 
#' @param loss the squared error loss.
#' @param DF degrees of freedom.
#' @param IF inflation factor.
#' @param mn the number of observations (univariate), or the number of responses times the number of observations (multivariate).
#' @param const1 the constant for DF.
#' @param const2 the constant for IF.
#' @export
PIC_inner <- function(loss, DF, IF, mn, const1 = 1.5, const2 = .8) {
  complexity <- (DF * const1 + IF * const2) / mn
  if (complexity >= 1 || any(is.na(complexity) || is.nan(complexity) || is.infinite(complexity))) {
    val <- Inf
  } else {
    val <- loss / (1 - complexity)
    # val <- switch(type,
    #               additive = log(loss) + complexity,
    #               fractional = loss / (1 - complexity),
    #               GCV = loss / (1 - complexity)^2,
    #               plug = loss + complexity * loss)
  }
  return(val)
}
#' Generic method for PIC (fractional form)
#' 
#' @param object an object with a PIC method.
#' @param n the number of observations.
#' @param const1 the constant for DF.
#' @param const2 the constant for IF.
#' @param ... additional arguments passed to methods.
#' @export
PIC <- function(object, n, const1 = 1.5, const2 = .8, ...) {
  UseMethod('PIC', object)
}
#' SPAS method for PIC
#' 
#' @param object an object with a PIC method.
#' @param n the number of observations.
#' @param const1 the constant for DF.
#' @param const2 the constant for IF.
#' @param ... not used.
#' @export
PIC.SPAS <- function(object, n, const1 = 1.5, const2 = .8, ...) {
  loss <- object$loss_final
  DF <- object$DF
  IF <- object$IF
  PIC_inner(loss, DF, IF, n, const1, const2)
}
