# ============================================================
# Kernel functions
# ============================================================

#' Epanechnikov kernel
#'
#' Evaluates the Epanechnikov kernel
#'
#'   K(u) = 3/4 (1 - u^2) I(|u| <= 1).
#'
#' @param u Numeric vector of standardized distances.
#'
#' @return Numeric vector of kernel values.
epanechnikov_kernel <- function(u) {
  (3 / 4) * (1 - u^2) * (abs(u) <= 1)
}


#' Integrated squared Epanechnikov kernel
#'
#' Computes
#'
#'   R(K) = integral K(u)^2 du.
#'
#' For the Epanechnikov kernel the analytical value is 3/5.
RK_epanechnikov <- function() {

  3 / 5

}


#' Second moment of the Epanechnikov kernel
#'
#' Computes
#'
#'   mu_2(K) = integral u^2 K(u) du.
#'
#' For the Epanechnikov kernel the analytical value is 1/5.
mu2_epanechnikov <- function() {

  1 / 5

}