#' Calibration targets for OUD model
#'
#' A list with calibration targets
#' @format A list with three calibration targets:
#' \describe{
#'   \item{deaths}{Overall cohort mortality. A data frame with 3 rows and 8 variables:
#'     \itemize{\item Target: Target name
#'              \item Time: Time in months
#'              \item Num: Number of deaths
#'              \item Pop: Population at risk
#'              \item value: Target value
#'              \item se: Standard error
#'              \item lb: 95\% CI lower bound
#'              \item ub: 95\% CI lower bound} }
#'   \item{overdose}{Cumulative non-fatal overdoses. A data frame with 3 rows and 8 variables:
#'     \itemize{\item Target: Target name
#'              \item Time: Time in months
#'              \item Num: Number of individuals in overdose state 
#'              \item Pop: Population at risk
#'              \item value: Target value
#'              \item se: Standard error
#'              \item lb: 95\% CI lower bound
#'              \item ub: 95\% CI lower bound} }
#' }
#' @md
"calibration_targets"