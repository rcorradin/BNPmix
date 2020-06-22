#' Collaborative Perinatal Project data
#'
#' @description
#' A subset of the Collaborative Perinatal Project data set (Klebanoff, 2009)
#' focused on studying the effect of DDE exposure on pregnancies (Longnecker et al., 2001).
#' The dataset contains the following variables for each pregnant women enrolled in the study:
#'\itemize{
#'   \item hosp, factor denoting the hospital where the woman was hospitalized;
#'   \item smoke, factor. It takes value 2 if the woman is a smoker, 1 otherwise;
#'   \item gest, gestational age (in weeks);
#'   \item dde, Dichlorodiphenyldichloroethylene (DDE) concentration in maternal serum;
#'   \item weight, body weight of the baby at birth (in grams);
#' }
#'
#' @docType data
#' @keywords dataset internal
#' @name CPP
#'
#' @usage data(CPP)
#' @format A data.frame
#'
#' @examples
#' data(CPP)
#' str(CPP)
#'
#' @references
#'
#' Klebanoff M. A. (2009) The collaborative perinatal project: a 50-year retrospective.
#' Paediatric and perinatal epidemiology, 23, 2.
#'
#' Longnecker, M. P., Klebanof, M. A., Zhou, H., Brock, J. (2001)
#' Association between maternal serum concentration of the DDT metabolite
#' DDE and preterm and small-for-gestational-age babies at birth. The Lancet, 358, 110-114.
#'
"CPP"
