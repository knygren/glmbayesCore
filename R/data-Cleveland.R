#' Cleveland Heart Disease Dataset
#'
#' A cleaned version of the Cleveland heart disease dataset from the UCI Machine
#' Learning Repository. This version contains only complete cases and includes a
#' derived binary outcome variable \code{hd} indicating the presence (\code{"Yes"})
#' or absence (\code{"No"}) of heart disease.
#'
#' @format A data frame with 297 observations and 15 variables:
#' \describe{
#'   \item{age}{Age in years (numeric).}
#'   \item{sex}{Sex (0 = female, 1 = male).}
#'   \item{cp}{Chest pain type (numeric code 1--4).}
#'   \item{trestbps}{Resting blood pressure (mm Hg).}
#'   \item{chol}{Serum cholesterol (mg/dl).}
#'   \item{fbs}{Fasting blood sugar > 120 mg/dl (1 = true, 0 = false).}
#'   \item{restecg}{Resting electrocardiographic results (numeric code).}
#'   \item{thalach}{Maximum heart rate achieved.}
#'   \item{exang}{Exercise-induced angina (1 = yes, 0 = no).}
#'   \item{oldpeak}{ST depression induced by exercise relative to rest.}
#'   \item{slope}{Slope of the peak exercise ST segment.}
#'   \item{ca}{Number of major vessels colored by fluoroscopy (0--3).}
#'   \item{thal}{Thalassemia status (numeric code).}
#'   \item{num}{Original UCI disease score (0--4).}
#'   \item{hd}{Binary heart disease indicator: \code{"No"} (num = 0) or \code{"Yes"} (num > 0).}
#' }
#'
#' @source UCI Machine Learning Repository: Heart Disease Data Set.
#' \url{https://archive.ics.uci.edu/dataset/45/heart+disease}
#'
#' @usage data("Cleveland")
#'
#' @example inst/examples/Ex_Cleveland_opencl.R
#'
#' @keywords datasets
#' @concept Bayesian Binomial Regression
"Cleveland"
