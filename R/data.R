#' Fake questionnaire data
#'
#' A simulated data set containing observations of 100 individuals at 4 time
#' points. The data was simulated in two groups (50 individuals each) and
#' contains 2 questionnaires with 5 items each and a continuous variable.
#'
#' @format A data frame with 400 rows and 14 variables:
#' \describe{
#'   \item{ID}{patient ID}
#'   \item{visit}{time point of the observation}
#'   \item{group}{to which simulated group the observation belongs to}
#'   \item{continuous_variable}{a continuous variable}
#'   \item{questionnaire_A_1}{the first item of questionnaire A with categories
#'   1 to 5; the same for items 2-5}
#'   \item{questionnaire_B_1}{the first item of questionnaire B with categories
#'   1 to 5; the same for items 2-5}
#' }
#' @source simulated data
"fake_questionnaire_data"
