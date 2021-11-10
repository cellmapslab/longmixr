#' Fake questionnaire data
#'
#' A simulated data set containing observations of 100 individuals at four time
#' points. The data was simulated in two groups (50 individuals each) and
#' contains two questionnaires with five items each, one questionnaire with
#' three continuous variables and one additional continuous variable.
#'
#' @format A data frame with 400 rows and 18 variables:
#' \describe{
#'   \item{ID}{patient ID}
#'   \item{visit}{time point of the observation}
#'   \item{group}{to which simulated group the observation belongs to}
#'   \item{age_visit_1}{age of the patient at time point 1}
#'   \item{single_continuous_variable}{a continuous variable}
#'   \item{questionnaire_A_1}{the first item of questionnaire A with categories
#'   1 to 5; the same for items 2-5}
#'   \item{questionnaire_B_1}{the first item of questionnaire B with categories
#'   1 to 5; the same for items 2-5}
#'   \item{questionnaire_C_1}{the first continuous variable of questionnaire C;
#'   the same for variables 2-3}
#' }
#' @source simulated data
"fake_questionnaire_data"
