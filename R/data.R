#' Fake questionnaire data
#'
#' A simulated data set containing observations of 100 individuals at four time
#' points. The data was simulated in two groups (50 individuals each) and
#' contains two questionnaires with five items each, one questionnaire with
#' five continuous variables and one additional cross-sectional continuous
#' variable. In this data set the group variable from the simulation is
#' included. You typically don't have this group variable in your data.
#'
#' @format A data frame with 400 rows and 20 variables:
#' \describe{
#'   \item{ID}{patient ID}
#'   \item{visit}{time point of the observation}
#'   \item{group}{to which simulated group the observation belongs to}
#'   \item{age_visit_1}{age of the patient at time point 1}
#'   \item{single_continuous_variable}{a cross-sectional continuous variable,
#'   i.e. there is only one unique value per individual}
#'   \item{questionnaire_A_1}{the first item of questionnaire A with categories
#'   1 to 5}
#'   \item{questionnaire_A_2}{the second item of questionnaire A with categories
#'   1 to 5}
#'   \item{questionnaire_A_3}{the third item of questionnaire A with categories
#'   1 to 5}
#'   \item{questionnaire_A_4}{the fourth item of questionnaire A with categories
#'   1 to 5}
#'   \item{questionnaire_A_5}{the fifth item of questionnaire A with categories
#'   1 to 5}
#'   \item{questionnaire_B_1}{the first item of questionnaire B with categories
#'   1 to 5}
#'   \item{questionnaire_B_2}{the second item of questionnaire B with categories
#'   1 to 5}
#'   \item{questionnaire_B_3}{the third item of questionnaire B with categories
#'   1 to 5}
#'   \item{questionnaire_B_4}{the fourth item of questionnaire B with categories
#'   1 to 5}
#'   \item{questionnaire_B_5}{the fifth item of questionnaire B with categories
#'   1 to 5}
#'   \item{questionnaire_C_1}{the first continuous variable of questionnaire C}
#'   \item{questionnaire_C_2}{the second continuous variable of questionnaire C}
#'   \item{questionnaire_C_3}{the third continuous variable of questionnaire C}
#'   \item{questionnaire_C_4}{the fourth continuous variable of questionnaire C}
#'   \item{questionnaire_C_5}{the fifth continuous variable of questionnaire C}
#' }
#' @source simulated data
"fake_questionnaire_data"
