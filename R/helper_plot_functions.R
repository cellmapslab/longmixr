#' Spaghetti plot for longmixr clusterings
#'
#' A helper function to plot spaghetti plots of continuous variables separated
#' by the clusters found by \code{longmixr}.
#'
#' @param model \code{lcc} object (output from \code{\link{longitudinal_consensus_cluster}})
#' @param data a \code{data.frame} that contains the variables to be plotted and
#' the time and ID variable used in the longmixr clustering; typically the data
#' used for the clustering
#' @param variable_names character vector of the continuous variables to be plotted
#' @param time_variable the name of the variable that depicts the time point of
#' the measurements
#' @param show_mean_sd_ribbon \code{boolean} if the mean and SD per variable
#' should be shown, the default is \code{TRUE}
#' @param number_of_clusters the number of clusters that should be plotted, the
#' default is \code{2}
#' @param scales \code{scales} argument of \code{facet_wrap}, the default is
#' \code{fixed}
#'
#' @details
#' The spaghetti plot shows the longitudinal trajectory (defined by
#' \code{time_variable}) of continuous variables separated by the clusters found
#' by \code{\link{longitudinal_consensus_cluster}}. The provided \code{data.frame}
#' for \code{data} can either be the same as used in the clustering with
#' \code{\link{longitudinal_consensus_cluster}} or needs to contain the same
#' \code{id_column} as in the clustering and a \code{time_variable}.
#'
#' @return a \code{ggplot} object that is plotted
#'
#' @importFrom stats aggregate reshape sd
#' @importFrom ggplot2 .data
#'
#' @export
#'
#' @examples
#' set.seed(5)
#' test_data <- data.frame(patient_id = rep(1:10, each = 4),
#' visit = rep(1:4, 10),
#' var_1 = c(rnorm(20, -1), rnorm(20, 3)) +
#' rep(seq(from = 0, to = 1.5, length.out = 4), 10),
#' var_2 = c(rnorm(20, 0.5, 1.5), rnorm(20, -2, 0.3)) +
#' rep(seq(from = 1.5, to = 0, length.out = 4), 10))
#' model_list <- list(flexmix::FLXMRmgcv(as.formula("var_1 ~ .")),
#' flexmix::FLXMRmgcv(as.formula("var_2 ~ .")))
#' clustering <- longitudinal_consensus_cluster(
#' data = test_data,
#' id_column = "patient_id",
#' max_k = 2,
#' reps = 3,
#' model_list = model_list,
#' flexmix_formula = as.formula("~s(visit, k = 4) | patient_id"))
#'
#' plot_spaghetti(
#'   model = clustering,
#'   data = test_data,
#'   variable_names = "var_1",
#'   time_variable = "visit"
#' )
plot_spaghetti <- function(model,
                           data,
                           variable_names,
                           time_variable,
                           show_mean_sd_ribbon = TRUE,
                           number_of_clusters = 2,
                           scales = "fixed"
) {
  checkmate::assert_class(model, "lcc")
  checkmate::assert_class(data, "data.frame")
  checkmate::assert_character(variable_names)
  checkmate::assert_true(all(variable_names %in% colnames(data)))
  checkmate::assert_character(time_variable, len = 1)
  checkmate::assert_true(time_variable %in% colnames(data))
  checkmate::assert_logical(show_mean_sd_ribbon, len = 1)
  checkmate::assert_int(number_of_clusters, lower = 2,
                        upper = model$general_information$call$max_k)

  plot_data <- reshape(data = data,
                       varying = variable_names,
                       v.names = "value",
                       timevar = "variable",
                       times = variable_names,
                       direction = "long")

  # get the name of the id column
  id_column_name <- model$general_information$call$id_column

  # check if the provided data already contains columns with the name pattern
  # assignment_num_clus_
  # if yes -> remove them and use the ones from the lcc model
  vector_test_colnames <- grepl("^assignment_num_clus_[0-9]*", colnames(plot_data))
  if (any(vector_test_colnames)) {
    plot_data[, vector_test_colnames] <- NULL
    warning("The provided data contains already columns with the pattern 'assignment_num_clus_'. These columns were removed and replaced by the cluster information from the longmixr model.")
  }

  # add cluster information
  plot_data <- merge(plot_data, model$general_information$cluster_assignments,
                     by = id_column_name, all.x = TRUE)

  # this step to create a patient_id per variable is needed, otherwise
  # ggplot can't differentiate between the different variables
  plot_data[[id_column_name]] <- paste0(plot_data[[id_column_name]], "_",
                                        plot_data$variable)

  # determine which clustering solution should be shown
  clus_assign <- paste0("assignment_num_clus_", number_of_clusters)

  # create factors for correct plotting
  plot_data[[time_variable]] <- as.factor(plot_data[[time_variable]])
  plot_data$variable <- as.factor(plot_data$variable)
  plot_data[[id_column_name]] <- as.factor(plot_data[[id_column_name]])
  plot_data[[clus_assign]] <- as.factor(plot_data[[clus_assign]])

  # calculate mean and SD per variable/time point/cluster
  data_for_aggregation <- plot_data[, "value", drop = FALSE]
  by_variables <- list(plot_data[[time_variable]],
                       plot_data$variable,
                       plot_data[[clus_assign]])
  names(by_variables) <- c(time_variable, "variable", clus_assign)

  additional_data <- aggregate(
    x = data_for_aggregation,
    by = by_variables,
    FUN = function(x) {
      cbind(
        mean = mean(x),
        sd = sd(x)
      )
    }
  )
  # the above returns a data.frame with the results as a matrix within one
  # data.frame column. Therefore, I need to expand it
  additional_data <- do.call(data.frame, additional_data)

  colnames(additional_data)[colnames(additional_data) == "value.1"] <- "mean"
  colnames(additional_data)[colnames(additional_data) == "value.2"] <- "sd"

  additional_data[[id_column_name]] <- 1

  p <- ggplot2::ggplot(data = plot_data) +
    ggplot2::geom_line(mapping = ggplot2::aes(x = .data[[time_variable]],
                                              y = .data[["value"]],
                                              col = .data[["variable"]],
                                              group = .data[[id_column_name]]),
                       alpha = 0.4) +
    ggplot2::facet_wrap(~.data[[clus_assign]], scales = scales)

  if (show_mean_sd_ribbon) {
    p <- p +
      ggplot2::geom_ribbon(data = additional_data,
                           mapping = ggplot2::aes(x = .data[[time_variable]],
                                                  y = .data[["mean"]],
                                                  ymin = .data[["mean"]] - .data[["sd"]],
                                                  ymax = .data[["mean"]] + .data[["sd"]],
                                                  fill = .data[["variable"]],
                                                  group = .data[["variable"]]),
                           alpha = 0.3) +
      ggplot2::geom_line(data = additional_data,
                         mapping = ggplot2::aes(x = .data[[time_variable]],
                                                y = .data[["mean"]],
                                                col = .data[["variable"]],
                                                group = .data[["variable"]]),
                         size = 2)
  }
  p
}

#' Alluvial plot for longmixr clusterings
#'
#' A helper function to plot alluvial plots of a categorical variable separated
#' by the clusters found by \code{longmixr}. You need to have
#' \code{ggalluvial} installed to use this function.
#'
#' @param model model \code{lcc} object (output from \code{\link{longitudinal_consensus_cluster}})
#' @param data a \code{data.frame} that contains the variables to be plotted and
#' the time and ID variable used in the longmixr clustering; typically the data
#' used for the clustering
#' @param variable_name name of the categorical variable to be plotted as character
#' @param time_variable the name of the variable that depicts the time point of
#' the measurements
#' @param number_of_clusters the number of clusters that should be plotted, the
#' default is \code{2}
#'
#' @return a \code{ggplot} object that is plotted
#'
#' @importFrom ggplot2 .data
#'
#' @export
#'
#' @examples
#' library(ggalluvial)
#' set.seed(5)
#' test_data <- data.frame(patient_id = rep(1:10, each = 4),
#' visit = rep(1:4, 10),
#' var_1 = c(rnorm(20, -1), rnorm(20, 3)) +
#' rep(seq(from = 0, to = 1.5, length.out = 4), 10),
#' var_2 = c(rnorm(20, 0.5, 1.5), rnorm(20, -2, 0.3)) +
#' rep(seq(from = 1.5, to = 0, length.out = 4), 10))
#' model_list <- list(flexmix::FLXMRmgcv(as.formula("var_1 ~ .")),
#' flexmix::FLXMRmgcv(as.formula("var_2 ~ .")))
#' clustering <- longitudinal_consensus_cluster(
#' data = test_data,
#' id_column = "patient_id",
#' max_k = 2,
#' reps = 3,
#' model_list = model_list,
#' flexmix_formula = as.formula("~s(visit, k = 4) | patient_id"))
#'
#' # add categorical variable for test plotting
#' test_data$cat <- sample(LETTERS[1:3], 40, replace = TRUE)
#'
#' plot_alluvial(
#'   model = clustering,
#'   data = test_data,
#'   variable_name = "cat",
#'   time_variable = "visit"
#' )
plot_alluvial <- function(model,
                          data,
                          variable_name,
                          time_variable,
                          number_of_clusters = 2
) {
  checkmate::assert_class(model, "lcc")
  checkmate::assert_class(data, "data.frame")
  checkmate::assert_character(variable_name, len = 1)
  checkmate::assert_true(all(variable_name %in% colnames(data)))
  checkmate::assert_character(time_variable, len = 1)
  checkmate::assert_true(time_variable %in% colnames(data))
  checkmate::assert_int(number_of_clusters, lower = 2,
                        upper = model$general_information$call$max_k)

  # get the name of the id column
  id_column_name <- model$general_information$call$id_column

  # check if the provided data already contains columns with the name pattern
  # assignment_num_clus_
  # if yes -> remove them and use the ones from the lcc model
  vector_test_colnames <- grepl("^assignment_num_clus_[0-9]*", colnames(data))
  if (any(vector_test_colnames)) {
    data[, vector_test_colnames] <- NULL
    warning("The provided data contains already columns with the pattern 'assignment_num_clus_'. These columns were removed and replaced by the cluster information from the longmixr model.")
  }

  # add cluster information
  plot_data <- merge(data, model$general_information$cluster_assignments,
                     by = id_column_name, all.x = TRUE)

  # make sure that the variable to be plotted is a factor
  plot_data[[variable_name]] <- as.factor(plot_data[[variable_name]])

  # determine which clustering solution should be shown
  clus_assign <- paste0("assignment_num_clus_", number_of_clusters)

  ggplot2::ggplot(plot_data) +
    ggplot2::aes(x = .data[[time_variable]],
                 stratum = .data[[variable_name]],
                 alluvium = .data[[id_column_name]],
                 fill = .data[[variable_name]],
                 label = .data[[variable_name]]) +
    ggplot2::scale_fill_brewer(type = "qual", palette = "Set2") +
    ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "frontback",
                          color = "darkgray") +
    ggalluvial::geom_stratum() +
    ggplot2::facet_wrap(~as.factor(.data[[clus_assign]])) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ylab("n") +
    ggplot2::ggtitle(paste0("Distribution of ", variable_name, " across clusters"))
}
