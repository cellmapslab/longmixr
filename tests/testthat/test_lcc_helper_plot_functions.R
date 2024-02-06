set.seed(5)
test_data <- data.frame(patient_id = rep(1:15, each = 4),
                        visit = rep(1:4, 15),
                        var_1 = c(rnorm(20, -1), rnorm(20, 3), rnorm(20, 0.5)) +
                          rep(seq(from = 0, to = 1.5, length.out = 4), 15),
                        var_2 = c(rnorm(20, 0.5, 1.5), rnorm(20, -2, 0.3),
                                  rnorm(20, 4, 0.5)) +
                          rep(seq(from = 1.5, to = 0, length.out = 4), 15),
                        cat_1 = sample(LETTERS[1:3], 60, replace = TRUE))
model_list <- list(flexmix::FLXMRmgcv(as.formula("var_1 ~ .")),
                   flexmix::FLXMRmgcv(as.formula("var_2 ~ .")))

clustering <- longitudinal_consensus_cluster(
  data = test_data,
  id_column = "patient_id",
  max_k = 3,
  reps = 3,
  model_list = model_list
)

# helper function to save the plots
save_png <- function(code, width = 800, height = 800) {
  path <- tempfile(fileext = ".png")
  png(path, width = width, height = height)
  on.exit(dev.off())
  code

  path
}

################################################################################
# plot_spaghetti

test_that("plot_spaghetti generates a plot in the CI", {
  skip_if_not(isTRUE(as.logical(Sys.getenv("CI"))),
              message = "plot_spaghetti simple plot generation only tested on CI")
  expect_s3_class(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit"
  ),
  c("gg", "ggplot"))
})

test_that("the plots stay the same", {
  skip_on_ci()
  expect_snapshot_file(save_png(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit"
  )), "plot_spaghetti_output.png")
})

test_that("plot_spaghetti only uses lcc objects", {
  expect_error(plot_spaghetti(
    model = 1:10,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit"
  ))
})

test_data_2 <- test_data
test_data_2$assignment_num_clus_2 <- 1
test_data_3 <- test_data
test_data_3$var_1 <- NULL
test_data_3$cat_1 <- NULL
test_data_4 <- test_data
test_data_4$visit <- NULL

test_that("the data argument is correct and checks for required columns", {
  expect_error(plot_spaghetti(
    model = clustering,
    data = cbind(1:10, 11:20),
    variable_names = "var_1",
    time_variable = "visit"
  ))

  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data_3,
    variable_names = "var_1",
    time_variable = "visit"
  ))

  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data_4,
    variable_names = "var_1",
    time_variable = "visit"
  ))
})

test_that("the variable_names argument is correct", {
  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_10",
    time_variable = "visit"
  ))

  expect_s3_class(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = c("var_1", "var_2"),
    time_variable = "visit"),
    c("gg", "ggplot"))

  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = c("var_1", "var_10"),
    time_variable = "visit"
  ))
})

test_that("the time_variable argument is correct", {
  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "time"
  ))

  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = c("visit", "time")
  ))
})

test_that("the show_mean_sd_ribbon argument is correct", {
  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit",
    show_mean_sd_ribbon = "no"
  ))

  test_plot_1 <- plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit",
    show_mean_sd_ribbon = TRUE
  )

  test_plot_2 <- plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit",
    show_mean_sd_ribbon = FALSE
  )

  expect_equal(all.equal(test_plot_1$labels, test_plot_2$labels),
               "Length mismatch: comparison on first 4 components")
})

test_that("the number_of_clusters argument is correct", {
  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit",
    number_of_clusters = "three"
  ))

  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit",
    number_of_clusters = 5
  ))

  expect_error(plot_spaghetti(
    model = clustering,
    data = test_data,
    variable_names = "var_1",
    time_variable = "visit",
    number_of_clusters = c(2, 3)
  ))
})

test_that("existing clustering result columns in the provided data are removed", {
  expect_warning(plot_spaghetti(
    model = clustering,
    data = test_data_2,
    variable_names = "var_1",
    time_variable = "visit",
    number_of_clusters = 2
  ))
})

################################################################################
# plot_alluvial

test_that("plot_alluvial generates a plot in the CI", {
  skip_if_not(isTRUE(as.logical(Sys.getenv("CI"))),
              message = "plot_spaghetti simple plot generation only tested on CI")
  expect_s3_class(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = "cat_1",
    time_variable = "visit"
  ),
  c("gg", "ggplot"))
})

test_that("the plots stay the same", {
  skip_on_ci()
  expect_snapshot_file(save_png(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = "cat_1",
    time_variable = "visit"
  )), "plot_alluvial_output.png")
})

test_that("plot_alluvial only uses lcc objects", {
  expect_error(plot_alluvial(
    model = 1:10,
    data = test_data,
    variable_name = "cat_1",
    time_variable = "visit"
  ))
})

test_that("the data argument is correct and checks for required columns", {
  expect_error(plot_alluvial(
    model = clustering,
    data = cbind(1:10, 11:20),
    variable_name = "cat_1",
    time_variable = "visit"
  ))

  expect_error(plot_alluvial(
    model = clustering,
    data = test_data_3,
    variable_name = "cat_1",
    time_variable = "visit"
  ))

  expect_error(plot_alluvial(
    model = clustering,
    data = test_data_4,
    variable_name = "cat_1",
    time_variable = "visit"
  ))
})

test_that("the variable_names argument is correct", {
  expect_error(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = "cat_10",
    time_variable = "visit"
  ))

  expect_error(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = c("cat_1", "var_2"),
    time_variable = "visit"
  ))
})

test_that("the time_variable argument is correct", {
  expect_error(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = "cat_1",
    time_variable = "time"
  ))

  expect_error(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = "cat_1",
    time_variable = c("visit", "time")
  ))
})

test_that("the number_of_clusters argument is correct", {
  expect_error(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = "cat_1",
    time_variable = "visit",
    number_of_clusters = "three"
  ))

  expect_error(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = "cat_1",
    time_variable = "visit",
    number_of_clusters = 5
  ))

  expect_error(plot_alluvial(
    model = clustering,
    data = test_data,
    variable_name = "cat_1",
    time_variable = "visit",
    number_of_clusters = c(2, 3)
  ))
})

test_that("existing clustering result columns in the provided data are removed", {
  expect_warning(plot_alluvial(
    model = clustering,
    data = test_data_2,
    variable_name = "cat_1",
    time_variable = "visit",
    number_of_clusters = 2
  ))
})
