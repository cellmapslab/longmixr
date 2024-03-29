set.seed(5)
test_data <- data.frame(patient_id = rep(1:15, each = 4),
                        visit = rep(1:4, 15),
                        var_1 = c(rnorm(20, -1), rnorm(20, 3), rnorm(20, 0.5)) +
                          rep(seq(from = 0, to = 1.5, length.out = 4), 15),
                        var_2 = c(rnorm(20, 0.5, 1.5), rnorm(20, -2, 0.3),
                                  rnorm(20, 4, 0.5)) +
                          rep(seq(from = 1.5, to = 0, length.out = 4), 15))
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
# because plot.lcc outputs several plots at once, I store them all in one big
# png and compare this png
save_png <- function(code, width = 800, height = 800) {
  path <- tempfile(fileext = ".png")
  png(path, width = width, height = height)
  op <- par(mfrow = c(4, 2))
  on.exit(dev.off())
  on.exit(par(op), add = TRUE)
  code

  path
}

test_that("a plot is generated in the CI", {
  skip_if_not(isTRUE(as.logical(Sys.getenv("CI"))),
              message = "simple plot generation only tested on CI")
  test_plot <- plot(clustering)
  expect_equal(test_plot, matrix(c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 7.9),
                                 nrow = 7))
})

test_that("the plots stay the same", {
  skip_on_ci()
  # because plot.lcc changes the graphics parameters, currently only the last
  # plot is recorded. Therefore, I can only check the last plot.
  expect_snapshot_file(save_png(plot(clustering)), "plot_lcc_output.png")
})

test_that("plot.lcc only uses lcc objects", {
  expect_error(plot.lcc(1:10))
})

test_that("the color_palette argument is correct", {
  expect_error(plot(clustering, color_palette = 1:10),
               regexp = "Assertion on 'color_palette' failed")

  expect_silent(plot(clustering,
                     color_palette = c("#9999FF", "#7F7FFF", "#6666FF")))
})

test_that("the which_plots argument is correct", {
  expect_error(plot(clustering, which_plots = 3))
  expect_error(plot(clustering, which_plots = "consensusmatrix_4"),
               regexp = "which_plot must be one of all, consensusmatrix_legend, consensusmatrix, consensusmatrix_2, consensusmatrix_3, CDF, delta, cluster_tracking, item_consensus, cluster_consensus.")
  skip_on_ci()
  # because plot.lcc changes the graphics parameters, currently only the last
  # plot is recorded. Therefore, I can only check the last plot for the
  # consensus plots, even though I specified several plots.
  expect_snapshot_file(save_png(plot(clustering, which_plots = "consensusmatrix")),
                       "plot_lcc_output_consensusmatrix.png")
  expect_snapshot_file(save_png(plot(clustering, which_plots = "delta")),
                       "plot_lcc_output_delta.png")
  expect_snapshot_file(save_png(plot(clustering, which_plots = c("CDF", "cluster_tracking"))),
                       "plot_lcc_output_cdf_cluster_tracking.png")
})

test_that("the n_item_consensus argument is correct", {
  expect_error(plot(clustering, n_item_consensus = "test"))
  expect_error(plot(clustering, n_item_consensus = -1))

  skip_on_ci()
  # only the last plot (with k=3) is recorded
  expect_snapshot_file(save_png(plot(clustering, which_plots = "item_consensus",
                                     n_item_consensus = 1)),
                       "plot_lcc_output_item_consenus_1.png")
})
