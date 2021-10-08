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
save_png <- function(code, width = 400, height = 3200) {
  path <- tempfile(fileext = ".png")
  png(path, width = width, height = height)
  op <- par(mfrow = c(8, 1))
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
