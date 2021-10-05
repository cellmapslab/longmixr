set.seed(5)
test_data <- data.frame(patient_id = rep(1:10, each = 4),
                        visit = rep(1:4, 10),
                        var_1 = c(rnorm(20, -1), rnorm(20, 3)) +
                          rep(seq(from = 0, to = 1.5, length.out = 4), 10),
                        var_2 = c(rnorm(20, 0.5, 1.5), rnorm(20, -2, 0.3)) +
                          rep(seq(from = 1.5, to = 0, length.out = 4), 10))
model_list <- list(flexmix::FLXMRmgcv(as.formula("var_1 ~ .")),
                   flexmix::FLXMRmgcv(as.formula("var_2 ~ .")))

clustering <- longitudinal_consensus_cluster(
  data = test_data,
  id_column = "patient_id",
  max_k = 2,
  reps = 3,
  model_list = model_list
)

test_that("the base case works", {
  expect_type(test_clustering_methods(clustering, c("average", "ward.D")),
                  "list")

  expect_length(test_clustering_methods(clustering, c("average", "ward.D")),
                2)

  expect_s3_class(test_clustering_methods(clustering, c("average", "ward.D"))[[1]],
                  "lcc")

  expect_s3_class(test_clustering_methods(clustering, c("average", "ward.D"))[[2]],
                  "lcc")
})

test_that("the input arguments are correct", {
  expect_error(test_clustering_methods(1:10, c("average", "ward.D")),
               regexp = "Assertion on 'results' failed")

  expect_error(test_clustering_methods(clustering, c("superman", "superwoman")),
               regexp = "invalid clustering method")

  expect_error(test_clustering_methods(clustering, 1:4),
               regexp = "Assertion on 'use_methods' failed")
})
