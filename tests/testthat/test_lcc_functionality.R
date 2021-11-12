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
  maxK = 3,
  reps = 3,
  model_list = model_list
)

test_that("the base case works", {
  # base case, everything should work
  expect_s3_class(clustering, "lcc")
})

test_that("the output object has the correct format", {
  expect_length(clustering, 3)
  expect_named(clustering, c("general_information", "cluster_2", "cluster_3"))
})

test_that("the general_information/consensus_matrices/call object has the correct format", {
  expect_length(clustering[["general_information"]], 3)
  expect_named(clustering[["general_information"]],
               c("consensus_matrices", "cluster_assignments", "call"))
  expect_length(clustering[["general_information"]][["consensus_matrices"]], 3)
  expect_type(clustering[["general_information"]][["consensus_matrices"]][[1]],
              "NULL")
  expect_type(clustering[["general_information"]][["consensus_matrices"]][[2]],
              "double")
  expect_type(clustering[["general_information"]][["consensus_matrices"]][[3]],
              "double")
  # 15 patients
  expect_equal(dim(clustering[["general_information"]][["consensus_matrices"]][[2]]),
               c(15, 15))
  expect_type(clustering[["general_information"]][["call"]], "language")

})

test_that("the general_information/cluster_assignments object has the correct format", {
  expect_s3_class(clustering[["general_information"]][["cluster_assignments"]],
                  "data.frame")
  # 15 patients and id_column + 2 different numbers of clusters
  expect_equal(dim(clustering[["general_information"]][["cluster_assignments"]]),
               c(15, 3))
  expect_true(all(c("patient_id", "assignment_num_clus_2", "assignment_num_clus_3")
                  %in% colnames(clustering[["general_information"]][["cluster_assignments"]])))
})

test_that("the entries for every number of clusters have the correct format", {
  expect_named(clustering[[2]], c("consensusMatrix", "consensusTree",
                                  "consensusClass", "found_flexmix_clusters"))
  expect_identical(clustering[[2]][["consensusMatrix"]],
                   clustering[["general_information"]][["consensus_matrices"]][[2]])
  expect_s3_class(clustering[[2]][["consensusTree"]], "hclust")
  expect_vector(clustering[[2]][["consensusClass"]], integer(), 15)
  expect_true(max(clustering[[2]][["consensusClass"]]) <= 2)
  expect_true(max(clustering[[3]][["consensusClass"]]) <= 3)
  expect_vector(clustering[[2]][["found_flexmix_clusters"]], integer(), 3)
  expect_true(max(clustering[[2]][["found_flexmix_clusters"]]) <= 2)
  expect_true(max(clustering[[3]][["found_flexmix_clusters"]]) <= 3)
})

test_that("the cluster assignments can be retrievedby the getter", {
  cluster_assignments <- get_clusters(clustering)
  expect_s3_class(cluster_assignments, "data.frame")
  expect_equal(dim(cluster_assignments), c(15, 3))
  expect_equal(colnames(get_clusters(clustering, number_clusters = 3)),
               c("patient_id", "assignment_num_clus_3"))
  expect_error(get_clusters(data.frame(a = 1)))
  expect_error(get_clusters(clustering, number_clusters = 1))
  expect_error(get_clusters(clustering, number_clusters = 4))
})
