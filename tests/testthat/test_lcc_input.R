set.seed(5)
test_data <- data.frame(patient_id = rep(1:10, each = 4),
                        visit = rep(1:4, 10),
                        var_1 = c(rnorm(20, -1), rnorm(20, 3)) +
                          rep(seq(from = 0, to = 1.5, length.out = 4), 10),
                        var_2 = c(rnorm(20, 0.5, 1.5), rnorm(20, -2, 0.3)) +
                          rep(seq(from = 1.5, to = 0, length.out = 4), 10))
model_list <- list(flexmix::FLXMRmgcv(as.formula("var_1 ~ .")),
                   flexmix::FLXMRmgcv(as.formula("var_2 ~ .")))

test_that("the data argument is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = as.matrix(test_data),
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'data' failed")

  test_data_2 <- test_data
  # change the ID column name
  colnames(test_data_2)[1] <- "some_other_id"
  expect_error(longitudinal_consensus_cluster(
    data = test_data_2,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'id_column' failed")

  expect_error(longitudinal_consensus_cluster(
    data = c(test_data, test_data),
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'data' failed")
})

test_that("the id_column argument is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "some_other_id",
    maxK = 2,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'id_column' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = 3,
    maxK = 2,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'id_column' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = c("patient_id", "patient_id"),
    maxK = 2,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'id_column' failed")
})

test_that("the maxK argument is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 1,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'maxK' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = -3,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'maxK' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2.5,
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'maxK' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = "3",
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'maxK' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = c(2, 3),
    reps = 3,
    model_list = model_list
  ),
  regexp = "Assertion on 'maxK' failed")
})

test_that("the reps argument is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3.5,
    model_list = model_list
  ),
  regexp = "Assertion on 'reps' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 0,
    model_list = model_list
  ),
  regexp = "Assertion on 'reps' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = -1,
    model_list = model_list
  ),
  regexp = "Assertion on 'reps' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = "3",
    model_list = model_list
  ),
  regexp = "Assertion on 'reps' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = c(3, 4),
    model_list = model_list
  ),
  regexp = "Assertion on 'reps' failed")
})

test_that("the pItem argument is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 1.1,
    model_list = model_list
  ),
  regexp = "Assertion on 'pItem' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = "0.8",
    model_list = model_list
  ),
  regexp = "Assertion on 'pItem' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = -0.1,
    model_list = model_list
  ),
  regexp = "Assertion on 'pItem' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 0,
    model_list = model_list
  ),
  regexp = "Assertion on 'pItem' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = c(0.5, 0.8),
    model_list = model_list
  ),
  regexp = "Assertion on 'pItem' failed")
})

test_that("the model_list argument is correct", {
  # in the base example, I use a list of flexmix drivers
  # now test that it also works with just one flexmix driver
  expect_s3_class(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 0.8,
    model_list = flexmix::FLXMRmgcv(as.formula("var_1 ~ ."))
  ),
  "lcc")

  expect_s3_class(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 0.8,
    model_list = list(flexmix::FLXMRmgcv(as.formula("var_1 ~ .")))
  ),
  "lcc")

  # the driver has to have a FLXMR class
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 0.8,
    model_list = flexmix::FLXMCmvnorm(as.formula("var_1 ~ ."), diagonal = TRUE)
  )
  )

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 0.8,
    model_list = list(flexmix::FLXMRmgcv(as.formula("var_1 ~ .")),
                      flexmix::FLXMCmvnorm(as.formula("var_2 ~ ."), diagonal = TRUE))
  )
  )
})

test_that("the flexmix_formula argument is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 0.8,
    model_list = model_list,
    flexmix_formula = "~s(visit, k = 4) | patient_id"
  ),
  regexp = "Assertion on 'flexmix_formula' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 0.8,
    model_list = model_list,
    flexmix_formula = c(as.formula("~s(visit, k = 4) | patient_id"),
                        as.formula("~s(visit, k = 4) | patient_id"))
  ),
  regexp = "Assertion on 'flexmix_formula' failed")

  # the formula musn't contain anything that is not specified in the data
  # argument
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    pItem = 0.8,
    model_list = model_list,
    flexmix_formula = as.formula("~s(studyday, k = 4) | patient_id")
  )
  )
})

test_that("the title argument is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    title = 345
  ),
  regexp = "Assertion on 'title' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    title = c("title1", "title2")
  ),
  regexp = "Assertion on 'title' failed")
})

test_that("the argument finalLinkage is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    finalLinkage = "superman"
  ),
  regexp = "'arg' should be one of")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    finalLinkage = 3
  ),
  regexp = "'arg' must be NULL or a character vector")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    finalLinkage = c("single", "complete")
  ),
  regexp = "'arg' must be of length 1")
})

test_that("the argument seed is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    seed = 4.5
  ),
  regexp = "Assertion on 'seed' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    seed = -1
  ),
  regexp = "Assertion on 'seed' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    seed = "3"
  ),
  regexp = "Assertion on 'seed' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    seed = c(3, 4)
  ),
  regexp = "Assertion on 'seed' failed")
})

test_that("the argument verbose is correct", {
  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    verbose = "TRUE"
  ),
  regexp = "Assertion on 'verbose' failed")

  expect_error(longitudinal_consensus_cluster(
    data = test_data,
    id_column = "patient_id",
    maxK = 2,
    reps = 3,
    model_list = model_list,
    verbose = c(TRUE, FALSE)
  ),
  regexp = "Assertion on 'verbose' failed")
})
