dc <- mtcars
# scale continuous variables
dc <- sapply(mtcars[, 1:7], scale)
# code factor variables
dc <- cbind(as.data.frame(dc),
            vs = as.factor(mtcars$vs),
            am = as.factor(mtcars$am),
            gear = as.factor(mtcars$gear),
            carb = as.factor(mtcars$carb))

# as here I basically only use functions from other packages, I only test if
# the d argument is set correctly

# set up a function to run the function that calls ConsensusClusterPlus
# to manually turn off any device or changed plotting codes
# set pdf(NULL) because otherwise the test leaves a Rplots.pdf file
run_ccc <- function(...) {
  pdf(NULL)
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  on.exit(dev.off(), add = TRUE)

  crosssectional_consensus_cluster(...)
}

test_that("the d argument is set correctly", {
  # this base case should work
  expect_type(run_ccc(
    data = dc,
    reps = 3,
    seed = 1
  ),
  "list")

  # no extra d argument is allowed
  expect_error(run_ccc(
    data = dc,
    d = dc,
    reps = 3,
    seed = 1
  ))
})
