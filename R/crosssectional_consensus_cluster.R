#' Cross-sectional clustering with categorical variables
#'
#' This function uses the \code{ConsensusClusterPlus} function from the package
#' with the same name with defaults for clustering data with categorical
#' variables.
#'
#' @param data a matrix or data.frame containing variables that should be used
#' for computing the distance. This argument is passed to \code{StatMatch::gower.dist}
#' @param reps number of repetitions, same as in \code{ConsensusClusterPlus}
#' @param finalLinkage linkage method for final clustering,
#' same as in \code{ConsensusClusterPlus}same as in \code{ConsensusClusterPlus}
#' @param innerLinkage linkage method for clustering steps,
#' same as in \code{ConsensusClusterPlus}
#' @param ... other arguments passed to \code{ConsensusClusterPlus}, attention:
#' the \code{d} argument can \bold{not} be set as it is directly computed by
#' \code{crosssectional_consensus_cluster}
#'
#' @return The output is produced by \code{ConsensusClusterPlus}
#'
#' @importFrom stats as.dist
#' @export
#'
#' @examples
#' dc <- mtcars
#' # scale continuous variables
#' dc <- sapply(mtcars[, 1:7], scale)
#' # code factor variables
#' dc <- cbind(as.data.frame(dc),
#'             vs = as.factor(mtcars$vs),
#'             am = as.factor(mtcars$am),
#'             gear = as.factor(mtcars$gear),
#'             carb = as.factor(mtcars$carb))
#' cc <- crosssectional_consensus_cluster(
#'   data = dc,
#'   reps = 10,
#'   seed = 1
#' )
crosssectional_consensus_cluster <- function(
  data,
  reps = 1000,
  finalLinkage = "ward.D2",
  innerLinkage = "ward.D2",
  ...) {

  data_as_distance <- as.dist(StatMatch::gower.dist(data))
  consensus_cluster <- ConsensusClusterPlus::ConsensusClusterPlus(
    d = data_as_distance,
    reps = reps,
    # this argument is actually not needed as it doesn't do anything when d is
    # a distance object (as it is the case here)
    distance = "StatMatch::gower.dist",
    finalLinkage = finalLinkage,
    innerLinkage = innerLinkage,
    ...)

  consensus_cluster
}
