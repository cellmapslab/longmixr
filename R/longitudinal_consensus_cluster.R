#' Longitudinal consensus clustering with flexmix
#'
#' This function performs longitudinal clustering with flexmix. To get robust
#' results, the data is subsampled and the clustering is performed on this
#' subsample. The results are combined in a consensus matrix and a final
#' hierarchical clustering step performed on this matrix. In this, it follows
#' the approach from the \code{ConsensusClusterPlus} package.
#'
#' @param data a \code{data.frame} with one or several observations per subject.
#' It needs to contain one column that specifies to which subject the entry (row)
#' belongs to. This ID column is specified in \code{id_column}. Otherwise, there
#' are no restrictions on the column names, as the model is specified in
#' \code{flexmix_formula}.
#' @param id_column name (character vector) of the ID column in \code{data} to
#' identify all observations of one subject
#' @param max_k maximum number of clusters, default is \code{3}
#' @param reps number of repetitions, default is \code{10}
#' @param p_item fraction of samples contained in subsampled sample, default is
#' \code{0.8}
#' @param model_list either one \code{flexmix} driver or a list of \code{flexmix}
#' drivers of class \code{FLXMR}
#' @param flexmix_formula a \code{formula} object that describes the \code{flexmix}
#' model relative to the formula in the flexmix drivers (the dot in the flexmix
#' drivers is replaced, see the example). That means that you usually only
#' specify the right-hand side of the formula here. However, this is not enforced
#' or checked to give you more flexibility over the \code{flexmix} interface
#' @param title name of the clustering; used if \code{writeTable = TRUE}
#' @param final_linkage linkage used for the last hierarchical clustering step on
#' the consensus matrix; has to be \code{average, ward.D, ward.D2, single, complete, mcquitty, median}
#' or \code{centroid}. The default is \code{average}
#' @param seed seed for reproducibility
#' @param verbose \code{boolean} if status messages should be displayed.
#' Default is \code{FALSE}
#'
#' @details
#' The data types \code{longitudinal_consensus_cluster} can handle depends on
#' how the \code{flexmix} models are set up, in principle all data types are
#' supported for which there is a \code{flexmix} driver with the desired
#' outcome variable.
#'
#' If you follow the dimension reduction approach outlined in
#' \code{vignette("Example clustering analysis", package = "longmixr")}, the
#' input data types depend on what \code{FAMD} from the \code{FactoMineR}
#' package can handle. \code{FAMD} accepts \code{numeric} variables and treats
#' all other variables as \code{factor} variables which it can handle as well.
#'
#' @return An object (list) of class \code{lcc} with length \code{maxk}.
#' The first entry \code{general_information} contains the entries:\tabular{ll}{
#'    \code{consensus_matrices} \tab a list of all consensus matrices (for all specified clusters) \cr
#'    \tab \cr
#'    \code{cluster_assignments} \tab a \code{data.frame} with an ID column named after \code{id_column} and a column for every specified number of clusters, e.g. \code{assignment_num_clus_2} \cr
#'    \tab \cr
#'    \code{call} \tab the call/all arguments how \code{longitudinal_consensus_cluster} was called
#' }
#'
#' The other entries correspond to the number of specified clusters (e.g. the
#' second entry corresponds to 2 specified clusters) and each contains a list with the
#' following entries:\tabular{ll}{
#'    \code{consensus_matrix} \tab the consensus matrix \cr
#'    \tab \cr
#'    \code{consensus_tree} \tab the result of the hierarchical clustering on the consensus matrix \cr
#'    \tab \cr
#'    \code{consensus_class} \tab the resulting class for every observation \cr
#'    \tab \cr
#'    \code{found_flexmix_clusters} \tab a vector of the actual found number of clusters by \code{flexmix} (which can deviate from the specified number)
#' }
#'
#' @importFrom stats as.formula hclust as.dist cutree
#' @importFrom utils write.csv write.table
#'
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
#' # not run
#' # plot(clustering)
#' # end not run
longitudinal_consensus_cluster <- function(data = NULL,
                                           id_column = NULL,
                                           max_k = 3,
                                           reps = 10,
                                           p_item = 0.8,
                                           model_list = NULL,
                                           flexmix_formula = as.formula("~s(visit, k = 4) | patient_id"),
                                           title = "untitled_consensus_cluster",
                                           final_linkage = c("average", "ward.D", "ward.D2", "single", "complete",
                                                             "mcquitty", "median", "centroid"),
                                           seed = 3794,
                                           verbose = FALSE) {

  # check variables
  checkmate::assert_count(seed, positive = TRUE)
  set.seed(seed)

  checkmate::assert_character(id_column, len = 1, any.missing = FALSE)
  checkmate::assert_data_frame(data, all.missing = FALSE, min.rows = 1, min.cols = 1)
  checkmate::assert_choice(id_column, colnames(data))
  checkmate::assert_integerish(max_k, len = 1, lower = 2)
  checkmate::assert_count(reps, positive = TRUE)
  checkmate::assert_number(p_item, lower = 1 / nrow(data), upper = 1)
  # model_list can be either a list of flexmix drivers or a single flexmix
  # driver
  checkmate::assert(checkmate::check_list(model_list, types = "FLXMR"),
                    checkmate::check_class(model_list, "FLXMR"))
  checkmate::assert_class(flexmix_formula, "formula")
  checkmate::assert_character(title, len = 1)
  final_linkage <- match.arg(final_linkage)
  checkmate::assert_logical(verbose, len = 1)

  call <- match.call()

  # perform the longitudinal clustering to create the consensus matrices
  results <- lcc_run(data = data,
                     id_column = id_column,
                     max_k = max_k,
                     reps = reps,
                     p_item = p_item,
                     model_list = model_list,
                     flexmix_formula = flexmix_formula,
                     verbose = verbose)

  consensus_matrices <- results[["consensus_matrices"]]
  flexmix_found_clusters <- results[["found_number_clusters"]]

  # generate the consensus clustering on the consensus matrices
  # (last clustering step)
  res <- vector(mode = "list", length = max_k)
  # do the final clustering for every specified number of clusters
  for (cluster_index in 2:max_k) {
    if (verbose) {
      message(paste("consensus ", cluster_index))
    }
    consensus_matrix <- consensus_matrices[[cluster_index]]
    hc_tree <- hclust(as.dist(1 - consensus_matrix), method = final_linkage)
    if (verbose) {
      message("clustered")
    }
    cut_tree_groups <- cutree(hc_tree, cluster_index)
    names(cut_tree_groups) <- sort(unique(data[, id_column]))

    res[[cluster_index]] <- list(consensus_matrix = consensus_matrix,
                                 consensus_tree = hc_tree,
                                 consensus_class = cut_tree_groups,
                                 found_flexmix_clusters = flexmix_found_clusters[[cluster_index]])
  }

  # generate a data.frame with all cluster assignments for the subjects
  assignment_table <- extract_assignment(results = res,
                                         id_column = id_column)

  res[[1]] <- list(consensus_matrices = consensus_matrices,
                   cluster_assignments = assignment_table,
                   call = call)
  names(res)[1] <- "general_information"

  # create meaningful names for the different clustering solutions
  if (length(res) > 1) {
    names(res)[-1] <- paste0("cluster_", seq(from = 2, to = length(res), by = 1))
  }

  class(res) <- c("lcc", class(res))
  return(res)
}

#' Main function to run the longitudinal clustering
#'
#' Internal function to actually perform the clustering
#'
#' @inheritParams longitudinal_consensus_cluster
#'
#' @return Returns a list with the following entries:\tabular{ll}{
#'    \code{consensus_matrices} \tab a list of all consensus matrices where the kth entry is the consensus matrix for k specified numbers of clusters. The first entry is \code{NULL} \cr
#'    \tab \cr
#'    \code{found_number_clusters} \tab a list of vectors with the actual number of clusters found by \code{flexmix} where the kth entry is for k specified numbers of clusters. The first entry is \code{NULL}
#' }
#' @noRd
lcc_run <- function(data,
                    id_column,
                    max_k,
                    reps,
                    p_item,
                    model_list,
                    flexmix_formula,
                    verbose) {
  # internal function, therefore no input checks

  connectivity_results <- vector(mode = "list", max_k)
  # get the number of unique patients
  unique_patient_ids <- sort(unique(data[, id_column]))
  n <- length(unique_patient_ids)
  # initialise the connectivity matrix
  sample_count_matrix <- matrix(0, ncol = n, nrow = n)
  # be lazy and initialise the list with 0, because the first element corresponds
  # to a cluster size of 1
  connectivity_results[[1]] = 0
  # initialise list for found number of clusters
  found_number_clusters <- vector(mode = "list", max_k)

  for (i in 1:reps) {
    if (verbose) {
      message(paste("random subsample", i))
    }
    sample_x = sample_patients(data = data,
                               p_samp = p_item,
                               id_column = id_column)

    # save the matrix how often a combination of 2 samples occurred in this
    # subsampled data set to correct that not all samples are contained in the
    # subsample
    sample_count_matrix <- connectivity_matrix(cluster_assignments =
                                                 rep(1, length(sample_x[["sample_patients"]])),
                                               current_matrix = sample_count_matrix,
                                               matrix_names = unique_patient_ids,
                                               sample_key = sample_x[["sample_patients"]])

    # cluster the subsampled data
    fitted_models <- flexmix::stepFlexmix(formula = flexmix_formula,
                                          data = sample_x[["subsample"]],
                                          k = seq(from = 2, to = max_k, by = 1),
                                          model = model_list,
                                          nrep = 1)

    for (cluster_index in seq(from = 2, to = max_k, by = 1)) {

      if (i == 1) {
        # initialise the connectivity matrix during the first time
        connectivity_results[[cluster_index]] <- matrix(0, ncol = n, nrow = n)
      }

      # extract the cluster assignments from flexmix
      if (inherits(fitted_models, "stepFlexmix")) {
        # when more than 1 cluster values are specified, the returned value is
        # a stepFlexmix object and the models can be accessed by the cluster
        # number
        assignment_all_values <- fitted_models@models[[as.character(cluster_index)]]@cluster
        # extract the actual found number of clusters (as flexmix might find
        # a cluster solution with less clusters than specified)
        actual_number_found_clusters <- fitted_models@models[[as.character(cluster_index)]]@k

      } else if (inherits(fitted_models, "flexmix") && fitted_models@k == as.integer(cluster_index)) {
        # if only one cluster value was used during fitting (i.e. only k=2)
        assignment_all_values <- fitted_models@cluster
        actual_number_found_clusters <- fitted_models@k
      } else {
        stop("Something went wrong with extracting the cluster assignments from the flexmix models")
      }

      merged_data <- data.frame(patient_id = sample_x[["subsample"]][, id_column],
                                cluster = assignment_all_values)
      # get the cluster assignment together with the corresponding ID
      # only one value for every patient (the values should be the same for every
      # repeated measurement)
      merged_data <- merged_data[!duplicated(merged_data[, "patient_id"]), ]

      # update the connectivity matrix
      connectivity_results[[cluster_index]] <- connectivity_matrix(cluster_assignments = merged_data$cluster,
                                                                   current_matrix = connectivity_results[[cluster_index]],
                                                                   matrix_names = unique_patient_ids,
                                                                   sample_key = merged_data$patient_id)
      # store the found number of clusters
      found_number_clusters[[cluster_index]] <- c(found_number_clusters[[cluster_index]],
                                                  actual_number_found_clusters)
    }
  }

  # calculate the final consensus matrices/fraction
  res <- vector(mode = "list", max_k)
  for (cluster_index in 2:max_k) {
    tmp <- triangle(connectivity_results[[cluster_index]], mode = 3)
    tmp_count <- triangle(sample_count_matrix, mode = 3)
    # tmp is a matrix of counts how often each pair of samples is in one cluster
    # temp_count is a matrix of counts how often (maximum) each pair of samples
    # could have been in a cluster
    res[[cluster_index]] <- tmp / tmp_count
    res[[cluster_index]][which(tmp_count == 0)] <- 0
  }

  list(consensus_matrices = res,
       found_number_clusters = found_number_clusters)
}

#' Try out different linkage methods
#'
#' In the final step, the consensus clustering performs a hierarchical clustering
#' step on the consensus cluster. This function tries out different linkage
#' methods and returns the corresponding clusterings. The outputs can be plotted
#' like the results from \code{\link{longitudinal_consensus_cluster}}.
#'
#' @param results clustering result of class \code{lcc}
#' @param use_methods character vector of one or several items of \code{average},
#' \code{ward.D}, \code{ward.D2}, \code{single}, \code{complete}, \code{mcquitty},
#' \code{median} or \code{centroid}
#'
#' @return a list of elements, each element of class \code{lcc}. The entries are
#' named after the used linkage method.
#'
#' @importFrom stats hclust as.dist cutree
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
#' clustering_linkage <- test_clustering_methods(results = clustering,
#' use_methods = c("average", "single"))
#' # not run
#' # plot(clustering_linkage[["single"]])
#' # end not run
test_clustering_methods <- function(results,
                                    use_methods = c("average", "ward.D", "ward.D2", "single", "complete",
                                                    "mcquitty", "median", "centroid")) {
  checkmate::assert_class(results, "lcc")
  checkmate::assert_character(use_methods, unique = TRUE)

  # try out all specified linkage methods
  new_results <- lapply(use_methods, function(current_method) {
    # generate the consensus clustering on the consensus matrices
    curr_res <- list()
    for (tk in seq(from = 2, to = length(results), by = 1)) {

      fm <- results[[tk]][["consensus_matrix"]]
      hc <- hclust(as.dist(1 - fm), method = current_method)
      ct <- cutree(hc, tk)
      names(ct) <- names(results[[tk]][["consensus_class"]])

      curr_res[[tk]] <- list(consensus_matrix = fm,
                             consensus_tree = hc,
                             consensus_class = ct,
                             found_flexmix_clusters = results[[tk]][["found_flexmix_clusters"]])
    }

    # add the correct general information
    original_assignments <- results[["general_information"]][["cluster_assignments"]]
    id_column <- colnames(original_assignments)[1]
    assignment_table <- extract_assignment(results = curr_res,
                                           id_column = id_column)

    # gather all consensus matrices to one list
    consensus_matrices <- list()
    for (i in seq(from = 2, to = length(results), by = 1)) {
      consensus_matrices[[i]] <- curr_res[[i]][["consensus_matrix"]]
    }
    curr_res[[1]] <- list(consensus_matrices = consensus_matrices,
                          cluster_assignments = assignment_table)
    names(curr_res)[1] <- "general_information"

    class(curr_res) <- c("lcc", class(curr_res))
    curr_res
  })

  names(new_results) <- use_methods
  new_results
}
