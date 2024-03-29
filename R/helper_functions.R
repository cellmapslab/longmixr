#' Subsample subjects
#'
#' Internal function to subsample a fraction of \code{pSamp} subjects and
#' return all observations belonging to these subjects. Adapted from
#' \code{ConsensusClusterPlus}.
#'
#' @inheritParams longitudinal_consensus_cluster
#' @param p_samp subsampling fraction, same as \code{p_item}
#'
#' @return List with the following entries:\tabular{ll}{
#'    \code{subsample} \tab \code{data.frame} of the subsample \cr
#'    \tab \cr
#'    \code{sample_patients} \tab names from the \code{id_column} of the subsampled subjects
#' }
#' @noRd
sample_patients <- function(data,
                            p_samp,
                            id_column) {
  # as there are repeated measurements in the data, sample patients and then
  # take all measurements from the sampled patients
  space <- unique(data[, id_column])
  space_dim <- length(space)
  sample_n <- floor(space_dim * p_samp)
  sample_patients <- sort(sample(space, sample_n, replace = FALSE))
  index <- data[, id_column] %in% sample_patients
  this_sample <- data[index, ]

  list(subsample = this_sample, sample_patients = sample_patients)
}

#' Update connectivity matrix
#'
#' Internal function to update a matrix that stores information how often two
#' samples have gotten the same cluster assignment.
#' Adapted from \code{ConsensusClusterPlus}.
#'
#' @param cluster_assignments vector of cluster assignment for every subject
#' @param current_matrix current connectivity matrix
#' @param matrix_names vector of the names of the rows/columns of the
#' connectivity matrix; usually the sorted unique names of the subjects
#' @param sample_key name of the subjects
#'
#' @return updated connectivity matrix
#' @noRd
connectivity_matrix <- function(cluster_assignments,
                                current_matrix,
                                matrix_names,
                                sample_key) {

  names(cluster_assignments) <- sample_key
  # for every cluster, get the list of patient ids that belong to this cluster
  cls <- lapply(unique(cluster_assignments), function(i) {
    names(cluster_assignments[cluster_assignments %in% i])
  })
  for (i in 1:length(cls)) {
    # check which rows/columns of the matrix (specified by matrix_names)
    # belong to the same cluster
    cl <- as.numeric(matrix_names %in% cls[[i]])
    # product of arrays with * function
    # with the 1/0 indicator (which sample were observed together), it updates
    # all cells to indicate the sample pair was observed or not
    updt <- outer(cl, cl)
    current_matrix <- current_matrix + updt
  }
  current_matrix
}

#' Generate triangle matrix
#'
#' Internal function to generate a triangle matrix;
#' adapted from \code{ConsensusClusterPlus}.
#'
#' @param input_matrix matrix
#' @param mode flag which matrix should be returned; can be \code{1, 2} or \code{3}
#'
#' @return matrix; if \code{mode = 1} return the lower triangle as vector;
#' if \code{mode = 2} return the transformed upper triangle (so that it is the lower one now)
#' as a matrix with the rest of the entries as \code{NA}; if \code{mode = 3}
#' return a matrix where the lower left triangle is replaced by the upper right triangle
#' @noRd
triangle <- function(input_matrix,
                     mode = 1) {
  n <- dim(input_matrix)[1]
  nm <- matrix(0, ncol = n, nrow = n)
  fm <- input_matrix

  # only use the upper half
  nm[upper.tri(nm)] <- input_matrix[upper.tri(input_matrix)]
  # in fm, the lower half is the same as the upper half
  fm <- t(nm) + nm
  diag(fm) <- diag(input_matrix)

  # after the above commands, fm is now a matrix where the lower half is the
  # same as the upper half and the diagonal is taken from the original matrix
  # I'm not sure why the original authors did it so complicated as the
  # connectivity matrices used as input should be symmetrical and therefore
  # fm should now be in the same format as input_matrix

  # generate a matrix (nm) where only the lower half is left over
  nm <- fm
  nm[upper.tri(nm)] <- NA
  diag(nm) <- NA
  vm <- input_matrix[lower.tri(nm)]

  if (mode == 1) {
    return(vm)
  }
  else if (mode == 3) {
    return(fm)
  }
  else if (mode == 2) {
    return(nm)
  }
}

#' Plot the CDFs of consensus clustering
#'
#' Internal function to plot the CDFs of the consensus solutions;
#' adapted from \code{ConsensusClusterPlus}.
#'
#' @param matrix_list list of all consensus matrices
#' @param breaks number of breaks
#' @param which_plots which plots should be plotted, can include \code{"CDF"} or
#' \code{"delta"}
#'
#' @importFrom graphics hist lines legend
#' @importFrom grDevices rainbow
#'
#' @return a CDF plot and/or delta CDF plot
#' @noRd
CDF <- function(matrix_list,
                breaks = 100,
                which_plots = c("CDF", "delta")) {
  if ("CDF" %in% which_plots) {
    # set up the plot
    plot(0, xlim = c(0, 1), ylim = c(0, 1), col = "white",
         bg = "white", xlab = "consensus index", ylab = "CDF",
         main = "consensus CDF", las = 2)
  }

  k <- length(matrix_list)
  this_colors <- rainbow(k - 1)
  area_k <- c()
  for (i in 2:length(matrix_list)) {
    # get the lower triangle of the connectivity matrix as a vector
    v <- triangle(matrix_list[[i]], mode = 1)

    # empirical CDF distribution
    h <- hist(v, plot = FALSE, breaks = seq(0, 1, by = 1 / breaks))
    h$counts <- cumsum(h$counts) / sum(h$counts)

    # calculate the area under the CDF curve by the histogram method
    this_area <- 0
    for (bi in 1:(length(h$breaks) - 1)) {
      this_area <- this_area + h$counts[bi] * (h$breaks[bi + 1] - h$breaks[bi])
    }
    area_k <- c(area_k, this_area)
    if ("CDF" %in% which_plots) {
      # add the CDF to the plot
      lines(h$mids, h$counts, col = this_colors[i - 1], lwd = 2,
            type = "l")
    }
  }
  if ("CDF" %in% which_plots) {
    legend(0.8, 0.5, legend = paste(rep("", k - 1), seq(2, k, by = 1), sep = ""),
           fill = this_colors)
  }

  # plot the area under the CDF change
  delta_k <- area_k[1]
  for (i in 2:(length(area_k))) {
    # proportional increase relative to the previous k
    delta_k <- c(delta_k, (area_k[i] - area_k[i - 1]) / area_k[i - 1])
  }
  if ("delta" %in% which_plots) {
    plot(1 + (1:length(delta_k)), y = delta_k, xlab = "k",
         ylab = "relative change in area under CDF curve",
         main = "Delta area", type = "b")
  }
}

#' Assign colours to cluster assignments
#'
#' Internal function to assign colours to the cluster assignments;
#' adapted from \code{ConsensusClusterPlus}.
#'
#' @param past_ct cluster assignments of the previous clustering (\code{k-1} numbers of clusters)
#' @param ct current cluster assignments (\code{k} numbers of clusters)
#' @param color_names vector of colour names
#' @param color_list results of a previous call to \code{set_cluster_colors }
#'
#' @return List with the following entries:\tabular{ll}{
#'    \code{1.} \tab vector of colours for every subject\cr
#'    \tab \cr
#'    \code{2.} \tab number; probably the number of different colours used (I haven't checked it)\cr
#'    \tab \cr
#'    \code{3.} \tab vector of unique colours in the first list entry
#' }
#' @noRd
set_cluster_colors <- function(past_ct,
                               ct,
                               color_names,
                               color_list) {

  new_colors <- c()
  if (length(color_list) == 0) {
    new_colors <- color_names[ct]
    color_i <- 2
  }
  else {
    new_colors <- rep(NULL, length(ct))
    color_i <- color_list[[2]]
    mo <- table(past_ct, ct)
    m <- mo / apply(mo, 1, sum)
    # do the following for each cluster
    for (tci in 1:ncol(m)) {
      max_c <- max(m[, tci])
      pci <- which(m[, tci] == max_c)
      # if the new column maximum is unique, the same cell is the row maximum
      # and is also unique
      if (sum(m[, tci] == max_c) == 1 &
          max(m[pci, ]) == max_c &
          sum(m[pci, ] == max_c) == 1) {
        # note: the greatest of the previous clusters' members are the greatest
        # in a current cluster's members
        new_colors[which(ct == tci)] <-
          unique(color_list[[1]][which(past_ct == pci)])
      }
      else {
        color_i <- color_i + 1
        new_colors[which(ct == tci)] <- color_names[color_i]
      }
    }
  }
  return(list(new_colors, color_i, unique(new_colors)))
}

#' Generate the tracking plot
#'
#' Internal function to generate the tracking plot;
#' adapted from \code{ConsensusClusterPlus}.
#'
#' @param assignment_matrix matrix of assigned colours (which here are equal
#' with cluster assignment) to every observation with one entry per specified
#' number of cluster (i.e. the rows are the different k and the columns the samples)
#'
#' @importFrom graphics rect segments text
#'
#' @return tracking plot
#' @noRd
cluster_tracking_plot <- function(assignment_matrix) {
  # set up the plot
  plot(NULL, xlim = c(-0.1, 1), ylim = c(0, 1), axes = FALSE,
       xlab = "samples", ylab = "k", main = "tracking plot")

  for (i in 1:nrow(assignment_matrix)) {
    rect(xleft = seq(0, 1 - 1 / ncol(assignment_matrix),
                     by = 1 / ncol(assignment_matrix)),
         ybottom = rep(1 - i / nrow(assignment_matrix), ncol(assignment_matrix)),
         xright = seq(1 / ncol(assignment_matrix), 1,
                      by = 1 / ncol(assignment_matrix)),
         ytop = rep(1 - (i - 1) / nrow(assignment_matrix),
                    ncol(assignment_matrix)),
         col = assignment_matrix[i, ], border = NA)
  }
  # hatch lines to indicate samples
  xl <- seq(0, 1 - 1 / ncol(assignment_matrix), by = 1 / ncol(assignment_matrix))
  segments(xl, rep(-0.1, ncol(assignment_matrix)), xl,
           rep(0, ncol(assignment_matrix)), col = "black")
  ypos = seq(1, 0, by = -1 / nrow(assignment_matrix)) -
    1 / (2 * nrow(assignment_matrix))
  text(x = -0.1, y = ypos[-length(ypos)],
       labels = seq(2, nrow(assignment_matrix) + 1, by = 1))
}

#' Generate colour palette
#'
#' Internal function to generate colour palette;
#' adapted from \code{ConsensusClusterPlus}.
#'
#' @param n number of colours
#'
#' @importFrom grDevices rgb
#'
#' @return vector of colour hex values where the first is red and the rest different
#' blue colours
#' @noRd
my_pal <- function(n = 10) {
  seq <- rev(seq(0, 255, by = 255 / (n)))
  pal_RGB <- cbind(seq, seq, 255)
  rgb(pal_RGB, maxColorValue = 255)
}

#' Extract one overview assignment table
#'
#' This internal function generates a \code{data.frame} with the cluster
#' assignments for every subject across all specified numbers of clusters
#'
#' @param results result list where the first entry is ignored and the other
#' entries respond to the number of specified clusters (the kth entry means
#' k specified clusters); every entry needs to contain the entry
#' \code{consensus_class}. The cluster assignments need to have the subject
#' name as attribute
#' @param id_column character vector of the ID column in the dataset
#'
#' @return \code{data.frame} with a column named after the \code{id_column}
#' containing the subject name and one column with the cluster assignments for
#' every specified number of clusters, e.g. \code{assignment_num_clus_2}
#' @noRd
extract_assignment <- function(results,
                               id_column) {
  # extract the assignment for every specified number of clusters
  cluster_assignments <- lapply(seq(from = 2, to = length(results), by = 1),
                                function(i) {
                                  results[[i]][["consensus_class"]]
                                })
  # generate a data.frame; the rownames are the subject names
  merged_data <- do.call("cbind", cluster_assignments)
  patient_id <- rownames(merged_data)
  merged_data <- as.data.frame(merged_data)
  rownames(merged_data) <- NULL
  merged_data <- cbind(patient_id, merged_data)
  colnames(merged_data) <- c(id_column, paste0("assignment_num_clus_",
                                               seq(from = 2, to = length(results),
                                                   by = 1)))
  merged_data

}

#' Generate item/subject consensus plot
#'
#' Internal function to generate the item consensus plot;
#' adapted from \code{ConsensusClusterPlus}.
#'
#' @param item_consensus_matrix matrix where the rows are different clusters,
#' the columns the samples and the value the item consensus of the sample with
#' the different clusters
#' @param cluster_colors colours for the clusters
#' @param item_order order of the items/subjects as in the solution for 2 clusters
#' @param title title of the plot
#'
#' @return item/subject consensus plot
#' @noRd
ranked_bar_plot <- function(item_consensus_matrix,
                            cluster_colors,
                            item_order,
                            title) {
  colors <- rbind()
  by_rank <- cbind()
  # space between bars
  spaceh <- 0.1
  # for every item, bring the data into the correct format
  for (i in 1:ncol(item_consensus_matrix)) {
    by_rank <- cbind(by_rank, sort(item_consensus_matrix[, i], na.last = FALSE))
    colors <- rbind(colors, order(item_consensus_matrix[, i], na.last = FALSE))
  }
  # maximum height of the graph
  max_h <- max(c(1.5, apply(by_rank, 2, sum)), na.rm = TRUE)
  # barplot largest to smallest so that the smallest is in front
  barplot(apply(by_rank, 2, sum), col = cluster_colors[colors[, 1]],
          space = spaceh, ylim = c(0, max_h),
          main = paste("item-consensus", title), border = NA, las = 1)
  for (i in 2:nrow(by_rank)) {
    barplot(apply(matrix(by_rank[i:nrow(by_rank), ], ncol = ncol(by_rank)),
                  2, sum), space = spaceh, col = cluster_colors[colors[, i]],
            ylim = c(0, max_h), add = TRUE, border = NA, las = 1)
  }
  xr <- seq(spaceh, ncol(item_consensus_matrix) +
              ncol(item_consensus_matrix) * spaceh,
            (ncol(item_consensus_matrix) + ncol(item_consensus_matrix) * spaceh) /
              ncol(item_consensus_matrix))
  text("*", x = xr + 0.5, y = max_h, col = cluster_colors[item_order], cex = 1.4)
}

#' Extract the cluster assignments
#'
#' This functions extracts the cluster assignments from an \code{lcc} object.
#' One can specify which for which number of clusters the assignments
#' should be returned.
#'
#' @param cluster_solution an \code{lcc} object
#' @param number_clusters default is \code{NULL} to return all assignments.
#' Otherwise specify a numeric vector with the number of clusters for which the
#' assignments should be returned, e.g. \code{2:4}
#'
#' @return a \code{data.frame} with an ID column (the name of the ID column
#' was specified by the user when calling the
#' \code{longitudinal_consensus_cluster}) function and one column with cluster
#' assignments for every specified number of clusters. Only the assignments
#' included in \code{number_clusters} are returned in the form of columns with
#' the names \code{assignment_num_clus_x}
#' @export
#'
#' @examples
#' # not run
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
#' cluster_assignments <- get_clusters(clustering, number_clusters = 2)
#' # end not run
get_clusters <- function(cluster_solution,
                         number_clusters = NULL) {
  checkmate::assert_class(cluster_solution, "lcc")
  cluster_assignments <-
    cluster_solution$general_information$cluster_assignments
  checkmate::assert_numeric(number_clusters,
                            lower = 2,
                            upper = ncol(cluster_assignments),
                            null.ok = TRUE,
                            any.missing = FALSE,
                            all.missing = FALSE,
                            unique = TRUE)

  if (is.null(number_clusters)) {
    cluster_assignments
  } else {
    # this assumes that the columns are ordered with increasing number of
    # clusters and that the first column is the ID column (with the ID column
    # name provided by the user to the longitudinal_consensus_cluster function)
    # as it is currently the case
    cluster_assignments[, c(1, number_clusters)]
  }
}
