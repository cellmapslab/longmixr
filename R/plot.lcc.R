#' Plot a longitudinal consensus clustering
#'
#' @param x \code{lcc} object (output from \code{\link{longitudinal_consensus_cluster}})
#' @param color_palette optional character vector of colors for consensus matrix
#' @param which_plots determine which plots should be plotted; the default is \code{"all"}.
#' Alternatively, a combination of the following values can be specified to plot
#' only some of the below mentioned plots: \code{"consensusmatrix_legend"},
#' \code{"consensusmatrix_x"} where \code{x} is replaced by the corresponding number
#' of clusters, \code{"CDF"}, \code{"delta"}, \code{"cluster_tracking"},
#' \code{"item_consensus"} or \code{"cluster_consensus"}. When you want to plot
#' all consensus matrices and the legend, you can just use \code{"consensusmatrix"}.
#' @param ... additional parameters for plotting; currently not used
#'
#' @return Plots the following plots (when selected):\tabular{ll}{
#'    \code{consensus matrix legend} \tab the legend for the following consensus matrix plots (select with \code{"consensusmatrix_legend"}) \cr
#'    \tab \cr
#'    \code{consensus matrix plot} \tab for every specified number of clusters, a heatmap of the consensus matrix and the result of the final clustering is shown (select with \code{"consensusmatrix_x"} where \code{x} is replaced by the corresponding number
#' of clusters) \cr
#'    \tab \cr
#'    \code{consensus CDF} \tab a line plot of the CDFs for all different specified numbers of clusters (select with \code{"CDF"})\cr
#'    \tab \cr
#'    \code{Delta area} \tab elbow plot of the difference in the CDFs between the different numbers of clusters (select with \code{"delta"}) \cr
#'    \tab \cr
#'    \code{tracking plot} \tab cluster assignment of the subjects throughout the different cluster solutions (select with \code{"cluster_tracking"}) \cr
#'    \tab \cr
#'    \code{item-consensus} \tab for every item (subject), calculate the average consensus value with all items that are assigned to one consensus cluster. This is repeated for every cluster and for all different numbers of clusters (select with \code{"item_consensus"}) \cr
#'    \tab \cr
#'    \code{cluster-consensus} \tab every bar represents the average pair-wise item-consensus within one consensus cluster (select with \code{"cluster_consensus"})
#' }
#'
#' @importFrom stats as.dendrogram heatmap median
#' @importFrom graphics barplot par
#'
#' @export
plot.lcc <- function(x,
                     color_palette = NULL,
                     which_plots = "all",
                     ...) {

  checkmate::assert_class(x, "lcc")
  checkmate::assert_character(color_palette, null.ok = TRUE)

  # determine the possible consensus matrices
  possible_consensusmatrix <- paste0("consensusmatrix_",
                                     seq(from = 2, to = length(x), by = 1))
  possible_plots <- c("all", "consensusmatrix_legend", "consensusmatrix",
                      possible_consensusmatrix, "CDF", "delta",
                      "cluster_tracking", "item_consensus", "cluster_consensus")

  checkmate::assert_character(which_plots)
  if (!all(which_plots %in% possible_plots)) {
    stop(paste0("which_plot must be one of ",
                paste0(possible_plots, collapse = ", "), "."))
  }

  # if "all" was selected for the plots, set which_plot to all possible plots,
  # so that I don't need to include "all" in every if condition when checking
  # which plots to plot
  if ("all" %in% which_plots) {
    which_plots <- possible_plots
  }

  # set up the colour palette
  color_list <- list()
  color_matrix <- NULL
  this_pal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
                "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca",
                "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106",
                "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776",
                "#ffffff")

  # set up the plot scale
  col_breaks <- NA
  if (is.null(color_palette)) {
    col_breaks <- 10
    color_palette <- my_pal(col_breaks)
  }
  else {
    col_breaks <- length(color_palette)
  }

  ##############################################################################
  # plot the consensus matrices
  ##############################################################################

  # plot the legend
  sc <- cbind(seq(0, 1, by = 1 / col_breaks))
  rownames(sc) <- sc[, 1]
  sc <- cbind(sc, sc)
  if ("consensusmatrix_legend" %in% which_plots ||
      "consensusmatrix" %in% which_plots) {
    heatmap(sc,
            Colv = NA,
            Rowv = NA,
            symm = FALSE,
            scale = "none",
            col = color_palette,
            na.rm = TRUE,
            labRow = rownames(sc),
            labCol = FALSE,
            main = "consensus matrix legend")
  }


  # plot the consensus matrices for every number of clusters
  # for every cluster, calculate the correct colours for every observation
  for (tk in seq(from = 2, to = length(x), by = 1)) {

    c_matrix <- x[[tk]][["consensus_matrix"]]
    c_tree <- x[[tk]][["consensus_tree"]]
    c_class <- x[[tk]][["consensus_class"]]

    found_flexmix_clusters <- x[[tk]][["found_flexmix_clusters"]]
    median_found_flexmix_clusters <- median(found_flexmix_clusters)

    # for every cluster solution except the first define the previous consensus
    # class so that the colours are assigned correctly across plots
    if (tk == 2) {
      previous_c_class <- NULL
    } else {
      previous_c_class <- x[[tk - 1]][["consensus_class"]]
    }
    color_list <- set_cluster_colors(previous_c_class,
                                     c_class,
                                     this_pal,
                                     color_list)

    # row ordered matrix for plotting with additional row of 0s (as in the
    # original ConsensusClusterPlus code)
    plot_c_matrix <- rbind(c_matrix[c_tree$order, ], 0)

    if (paste0("consensusmatrix_", tk) %in% which_plots ||
        "consensusmatrix" %in% which_plots) {
      heatmap(plot_c_matrix,
              Colv = as.dendrogram(c_tree),
              Rowv = NA,
              symm = FALSE,
              scale = "none",
              col = color_palette,
              na.rm = TRUE,
              labRow = FALSE,
              labCol = FALSE,
              margins = c(5, 5),
              main = paste("consensus matrix k=", tk, "; median flexmix clusters: ",
                           median_found_flexmix_clusters, sep = ""),
              ColSideColors = color_list[[1]])
      legend("topright", legend = unique(c_class), fill = unique(color_list[[1]]),
             horiz = FALSE)
    }
    color_matrix <- rbind(color_matrix, color_list[[1]])
  }

  ##############################################################################
  # plot the CDF, delta CDF and observation tracking plots
  ##############################################################################
  if ("CDF" %in% which_plots || "delta" %in% which_plots) {
    CDF(x[["general_information"]][["consensus_matrices"]],
        which_plots = which_plots)
  }

  n_last_element <- length(x)
  colour_tracking_matrix <- color_matrix[, x[[n_last_element]]$consensus_tree$order]
  # if only 2 clusters were specified, colour_tracking_matrix is not a matrix
  # but a vector -> then the plot doesn't work -> transform
  if (!is.matrix(colour_tracking_matrix) && n_last_element == 2) {
    colour_tracking_matrix <- matrix(colour_tracking_matrix, nrow = 1)
  }
  if ("cluster_tracking" %in% which_plots) {
    cluster_tracking_plot(colour_tracking_matrix)
  }

  ##############################################################################
  # plot the item-consensus
  ##############################################################################

  cluster_consensus <- rbind()
  cci <- rbind()
  sumx <- list()
  colors_arr <- c()
  old_par <- par(mfrow = c(3, 1), mar = c(4, 3, 2, 0))
  on.exit(par(old_par))
  # tk is the number of predefined clusters
  for (tk in seq(from = 2, to = length(x), by = 1)) {
    ei_cols <- c()
    c_matrix <- x[[tk]][["consensus_matrix"]]
    c_class <- x[[tk]][["consensus_class"]]

    # for every subject (item), calculate the average consensus value with all
    # subjects who are grouped into one cluster
    # do this for every cluster
    c_matrix <- triangle(c_matrix, mode = 2)
    # for each cluster in tk/the predefined number of clusters
    # e.g. for tk = 2, there should be 2 clusters and for both clusters the
    # mean item consensus is calculated
    for (cluster_i in sort(unique(c_class))) {
      items <- which(c_class == cluster_i)
      n_k <- length(items)
      mk <- sum(c_matrix[items, items], na.rm = TRUE) / ((n_k * (n_k - 1)) / 2)
      # cluster consensus
      cluster_consensus <- rbind(cluster_consensus, c(tk, cluster_i, mk))
      for (item_i in rev(x[[2]]$consensus_tree$order)) {
        denom <- if (item_i %in% items) {
          n_k - 1
        }
        else {
          n_k
        }
        # mean item consensus to a cluster
        mean_item_consensus <- sum(c(c_matrix[item_i, items],
                                     c_matrix[items, item_i]),
                                   na.rm = TRUE) / denom
        # add a new row with cluster, cluster index, item index, item consensus
        cci <- rbind(cci, c(tk, cluster_i, item_i, mean_item_consensus))
      }
      ei_cols <- c(ei_cols, rep(cluster_i, length(c_class)))
    }
    # only plot the new tk data
    cck <- cci[which(cci[, 1] == tk), ]
    # group by item, order by cluster i
    w <- lapply(split(cck, cck[, 3]), function(x) {
      y <- matrix(unlist(x), ncol = 4)
      y[order(y[, 2]), 4]
    })

    # set up the matrix for plotting
    q <- matrix(as.numeric(unlist(w)), ncol = length(w), byrow = FALSE)
    # order by leave order of tk = 2
    q <- q[, x[[2]]$consensus_tree$order]
    # this results in q: a matrix of tk rows and sample columns, values are
    # item consensus of sample to the cluster
    # so for a defined possible number of clusters (tk), the values in the rows
    # are the item consensus for the possible clusters

    # it needs to be colorM[tk - 1, ] because the first element in
    # colorM refers to tk (so for 2 clusters, the information is stored in the
    # first entry and not in the second)
    this_colors <- unique(cbind(x[[tk]]$consensus_class, color_matrix[tk - 1, ]))
    this_colors <- this_colors[order(as.numeric(this_colors[, 1])), 2]
    colors_arr <- c(colors_arr, this_colors)
    if ("item_consensus" %in% which_plots) {
      ranked_bar_plot(item_consensus_matrix = q,
                      cluster_colors = this_colors,
                      item_order = c_class[x[[2]]$consensus_tree$order],
                      title = paste("k=", tk, sep = ""))
    }
  }

  ##############################################################################
  # plot the cluster-consensus
  ##############################################################################

  cluster_consensus_y <- cluster_color <- number_clusters_lab <- NULL
  # bring the cluster consensus data into the correct format
  previous_number_cluster <- cluster_consensus[1, 1]
  for (i in seq_len(length(colors_arr))) {
    # if the current number of predefined clusters (in the previous loops called
    # tk) is not the same as the previous, then insert 0s as space between the
    # different numbers of clusters on the x axis
    if (previous_number_cluster != cluster_consensus[i, 1]) {
      cluster_consensus_y <- c(cluster_consensus_y, 0, 0)
      cluster_color <- c(cluster_color, NA, NA)
      previous_number_cluster <- cluster_consensus[i, 1]
      number_clusters_lab <- c(number_clusters_lab, NA, NA)
    }
    cluster_consensus_y <- c(cluster_consensus_y, cluster_consensus[i, 3])
    cluster_color <- c(cluster_color, colors_arr[i])
    number_clusters_lab <- c(number_clusters_lab, cluster_consensus[i, 1])
  }
  names(cluster_consensus_y) <- number_clusters_lab
  # no need to store the parameters here, as the original mfrow and mar
  # parameters are stored and restored on exit already earlier in this function
  par(mfrow = c(3, 1), mar = c(4, 3, 2, 0))
  if ("cluster_consensus" %in% which_plots) {
    barplot(cluster_consensus_y, col = cluster_color, border = cluster_color,
            main = "cluster-consensus", ylim = c(0, 1), las = 1)
  }
}
