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
#'
#' @importFrom graphics hist lines legend
#' @importFrom grDevices rainbow
#'
#' @return a CDF plot
CDF <- function(matrix_list,
                breaks = 100) {
  # set up the plot
  plot(0, xlim = c(0, 1), ylim = c(0, 1), col = "white",
       bg = "white", xlab = "consensus index", ylab = "CDF",
       main = "consensus CDF", las = 2)

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
    # add the CDF to the plot
    lines(h$mids, h$counts, col = this_colors[i - 1], lwd = 2,
          type = "l")
  }
  legend(0.8, 0.5, legend = paste(rep("", k - 1), seq(2, k, by = 1), sep = ""),
         fill = this_colors)

  # plot the area under the CDF change
  delta_k <- area_k[1]
  for (i in 2:(length(area_k))) {
    # proportional increase relative to the previous k
    delta_k <- c(delta_k, (area_k[i] - area_k[i - 1]) / area_k[i - 1])
  }
  plot(1 + (1:length(delta_k)), y = delta_k, xlab = "k",
       ylab = "relative change in area under CDF curve",
       main = "Delta area", type = "b")
}

#' Assign colours to cluster assignments
#'
#' Internal function to assign colours to the cluster assignments;
#' taken from \code{ConsensusClusterPlus}.
#'
#' @param past_ct cluster assignments of the previous clustering (\code{k-1} numbers of clusters)
#' @param ct current cluster assignments (\code{k} numbers of clusters)
#' @param colorU vector of colour names
#' @param colorList results of a previous call to \code{setClusterColors}
#'
#' @return List with the following entries:\tabular{ll}{
#'    \code{1.} \tab vector of colours for every subject\cr
#'    \tab \cr
#'    \code{2.} \tab number; probably the number of different colours used (I haven't checked it)\cr
#'    \tab \cr
#'    \code{3.} \tab vector of unique colours in the first list entry
#' }
setClusterColors <- function(past_ct,
                             ct,
                             colorU,
                             colorList) {

  newColors = c()
  if (length(colorList) == 0) {
    newColors = colorU[ct]
    colori = 2
  }
  else {
    newColors = rep(NULL, length(ct))
    colori = colorList[[2]]
    mo = table(past_ct, ct)
    m = mo/apply(mo, 1, sum)
    for (tci in 1:ncol(m)) {
      maxC = max(m[, tci])
      pci = which(m[, tci] == maxC)
      if (sum(m[, tci] == maxC) == 1 & max(m[pci, ]) ==
          maxC & sum(m[pci, ] == maxC) == 1) {
        newColors[which(ct == tci)] = unique(colorList[[1]][which(past_ct ==
                                                                    pci)])
      }
      else {
        colori = colori + 1
        newColors[which(ct == tci)] = colorU[colori]
      }
    }
  }
  return(list(newColors, colori, unique(newColors)))
}

#' Generate the tracking plot
#'
#' Internal function to generate the tracking plot; taken from \code{ConsensusClusterPlus}.
#'
#' @param m list of assigned colours to every observation with one entry per
#' specified number of cluster
#'
#' @importFrom graphics rect segments text
#'
#' @return tracking plot
clusterTrackingPlot <- function(m) {
  plot(NULL, xlim = c(-0.1, 1), ylim = c(0, 1), axes = FALSE,
       xlab = "samples", ylab = "k", main = "tracking plot")
  for (i in 1:nrow(m)) {
    rect(xleft = seq(0, 1 - 1/ncol(m), by = 1/ncol(m)), ybottom = rep(1 -
                                                                        i/nrow(m), ncol(m)), xright = seq(1/ncol(m), 1, by = 1/ncol(m)),
         ytop = rep(1 - (i - 1)/nrow(m), ncol(m)), col = m[i,
         ], border = NA)
  }
  xl = seq(0, 1 - 1/ncol(m), by = 1/ncol(m))
  segments(xl, rep(-0.1, ncol(m)), xl, rep(0, ncol(m)), col = "black")
  ypos = seq(1, 0, by = -1/nrow(m)) - 1/(2 * nrow(m))
  text(x = -0.1, y = ypos[-length(ypos)], labels = seq(2, nrow(m) +
                                                         1, by = 1))
}

#' Generate colour palette
#'
#' Internal function to generate colour palette; taken from \code{ConsensusClusterPlus}.
#'
#' @param n number of colours
#'
#' @importFrom grDevices rgb
#'
#' @return vector of colour hex values where the first is red and the rest different
#' blue colours
myPal <- function(n = 10) {
  seq <- rev(seq(0, 255, by = 255 / (n)))
  palRGB <- cbind(seq, seq, 255)
  rgb(palRGB, maxColorValue = 255)
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
#' Internal function to generate the item consensus plot; taken from \code{ConsensusClusterPlus}.
#'
#' @param d data matrix
#' @param myc colours for the clusters
#' @param cc order of the items/subjects as in the solution for 2 clusters
#' @param title title of the plot
#'
#' @return item/subject consensus plot
rankedBarPlot <- function (d,
                           myc,
                           cc,
                           title) {
  colors = rbind()
  byRank = cbind()
  spaceh = 0.1
  for (i in 1:ncol(d)) {
    byRank = cbind(byRank, sort(d[, i], na.last = F))
    colors = rbind(colors, order(d[, i], na.last = F))
  }
  maxH = max(c(1.5, apply(byRank, 2, sum)), na.rm = T)
  barp = barplot(apply(byRank, 2, sum), col = myc[colors[, 1]],
                 space = spaceh, ylim = c(0, maxH),
                 main = paste("item-consensus", title), border = NA, las = 1)
  for (i in 2:nrow(byRank)) {
    barplot(apply(matrix(byRank[i:nrow(byRank), ], ncol = ncol(byRank)),
                  2, sum), space = spaceh, col = myc[colors[, i]],
            ylim = c(0, maxH), add = T, border = NA, las = 1)
  }
  xr = seq(spaceh, ncol(d) + ncol(d) * spaceh, (ncol(d) + ncol(d) *
                                                  spaceh)/ncol(d))
  text("*", x = xr + 0.5, y = maxH, col = myc[cc], cex = 1.4)
}
