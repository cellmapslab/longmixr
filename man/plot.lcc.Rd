% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.lcc.R
\name{plot.lcc}
\alias{plot.lcc}
\title{Plot a longitudinal consensus clustering}
\usage{
\method{plot}{lcc}(x, color_palette = NULL, which_plots = "all", n_item_consensus = 3, ...)
}
\arguments{
\item{x}{\code{lcc} object (output from \code{\link{longitudinal_consensus_cluster}})}

\item{color_palette}{optional character vector of colors for consensus matrix}

\item{which_plots}{determine which plots should be plotted; the default is \code{"all"}.
Alternatively, a combination of the following values can be specified to plot
only some of the below mentioned plots: \code{"consensusmatrix_legend"},
\code{"consensusmatrix_x"} where \code{x} is replaced by the corresponding number
of clusters, \code{"CDF"}, \code{"delta"}, \code{"cluster_tracking"},
\code{"item_consensus"} or \code{"cluster_consensus"}. When you want to plot
all consensus matrices and the legend, you can just use \code{"consensusmatrix"}.}

\item{n_item_consensus}{determines how many item consensus plots are plotted
together in one plot before a new plot is used; the default is \code{3}.}

\item{...}{additional parameters for plotting; currently not used}
}
\value{
Plots the following plots (when selected):\tabular{ll}{
\code{consensus matrix legend} \tab the legend for the following consensus matrix plots (select with \code{"consensusmatrix_legend"}) \cr
\tab \cr
\code{consensus matrix plot} \tab for every specified number of clusters, a heatmap of the consensus matrix and the result of the final clustering is shown (select with \code{"consensusmatrix_x"} where \code{x} is replaced by the corresponding number
of clusters) \cr
\tab \cr
\code{consensus CDF} \tab a line plot of the CDFs for all different specified numbers of clusters (select with \code{"CDF"})\cr
\tab \cr
\code{Delta area} \tab elbow plot of the difference in the CDFs between the different numbers of clusters (select with \code{"delta"}) \cr
\tab \cr
\code{tracking plot} \tab cluster assignment of the subjects throughout the different cluster solutions (select with \code{"cluster_tracking"}) \cr
\tab \cr
\code{item-consensus} \tab for every item (subject), calculate the average consensus value with all items that are assigned to one consensus cluster. This is repeated for every cluster and for all different numbers of clusters (select with \code{"item_consensus"}) \cr
\tab \cr
\code{cluster-consensus} \tab every bar represents the average pair-wise item-consensus within one consensus cluster (select with \code{"cluster_consensus"})
}
}
\description{
Plot a longitudinal consensus clustering
}
