#' Weekly DJIA negative-correlation network statistics
#'
#' Prepared inputs for the DJIA application in Chapter 4 of the dissertation.
#' The series contains 158 observations from negative-correlation networks on
#' 29 stocks, labelled weekly from January 1, 2007 through January 4, 2010.
#'
#' @format A named list with two components:
#' \describe{
#'   \item{net_stats}{A numeric matrix with 158 rows and three columns:
#'     \code{Edges}, the number of undirected edges; \code{Triangles}, the
#'     number of triangles; and \code{Homophily}, the number of edges joining
#'     stocks with the same risk-orientation label. Rows are chronological.
#'     These are original counts, before centering or standardization.}
#'   \item{date_list}{A Date vector of 158 weekly labels aligned with the
#'     rows of \code{net_stats}, from 2007-01-01 through 2010-01-04.}
#' }
#' @details Each network uses correlations computed from four consecutive
#'   rows of weekly log returns. An undirected edge is present when the
#'   correlation is strictly negative; diagonal entries are zero.
#'
#'   Risk-orientation labels are assigned once using each stock's cumulative
#'   degree over all 158 networks. Stocks whose cumulative degree is strictly
#'   above the median receive label H; all other stocks receive label M.
#'   Homophily counts edges within either label, not the proportion of such
#'   edges. The data contain statistics and dates, not adjacency matrices.
#'
#'   The preparation preserves the original script's native-row convention:
#'   rows \code{i:(i + 3)} of the ecp returns are labelled by \code{dates[i]}.
#'   The source dates run in reverse chronological order, so each window
#'   extends toward earlier dates. The resulting network series is ordered
#'   chronologically. Labels are preserved without shifting the windows.
#'
#'   The documentation filename is \code{R/DJIA.R}; the dataset object remains
#'   \code{DJIA_networks} to preserve the loading call in Chapter 4 and avoid
#'   confusion with the stock-return dataset \code{ecp::DJIA}.
#' @source Yuhui Wang's supplied \code{DJIA_net_stats_matrix.rds} and DJIA
#'   network-preparation script. The underlying stock-return data are
#'   \code{DJIA} in the ecp package:
#'   \url{https://CRAN.R-project.org/package=ecp}.
#'   Network construction follows the DJIA workflow of Kei et al. (2025):
#'   \url{https://github.com/allenkei/CPDstergm_demo/blob/main/experiment_real.R}.
#'   The repository's \code{data-raw/DJIA.R} packages the supplied counts and
#'   can reconstruct them from the returns for comparison.
#' @examples
#' data("DJIA_networks", package = "bcmlr")
#' X <- DJIA_networks$net_stats
#' date_list <- DJIA_networks$date_list
#' dim(X)
#' colnames(X)
#' range(date_list)
#' @docType data
#' @keywords datasets
"DJIA_networks"
