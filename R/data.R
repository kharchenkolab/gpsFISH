#' scRNA-seq count data.
#'
#' A dataset containing the expression of 19,972 genes in 550 cells.
#'
#' @format A count matrix with 19,972 rows and 550 columns:
"sc_count"

#' scRNA-seq count data cell annotation.
#'
#' A dataset containing the cell type label of the 550 cells in the scRNA-seq data.
#'
#' @format A data frame with 550 rows and 2 columns:
#' \describe{
#'   \item{cell_name}{name of each cell}
#'   \item{class_label}{annotated cell type of each cell}
#'   ...
#' }
"sc_cluster"

#' Spatial transcriptomics count data.
#'
#' A dataset containing the expression of 33 genes in 537 cells.
#'
#' @format A count matrix with 33 rows and 537 columns:
"spatial_count"

#' Spatial transcriptomics count data cell annotation.
#'
#' A dataset containing the cell type label of the 537 cells in the spatial transcriptomics data.
#'
#' @format A data frame with 537 rows and 2 columns:
#' \describe{
#'   \item{cell_name}{name of each cell}
#'   \item{class_label}{annotated cell type of each cell}
#'   ...
#' }
"spatial_cluster"
