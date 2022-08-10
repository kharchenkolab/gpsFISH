#' @useDynLib gpsFish
NULL


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
#' }
"spatial_cluster"


#' Pre-trained Bayesian model based on the Codeluppi dataset.
#'
#' A list object returned by the simulation_training_ZINB function and then trimmed by the simulation_training_ZINB_trim function. The model is trained on the Codeluppi dataset.
#'
#' @format A list with 7 slots:
#' \describe{
#'   \item{data2fit}{A data frame containing the input data for fitting the Bayesian model. To reduce the size of simulation_params, we removed information stored in this slot since it is not used in downstream examples, e.g., simulate spatial transcriptomics data.}
#'   \item{model.summary}{A named list with elements \code{summary} and \code{c_summary}, which contain summaries of a stanfit object. To reduce the size of simulation_params, we removed information stored in this slot since it is not used in downstream examples, e.g., simulate spatial transcriptomics data.}
#'   \item{fixed.effect}{A named list with elements \code{summary} and \code{c_summary}, which contain summaries of specified parameters.}
#'   \item{c_i_full}{A matrix containing the estimated gene specific intercept \eqn{c_i} for each gene.}
#'   \item{gamma_i_full}{A matrix containing the estimated gene specific coefficient \eqn{\gamma_i} for each gene.}
#'   \item{lib.size}{A numeric vector containing the library size, i.e., total number of molecule per cell, of cells in the spatial transcriptomics data.}
#'   \item{posterior}{A list containing the extracted samples of variables in the Bayesian model from their posterior distribution.}
#' }
"simulation_params"


#' Probe count per gene.
#'
#' A data frame containing the number of probe count that can be designed to target a gene
#'
#' @format A data frame with 27308 rows and 3 columns. Each row corresponds to one gene. Of note, one gene could appear in multiple rows. For the 3 columns:
#' \describe{
#'   \item{gene_name}{Gene name.}
#'   \item{chr}{The chromosome that the gene is located.}
#'   \item{probe_count}{Number of probe count.}
#' }
"probe_count"


#' Cell type hierarchy.
#'
#' A data frame containing the cell type hierarchy information. Specifically, each column contains cell type annotation at a given granularity.
#' Columns are the left representing more broadly defined cell types and columns on the right representing more detailed subpopulations.
#' Each row represents one cell type.
#'
#' @format A data frame with 11 rows and 2 columns. Each row represents one cell type. For the 2 columns:
#' \describe{
#'   \item{class}{Cell type at the class level.}
#'   \item{subclass}{Cell type at the subclass level.}
#' }
"cell_type_hierarchy"
