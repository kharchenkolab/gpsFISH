#' Calculate relative expression of genes in a given cell type
#'
#' @description relative_freq calculates relative expression of genes in a given cell type
#'
#' @param count_matrix A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name
#' @param gene_list A character vector of gene names for which we want to calculate the relative expression. The gene name should match the row name of \code{count_matrix}
#' @param cluster_label A character vector with one element representing the cell type name
#' @param cell_cluster_conversion A data frame with each row representing information of one cell. First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#'
#' @return A numeric vector of relative expression of each gene in the cell type specified by cluster_label
#' @export
#'
#' @examples
#' data(spatial_count)
#' data(spatial_cluster)
#' cell.type = unique(spatial_cluster$class_label)[1]
#' rf = relative_freq(spatial_count, rownames(spatial_count), cell.type, spatial_cluster)
#' rf
relative_freq=function(count_matrix, gene_list, cluster_label, cell_cluster_conversion){
  cells_in_cluster=colnames(count_matrix)[which(as.character(cell_cluster_conversion[colnames(count_matrix),"class_label"]) %in% cluster_label)]
  total_count=sum(count_matrix[,cells_in_cluster])
  proportion=rowSums(as.matrix(count_matrix[gene_list,cells_in_cluster]))/total_count
  return(proportion)
}
