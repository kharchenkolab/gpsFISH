#' Simulation spatial transcriptomics data from scRNA-seq data
#'
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name.
#' @param cell_cluster_conversion A data frame with each row representing information of one cell.
#' First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#' @param relative_prop A list with two elements:
#'   * "cluster.average": A matrix containing the relative expression of each gene in each cell type with gene name as row name and cell type name as column name.
#'   The denominator for relative expression calculation needs to be all genes in the transcriptome before filtering out lowly expressed genes.
#'   * "cell.level": A matrix containing the relative expression of each gene in each cell with gene name as row name and cell name as column name.
#'   The denominator for relative expression calculation needs to be all genes in the transcriptome before filtering out lowly expressed genes.
#' @param sample_new_levels A character specifying how simulation is performed for genes we have observed in the data used to train the Bayesian model.
#' Specifically, during the training of the Bayesian model, we have estimations of the platform effect for genes in the training data.
#' When we simulate spatial transcriptomics data for genes we have already seen in the training data,
#' we can use their estimated platform effect \eqn{\gamma_i} and \eqn{c_i} ("old_levels"), or we can randomly sample \eqn{\gamma_i} and \eqn{c_i} from their posterior distribution ("random").
#' For genes not in the training data, we will randomly sample \eqn{\gamma_i} and \eqn{c_i} from their posterior distribution since we don't have their estimation.
#' @param use_average_cluster_profiles A logical value indicating if we want to use relative expression per cell type as the relative expression input for simulating spatial transcriptomics data.
#' If TRUE, then value in \code{relative_prop$cluster.average} is used. And for each gene, cells from the same cell type will have the same value.
#' If FALSE, then we use the relative expression per cell as input, i.e., \code{relative_prop$cell.level}.
#' Default is FALSE.
#' @param simulation_type A character specifying whether simulation is performed. There are two options. "No_simulation" means no simulation. "Simulation" means we do simulation.
#' @param simulation_parameter A simulation model returned by simulation_training_ZINB_trim.
#' @param simulation_model A character specifying the type of simulation model. Default is the Bayesian model ("ZINB").
#'
#' @return A matrix containing the simulated spatial transcriptomics data. Rows and columns are the same with \code{count_table}.
#' @export
#'
#' @examples
#' data(sc_count)
#' data(sc_cluster)
#' data(simulation_params)
#'
#' #calculate relative proportion
#' unique_cluster_label=as.character(unique(sc_cluster$class_label))
#' #cluster-wise relative proportion
#' relative_prop = list()
#' relative_prop[["cluster.average"]] = sapply(unique_cluster_label,
#'                                             gene_list=rownames(sc_count),
#'                                             relative_freq,
#'                                             count_matrix=sc_count,
#'                                             cell_cluster_conversion=sc_cluster)
#' #individual cell level relative proportion
#' relative_prop[["cell.level"]] = t(t(sc_count)/colSums(sc_count))
#'
#' #simulation spatial transcriptomics data
#' genes2simulate = sample(rownames(sc_count), 100)
#' simulate_sc_count = sc2spatial(count_table = sc_count,
#'                                cell_cluster_conversion = sc_clkuster,
#'                                relative_prop = relative_prop,
#'                                sample_new_levels = "old_levels",
#'                                use_average_cluster_profiles = FALSE,
#'                                simulation_type = "Simulation",
#'                                simulation_parameter = simulation_params,
#'                                simulation_model = "ZINB")
sc2spatial=function(count_table,
                    cell_cluster_conversion = NULL,
                    relative_prop = NULL,
                    sample_new_levels = NULL,
                    use_average_cluster_profiles = FALSE,
                    simulation_type,
                    simulation_parameter = NULL,
                    simulation_model = "ZINB"){
  if (simulation_type=="No_simulation"){
    result = count_table
    return(result)
  }
  if (simulation_type=="Simulation"){
    if (is.null(rownames(count_table))) stop("'count_table' should have gene name as row name")
    if (is.null(colnames(count_table))) stop("'count_table' should have cell name as column name")

    if (!identical(colnames(cell_cluster_conversion), c("cell_name", "class_label"))) stop("'cell_cluster_conversion' should have column name as 'cell_name' and 'class_label'")

    if (length(base::setdiff(colnames(count_table), rownames(cell_cluster_conversion)))>0) stop("There are cells in 'count_table' that are not in 'cell_cluster_conversion'")

    # if (simulation_model=="Naive_simulation_with_probe_failure" || simulation_model=="Naive_simulation_without_probe_failure"){
    #   result=simulation_naive(count_table = count_table, simulation_parameter = simulation_parameter, relative_prop = relative_prop)
    # }
    if (simulation_model=="ZINB"){
      result=simulation_ZINB(count_table = count_table, cell_cluster_conversion = cell_cluster_conversion, relative_prop = relative_prop, simulation_parameter = simulation_parameter, sample_new_levels = sample_new_levels, use_average_cluster_profiles = use_average_cluster_profiles)
    }
    return(result)
  }
}

