#' Simulation spatial transcriptomics data from scRNA-seq data
#'
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name.
#' @param cell_cluster_conversion A character vector containing the cell type of each cell.
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
#' class_label_per_cell = as.character(sc_cluster[colnames(sc_count),"class_label"])
#' simulate_sc_count = sc2spatial(count_table = sc_count[genes2simulate,],
#'                                cell_cluster_conversion = class_label_per_cell,
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

    if (dim(count_table)[2] != length(cell_cluster_conversion)) stop("Length of 'cell_cluster_conversion' should equal to the number of columns in 'count_table'")

    # if (simulation_model=="Naive_simulation_with_probe_failure" || simulation_model=="Naive_simulation_without_probe_failure"){
    #   result=simulation_naive(count_table = count_table, simulation_parameter = simulation_parameter, relative_prop = relative_prop)
    # }
    if (simulation_model=="ZINB"){
      result=simulation_ZINB(count_table = count_table, cell_cluster_conversion = cell_cluster_conversion, relative_prop = relative_prop, simulation_parameter = simulation_parameter, sample_new_levels = sample_new_levels, use_average_cluster_profiles = use_average_cluster_profiles)
    }
    return(result)
  }
}



#' Simulating spatial transcriptomics data from scRNA-seq data using Bayesian modelling
#'
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name.
#' @param cell_cluster_conversion A character vector containing the cell type of each cell.
#' @param relative_prop A list with two elements:
#'   * "cluster.average": A matrix containing the relative expression of each gene in each cell type with gene name as row name and cell type name as column name.
#'   The denominator for relative expression calculation needs to be all genes in the transcriptome before filtering out lowly expressed genes.
#'   * "cell.level": A matrix containing the relative expression of each gene in each cell with gene name as row name and cell name as column name.
#'   The denominator for relative expression calculation needs to be all genes in the transcriptome before filtering out lowly expressed genes.
#' @param simulation_parameter A simulation model returned by simulation_training_ZINB_trim.
#' @param sample_new_levels A character specifying how simulation is performed for genes we have observed in the data used to train the Bayesian model.
#' Specifically, during the training of the Bayesian model, we have estimations of the platform effect for genes in the training data.
#' When we simulate spatial transcriptomics data for genes we have already seen in the training data,
#' we can use their estimated platform effect \eqn{\gamma_i} and \eqn{c_i} ("old_levels"), or we can randomly sample \eqn{\gamma_i} and \eqn{c_i} from their posterior distribution ("random").
#' For genes not in the training data, we will randomly sample \eqn{\gamma_i} and \eqn{c_i} from their posterior distribution since we don't have their estimation.
#' @param use_average_cluster_profiles A logical value indicating if we want to use relative expression per cell type as the relative expression input for simulating spatial transcriptomics data.
#' If TRUE, then value in \code{relative_prop$cluster.average} is used. And for each gene, cells from the same cell type will have the same value.
#' If FALSE, then we use the relative expression per cell as input, i.e., \code{relative_prop$cell.level}.
#' Default is FALSE.
#'
#' @return A matrix containing the simulated spatial transcriptomics data. Rows and columns are the same with \code{count_table}.
#' @export
#'
simulation_ZINB=function(count_table, cell_cluster_conversion, relative_prop, simulation_parameter, sample_new_levels, use_average_cluster_profiles = FALSE){
  gene_list=rownames(count_table)
  cell_list=colnames(count_table)
  num_gene=length(gene_list)
  num_cell=length(cell_list)

  #get relative proportion
  if (use_average_cluster_profiles){                                         #for the same gene, all cells in the same cell type will have the same sc.prop
    sc_prop = relative_prop$cluster.average[gene_list, cell_cluster_conversion]
    colnames(sc_prop) = cell_list
  }else{
    sc_prop = relative_prop$cell.level[gene_list, cell_list]
  }

  prediction = ZINB_predict(sc_prop = sc_prop, simulation_parameter = simulation_parameter, sample_new_levels = sample_new_levels,
                            gene_list = gene_list, cell_list = cell_list, num_gene = num_gene, num_cell = num_cell)

  return(prediction$simu_count_matrix)
}




#' Simulating spatial transcriptomics data from scRNA-seq data using Bayesian modelling
#'
#' @param sc_prop A matrix containing the relative expression for each gene in each cell.
#' @param simulation_parameter A simulation model returned by simulation_training_ZINB_trim.
#' @param sample_new_levels A character specifying how simulation is performed for genes we have observed in the data used to train the Bayesian model.
#' Specifically, during the training of the Bayesian model, we have estimations of the platform effect for genes in the training data.
#' When we simulate spatial transcriptomics data for genes we have already seen in the training data,
#' we can use their estimated platform effect \eqn{\gamma_i} and \eqn{c_i} ("old_levels"), or we can randomly sample \eqn{\gamma_i} and \eqn{c_i} from their posterior distribution ("random").
#' For genes not in the training data, we will randomly sample \eqn{\gamma_i} and \eqn{c_i} from their posterior distribution since we don't have their estimation.
#' @param gene_list A character vector of gene names.
#' @param cell_list A character vector of cell names.
#' @param num_gene Number of genes.
#' @param num_cell Number of cells.
#'
#' @return A list with elements:
#'   \item{simu_count_matrix}{A matrix containing the simulated spatial transcriptomics data.}
#'   \item{alpha.global}{Randomly sampled \eqn{\alpha} from its posterior distribution.}
#'   \item{beta.global}{Randomly sampled \eqn{\beta} from its posterior distribution.}
#'   \item{zi.global}{Randomly sampled \eqn{\pi} from its posterior distribution.}
#'   \item{mu_gamma}{Randomly sampled \eqn{\mu_\gamma} from its posterior distribution.}
#'   \item{sigma_gamma}{Randomly sampled \eqn{\sigma_\gamma} from its posterior distribution.}
#'   \item{mu_c}{Randomly sampled \eqn{\mu_c} from its posterior distribution.}
#'   \item{sigma_c}{Randomly sampled \eqn{\signa_c} from its posterior distribution.}
#'   \item{gamma_i_gene}{Simulated \eqn{\gamma_i} for each gene.}
#'   \item{c_i_gene}{Simulated \eqn{c_i} for each gene.}
#'   \item{new.cell.size}{Simulated cell depth, i.e., total number of molecules per cell.}
#' @export
#'
ZINB_predict=function(sc_prop, simulation_parameter, sample_new_levels,
                      gene_list, cell_list, num_gene, num_cell){          #get estimated.random.intercept is the slowest step (profiling doesn't show this though)
  posterior = simulation_parameter$posterior
  lib.size = simulation_parameter$lib.size
  distortion.pattern = simulation_parameter$distortion_pattern

  #get posterior samples of parameters needed
  posterior.alpha = posterior$posterior.alpha
  posterior.beta = posterior$posterior.beta
  posterior.zi = posterior$posterior.zi
  posterior.mu_gamma = posterior$posterior.mu_gamma
  posterior.sigma_gamma = posterior$posterior.sigma_gamma
  posterior.mu_c = posterior$posterior.mu_c
  posterior.sigma_c = posterior$posterior.sigma_c

  #alpha, sigma, theta and pie are the same for all gene and all cells so we select one value for each
  alpha.global = base::sample(posterior.alpha, 1)
  beta.global = base::sample(posterior.beta, 1)
  zi.global = base::sample(posterior.zi, 1)
  mu_gamma = base::sample(posterior.mu_gamma, 1)
  sigma_gamma = base::sample(posterior.sigma_gamma, 1)
  mu_c = base::sample(posterior.mu_c, 1)
  sigma_c = base::sample(posterior.sigma_c, 1)

  #generate gamma_i and c_i for each data point
  if (sample_new_levels == "random"){
    gamma_i_gene = exp(stats::rnorm(n = num_gene, mean = mu_gamma, sd = sigma_gamma))
    c_i_gene = stats::rnorm(n = num_gene, mean = mu_c, sd = sigma_c)
    names(gamma_i_gene) = names(c_i_gene) = gene_list
  }

  if (sample_new_levels == "old_levels"){
    estimated_gamma_i = simulation_parameter$gamma_i_full[,"mean"]
    estimated_c_i = simulation_parameter$c_i_full[,"mean"]

    intercept.genes = base::intersect(gene_list, names(estimated_gamma_i))    #get the genes that appear in training dataset

    gamma_i_gene = exp(stats::rnorm(n = num_gene, mean = mu_gamma, sd = sigma_gamma))
    c_i_gene = stats::rnorm(n = num_gene, mean = mu_c, sd = sigma_c)
    names(gamma_i_gene) = names(c_i_gene) = gene_list

    gamma_i_gene[intercept.genes] = estimated_gamma_i[intercept.genes]      #replace randomly generated alpha_i by estimated alpha_i for those genes that show up in training dataset
    c_i_gene[intercept.genes] = estimated_c_i[intercept.genes]      #replace randomly generated alpha_i by estimated alpha_i for those genes that show up in training dataset
  }

  #calculate lambda
  lambda = boot::inv.logit(rep.col(gamma_i_gene, dim(sc_prop)[2]) * sqrt(sc_prop) + rep.col(c_i_gene, dim(sc_prop)[2]))          #for model 22.6

  #get theta for each gene and each cell
  shape = exp(beta.global + lambda)

  #random sample new cell size for each new cell from cell size distribution
  new.cell.size = base::sample(log(lib.size), num_cell, replace = T)                  #we do log transform here to save time

  #calculate mu for each gene and each cell
  mu = exp(alpha.global + log(lambda) + rep_row(new.cell.size, dim(sc_prop)[1]))

  #make predictions
  simu.count = suppressWarnings(zinb_generator(n=num_gene*num_cell,
                                               mu=mu,
                                               size=shape,     #1/dispersion, it is the same with the shape parameter in zinb model fitting
                                               pie=zi.global))      #zero inflation

  #change simu.count to matrix
  simu.count.matrix = matrix(simu.count, nrow = num_gene, ncol = num_cell)
  rownames(simu.count.matrix) = rownames(sc_prop)
  colnames(simu.count.matrix) = colnames(sc_prop)

  return(list(simu_count_matrix = simu.count.matrix,
              alpha.global = alpha.global,
              beta.global = beta.global,
              zi.global = zi.global,
              mu_gamma = mu_gamma,
              sigma_gamma = sigma_gamma,
              mu_c = mu_c,
              sigma_c = sigma_c,
              gamma_i_gene = gamma_i_gene,
              c_i_gene = c_i_gene,
              new.cell.size = new.cell.size))
}

