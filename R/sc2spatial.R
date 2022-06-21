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



#' Simulating spatial transcriptomics data from scRNA-seq data using Bayesian modelling
#'
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name.
#' @param cell_cluster_conversion A data frame with each row representing information of one cell.
#' First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
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
  ct_list=base::unique(cell_cluster_conversion)
  num_gene=length(gene_list)
  num_cell=length(cell_list)
  num_ct=length(ct_list)

  #get relative proportion
  if (use_average_cluster_profiles){                                         #for the same gene, all cells in the same cell type will have the same sc.prop
    sc_prop = relative_prop$cluster.average[gene_list, cell_cluster_conversion]
    colnames(sc_prop) = cell_list
  }else{
    sc_prop = relative_prop$cell.level[gene_list, cell_list]
  }

  prediction = ZINB_predict(sc_prop = sc_prop, simulation_parameter = simulation_parameter, sample_new_levels = sample_new_levels,
                            gene_list = gene_list, cell_list = cell_list, ct_list = ct_list, num_gene = num_gene, num_cell = num_cell, num_ct = num_ct, class_label = cell_cluster_conversion)

  return(prediction$simu_count_matrix)
}




ZINB_predict=function(sc_prop, simulation_parameter, sample_new_levels,
                      gene_list, cell_list, ct_list, num_gene, num_cell, num_ct, class_label){          #get estimated.random.intercept is the slowest step (profiling doesn't show this though)
  # fit.model = simulation_parameter$distortion_model
  posterior = simulation_parameter$posterior
  lib.size = simulation_parameter$lib.size
  distortion.pattern = simulation_parameter$distortion_pattern

  #get posterior samples of parameters needed
  # posterior.alpha = rstan::extract(fit.model, pars="alpha")$alpha
  # posterior.beta = rstan::extract(fit.model, pars="beta")$beta
  # posterior.zi = rstan::extract(fit.model, pars="zi")$zi
  # posterior.mu_gamma = rstan::extract(fit.model, pars="mu_gamma")$mu_gamma
  # posterior.sigma_gamma = rstan::extract(fit.model, pars="sigma_gamma")$sigma_gamma
  # posterior.xi_c = rstan::extract(fit.model, pars="xi_c")$xi_c
  # posterior.omega_c = rstan::extract(fit.model, pars="omega_c")$omega_c
  # posterior.kappa_c = rstan::extract(fit.model, pars="kappa_c")$kappa_c
  posterior.alpha = posterior$posterior.alpha
  posterior.beta = posterior$posterior.beta
  posterior.zi = posterior$posterior.zi
  posterior.mu_gamma = posterior$posterior.mu_gamma
  posterior.sigma_gamma = posterior$posterior.sigma_gamma
  # posterior.xi_c = posterior$posterior.xi_c
  # posterior.omega_c = posterior$posterior.omega_c
  # posterior.kappa_c = posterior$posterior.kappa_c
  posterior.mu_c = posterior$posterior.mu_c
  posterior.sigma_c = posterior$posterior.sigma_c

  #alpha, sigma, theta and pie are the same for all gene and all cells so we select one value for each
  alpha.global = sample(posterior.alpha, 1)
  beta.global = sample(posterior.beta, 1)
  zi.global = sample(posterior.zi, 1)
  mu_gamma = sample(posterior.mu_gamma, 1)
  sigma_gamma = sample(posterior.sigma_gamma, 1)
  # xi_c = sample(posterior.xi_c, 1)
  # omega_c = sample(posterior.omega_c, 1)
  # kappa_c = sample(posterior.kappa_c, 1)
  mu_c = sample(posterior.mu_c, 1)
  sigma_c = sample(posterior.sigma_c, 1)

  #generate gamma_i and c_i for each data point
  if (sample_new_levels == "random"){
    gamma_i_gene = exp(rnorm(n = num_gene, mean = mu_gamma, sd = sigma_gamma))
    # c_i_gene = rsn(n = num_gene, xi = xi_c, omega = omega_c, alpha = kappa_c)
    c_i_gene = rnorm(n = num_gene, mean = mu_c, sd = sigma_c)
    names(gamma_i_gene)=names(c_i_gene)=gene_list
  }

  if (sample_new_levels == "old_levels"){
    estimated_gamma_i = simulation_parameter$gamma_i_full[,"mean"]
    estimated_c_i = simulation_parameter$c_i_full[,"mean"]

    intercept.genes = intersect(gene_list, names(estimated_gamma_i))    #get the genes that appear in training dataset

    gamma_i_gene = exp(rnorm(n = num_gene, mean = mu_gamma, sd = sigma_gamma))
    # c_i_gene = rsn(n = num_gene, xi = xi_c, omega = omega_c, alpha = kappa_c)
    c_i_gene = rnorm(n = num_gene, mean = mu_c, sd = sigma_c)
    names(gamma_i_gene)=names(c_i_gene)=gene_list

    gamma_i_gene[intercept.genes] = estimated_gamma_i[intercept.genes]      #replace randomly generated alpha_i by estimated alpha_i for those genes that show up in training dataset
    c_i_gene[intercept.genes] = estimated_c_i[intercept.genes]      #replace randomly generated alpha_i by estimated alpha_i for those genes that show up in training dataset
  }

  #calculate lambda
  lambda = inv.logit(rep.col(gamma_i_gene, dim(sc_prop)[2]) * sqrt(sc_prop) + rep.col(c_i_gene, dim(sc_prop)[2]))          #for model 22.6

  #get theta for each gene and each cell
  shape = exp(beta.global + lambda)

  #random sample new cell size for each new cell from cell size distribution
  new.cell.size = sample(log(lib.size), num_cell, replace = T)                  #we do log transform here to save time

  #calculate miu for each gene and each cell
  miu = exp(alpha.global + log(lambda) + rep.row(new.cell.size, dim(sc_prop)[1]))

  #make predictions
  simu.count = suppressWarnings(zinb_generator(n=num_gene*num_cell,
                                               miu=miu,
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
              # xi_c = xi_c,
              # omega_c = omega_c,
              # kappa_c = kappa_c,
              mu_c = mu_c,
              sigma_c = sigma_c,
              gamma_i_gene = gamma_i_gene,
              c_i_gene = c_i_gene,
              alpha.global = alpha.global,
              new.cell.size = new.cell.size))
}

