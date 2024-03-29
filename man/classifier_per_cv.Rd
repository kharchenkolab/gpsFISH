% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpsFISH_optimize.R
\name{classifier_per_cv}
\alias{classifier_per_cv}
\title{Perform classification for each cross validation}
\usage{
classifier_per_cv(
  current_round,
  cvlabel,
  data4cv,
  class_label_per_cell,
  metric = "Accuracy",
  method = "NaiveBayes",
  relative_prop = NULL,
  sample_new_levels = NULL,
  use_average_cluster_profiles = FALSE,
  simulation_type,
  simulation_parameter,
  simulation_model = "ZINB",
  cell_cluster_conversion,
  weight_penalty = NULL
)
}
\arguments{
\item{current_round}{A character vector specifying the name of the current cross validation. Format should be "Fold1", "Fold2", etc.}

\item{cvlabel}{A list with row indices per fold returned by create_folds.}

\item{data4cv}{A numeric matrix containing the expression per gene per cell with gene name as row name and cell name as column name.}

\item{class_label_per_cell}{A character vector specifying the cell type of each cell.}

\item{metric}{A character specifying the metric to use for evaluating the gene panel's classification performance.
Default is "Accuracy", which is the overall accuracy of classification.
The other options is "Kappa", which is the Kappa statistics.}

\item{method}{A character specifying the classification method to use. Default is naive Bayes ("NaiveBayes"). The other option is random forest ("RandomForest").}

\item{relative_prop}{A list with two elements:
\itemize{
\item "cluster.average": A matrix containing the relative expression of each gene in each cell type with gene name as row name and cell type name as column name.
The denominator for relative expression calculation needs to be all genes in the transcriptome before filtering out lowly expressed genes.
\item "cell.level": A matrix containing the relative expression of each gene in each cell with gene name as row name and cell name as column name.
The denominator for relative expression calculation needs to be all genes in the transcriptome before filtering out lowly expressed genes.
}}

\item{sample_new_levels}{A character specifying how simulation is performed for genes we have observed in the data used to train the Bayesian model.
Specifically, during the training of the Bayesian model, we have estimations of the platform effect for genes in the training data.
When we simulate spatial transcriptomics data for genes we have already seen in the training data,
we can use their estimated platform effect \eqn{\gamma_i} and \eqn{c_i} ("old_levels"), or we can randomly sample \eqn{\gamma_i} and \eqn{c_i} from their posterior distribution ("random").
For genes not in the training data, we will randomly sample \eqn{\gamma_i} and \eqn{c_i} from their posterior distribution since we don't have their estimation.}

\item{use_average_cluster_profiles}{A logical value indicating if we want to use relative expression per cell type as the relative expression input for simulating spatial transcriptomics data.
If TRUE, then value in \code{relative_prop$cluster.average} is used. And for each gene, cells from the same cell type will have the same value.
If FALSE, then we use the relative expression per cell as input, i.e., \code{relative_prop$cell.level}.
Default is FALSE.}

\item{simulation_type}{A character specifying whether simulation is performed. There are two options. "No_simulation" means no simulation. "Simulation" means we do simulation.}

\item{simulation_parameter}{A simulation model returned by simulation_training_ZINB_trim.}

\item{simulation_model}{A character specifying the type of simulation model. Default is the Bayesian model ("ZINB").}

\item{cell_cluster_conversion}{A character vector containing the cell type of each cell.}

\item{weight_penalty}{Optional. A weighted penalty matrix specifying the partial credit and extra penalty for correct and incorrect classifications between pairs of cell types.
It should be a square matrix with cell types as both row and column name. Default is NULL.}
}
\value{
A list with elements:
\item{fitness_per_cv}{Fitness of the current gene panel in the current cross validation.}
\item{confusion.matrix}{Confusion matrix of the current gene panel in the current cross validation.}
\item{norm.confusion.matrix}{Normalized confusion matrix of the current gene panel in the current cross validation.}
\item{statsbyclass}{Classification statistics of the current gene panel in the current cross validation.}
\item{pred.prob}{Predicted probability of cells based on the current gene panel in the current cross validation.}
\item{AUC.byclass}{AUC per cell type of the current gene panel in the current cross validation.}
}
\description{
Perform classification for each cross validation
}
