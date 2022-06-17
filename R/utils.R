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
  if (is.null(rownames(count_matrix))) stop("'count_matrix' should have gene name as row name")
  if (is.null(colnames(count_matrix))) stop("'count_matrix' should have cell name as column name")

  if (length(setdiff(gene_list, rownames(count_matrix)))>0) stop("There are genes in 'gene_list' that are not in 'count_matrix'")

  if (!identical(colnames(cell_cluster_conversion), c("cell_name", "class_label"))) stop("'cell_cluster_conversion' should have column name as 'cell_name' and 'class_label'")

  if (length(cluster_label)>1) stop("length of 'cluster_label' should be 1")
  if (length(setdiff(cluster_label, unique(cell_cluster_conversion$class_label)))>0) stop("'cluster_label' is not a cell type in 'cell_cluster_conversion'")

  cells_in_cluster=colnames(count_matrix)[which(as.character(cell_cluster_conversion[colnames(count_matrix),"class_label"]) %in% cluster_label)]
  total_count=sum(count_matrix[,cells_in_cluster])
  proportion=rowSums(as.matrix(count_matrix[gene_list,cells_in_cluster]))/total_count
  return(proportion)
}


#' Calculate average expression of genes in a given cell type
#'
#' @description average_count calculates average expression of genes in a given cell type
#'
#' @param count_matrix A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name
#' @param gene_list A character vector of gene names for which we want to calculate the average expression. The gene name should match the row name of \code{count_matrix}
#' @param cluster_label A character vector with one element representing the cell type name
#' @param cell_cluster_conversion A data frame with each row representing information of one cell. First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#'
#' @return A numeric vector of average expression of each gene in the cell type specified by cluster_label
#' @export
#'
#' @examples
#' data(spatial_count)
#' data(spatial_cluster)
#' cell.type = unique(spatial_cluster$class_label)[1]
#' average_expr = average_count(spatial_count, rownames(spatial_count), cell.type, spatial_cluster)
#' average_expr
average_count=function(count_matrix, gene_list, cluster_label, cell_cluster_conversion){
  if (is.null(rownames(count_matrix))) stop("'count_matrix' should have gene name as row name")
  if (is.null(colnames(count_matrix))) stop("'count_matrix' should have cell name as column name")

  if (length(setdiff(gene_list, rownames(count_matrix)))>0) stop("There are genes in 'gene_list' that are not in 'count_matrix'")

  if (!identical(colnames(cell_cluster_conversion), c("cell_name", "class_label"))) stop("'cell_cluster_conversion' should have column name as 'cell_name' and 'class_label'")

  if (length(cluster_label)>1) stop("length of 'cluster_label' should be 1")
  if (length(setdiff(cluster_label, unique(cell_cluster_conversion$class_label)))>0) stop("'cluster_label' is not a cell type in 'cell_cluster_conversion'")

  cells_in_cluster=colnames(count_matrix)[which(as.character(cell_cluster_conversion[colnames(count_matrix),"class_label"]) %in% cluster_label)]
  ave_count=rowSums(as.matrix(count_matrix[gene_list,cells_in_cluster]))/length(cells_in_cluster)
  return(ave_count)
}


#' A scatterplot with color corresponding to density
#'
#' @description Visualize two dimensional data using scatterplot with density information
#' @param x A numeric vector containing the \code{x} coordinates of points in the plot
#' @param y A numeric vector containing the \code{y} coordinates of points in the plot
#' @param xlab a title for the x axis. Default is NULL
#' @param ylab a title for the y axis. Default is NULL
#' @param main an overall title for the plot. Default is NULL
#' @param add.diagonal a logical value indicating whether to add \code{x = y} line to the plot. Default is FALSE
#' @param x_low A numeric value representing the lower limit of \code{x} coordinate. Default is NULL
#' @param x_high A numeric value representing the upper limit of \code{x} coordinate. Default is NULL
#' @param y_low A numeric value representing the lower limit of \code{y} coordinate. Default is NULL
#' @param y_high A numeric value representing the upper limit of \code{y} coordinate. Default is NULL
#'
#' @export
#'
#' @examples
#' require(ggplot2)
#' x = rnorm(1000, mean = 0, sd = 1)
#' y = rnorm(1000, mean = 0, sd = 1)
#' densityscatter(x, y, xlab = "x axis", ylab = "y axis", main = "example plot", add.diagonal = TRUE,
#'                x_low = 0, x_high = 3, y_low = 0, y_high = 3)
densityscatter=function(x, y, xlab = NULL, ylab = NULL, main = NULL, add.diagonal = F,
                        x_low = NULL, x_high = NULL, y_low = NULL, y_high = NULL){
  if (!is.numeric(x)) stop("'x' should be a numeric vector")
  if (!is.numeric(y)) stop("'y' should be a numeric vector")
  if (length(x) != length(y)) stop("'x' and 'y' should have the same length")

  df=data.frame(x=x, y=y)
  if (is.null(x_low)){
    x_low = floor(min(df$x))
  }
  if (is.null(x_high)){
    x_high = ceiling(max(df$x))
  }
  if (is.null(y_low)){
    y_low = floor(min(df$y))
  }
  if (is.null(y_high)){
    y_high = ceiling(max(df$y))
  }
  # y_high = 0
  p = ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y)) +
    ggpointdensity::geom_pointdensity() +
    viridis::scale_color_viridis() +
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(main) + ggplot2::xlim(x_low, x_high) + ggplot2::ylim(y_low, y_high) +
    ggplot2::theme(
      # axis
      axis.title.x = ggplot2::element_text(face="bold", size=16, margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = ggplot2::element_text(face="bold", size=16, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text.x = ggplot2::element_text(hjust = 0.5, vjust = 0, size=18),
      axis.text.y = ggplot2::element_text(vjust = 0.5, hjust = 0, size=18),
      axis.line = ggplot2::element_line(colour = "black"),
      #background
      panel.background = ggplot2::element_blank()
    )
  if (add.diagonal){
    p = p + ggplot2::geom_abline(intercept = 0, slope = 1)
  }
  print(p)
}

#' Data transformation
#'
#' @description data_transformation performs basic data transformations, e.g., log transformation
#' @param input_data A numeric vector of values to transform
#' @param trans_type Type of transformation:
#'
#' * \code{log}: log transformation
#' * \code{square_root}: square root transformation
#' * \code{logit}: logit transformation
#' * \code{none}: no transformation
#' @param base Base of log transformation. Default is 10
#'
#' @details For log and logit transofmration, if the value before transformation is 0, we will replace it by a small value before the transformation to avoid negative infinity. The small value is calculated as the minimum of all non-zero values from the input divided by 2
#' @return A numeric vector of transformed values
#' @export
#'
#' @examples
#' x = sample(1:100, size = 10)
#' data_transformation(x, trans_type = "log", base=10)
data_transformation=function(input_data ,trans_type, base=10){
  if (!is.numeric(input_data)) stop("'input_data' should be a numeric vector")
  if (length(setdiff(trans_type, c("log","square_root","logit","none")))>0) stop("'trans_type' should be one of 'log', 'logit', 'square_root', or 'none'")

  if (trans_type=="none"){
    return(input_data)
  }
  if (trans_type=="log"){
    if (sum(input_data==0)>0){             #if we have zeros
      #find the minimum value besides 0
      temp=min(input_data[input_data!=0], na.rm=T)    #find a pseudocount
      if (base==10){
        return(log10(input_data+temp/2))
      }
      if (base==exp(1)){
        return(log(input_data+temp/2))
      }
    }
    if (sum(input_data==0)==0){            #if there is no zero
      if (base==10){
        return(log10(input_data))            #direct log transformation without adding pseudocount
      }
      if (base==exp(1)){
        return(log(input_data))
      }
    }
  }
  if (trans_type=="square_root"){
    return(sqrt(input_data))
  }
  if (trans_type=="logit"){
    if (sum(input_data==0)>0){
      #find the minimum value besides 0
      temp=min(input_data[input_data!=0], na.rm=T)    #find a pseudocount

      return(log((input_data+temp/2)/(1-(input_data+temp/2))))
    }
    if (sum(input_data==0)==0){
      return(log(input_data/(1-input_data)))
    }
  }
}


#' Check distortion of gene expression between scRNA-seq and spatial transcriptomics
#'
#' @description distortion_test performs linear regression and deming regression on gene expression from scRNA-seq and spatial transcriptomics to check their agreement
#' @param prop_vector A numeric vector containing the expression level of one gene from two data modalities
#' @param group_label A character vector indicating which value is from which data modality (\code{sc} for scRNA-seq and \code{spatial} for spatial transcriptomics)
#'
#' @return Coefficients from linear and deming regression
#' @export
#'
#' @examples
#' x = sample(1:10, 100, replace=TRUE)
#' y = sample(100:1000, 100, replace=TRUE)
#' label = c(rep("sc", 100), rep("spatial", 100))
#' distortion_test(c(x,y), label)
distortion_test=function(prop_vector, group_label){
  if (!is.numeric(prop_vector)) stop("'prop_vector' should be a numeric vector")
  if (!is.character(group_label)) stop("'group_label' should be a character vector")
  if (length(prop_vector) != length(group_label)) stop("'prop_vector' and 'group_label' should have the same length")
  if (length(unique(group_label)) != 2) stop("'group_label' should contain only two unique labels")

  prop_sc=as.numeric(prop_vector)[which(group_label=="sc")]
  prop_spatial=as.numeric(prop_vector)[which(group_label=="spatial")]
  dat=data.frame(x=prop_sc, y=prop_spatial)

  #fit a linear regression model
  lr = stats::lm(y ~ 0 + x, data = dat)
  #lr=lm(dat$y ~ dat$x, data=dat)
  beta.regular = as.numeric(lr$coefficients[1])
  #intercept=as.numeric(lr$coefficients[1])

  #fit a Deming regression model to account for regression dilution
  fit = try(deming::deming(y ~ 0 + x, data=dat), silent = T)
  if (class(fit)=="try-error"){
    beta.deming=NA
  }else{
    beta.deming = as.numeric(fit$coefficients[2])
  }
  return(c(beta.regular, beta.deming))
}




#' Clustering of cell types
#'
#' @description cluster_dis first calculates distance between cell types and then performs hierarchical clustering.
#' Specifically, expression of cells from the same cell type is averaged to represent the expression level of a cell type.
#' Then, highly variable genes are selected based on each gene's average expression across cell types.
#' Expression information of top highly expressed genes is used to calculate distance between cell types.
#' The calculated distance is finally used for hierarchical clustering.
#'
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name
#' @param cell_cluster_conversion A data frame with each row representing information of one cell. First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#' @param dist_metric A character specifying the metric for distance calculation. Default is "correlation", which uses 1 - Pearson correlation coefficient as distance. Other options are "euclidean" for Euclidean distance and "manhattan" for Manhattan distance.
#' @param cluster_metric A character specifying the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete" (default), "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param top.var.gene A numeric value specifying the number of top highly variable genes to use for distance calculation. Default is 1000
#' @param log.transform A logical value indicating whether to perform log transformation on input count data. Default is TRUE.
#'
#' @return A list containing 3 slots:
#' \item{distance_matrix}{A \code{dist} object containing the distance between all pairs of cell types.}
#' \item{hierarchy}{Hierarchical clustering result returned by \code{hclust}.}
#' \item{dendrogram}{Clustering dendrogram.}
#' @export
#'
#' @examples
#' data(norm_sc_count)
#' data(sc_cluster)
#' cluster_distance=cluster_dis(count_table = norm_sc_count, cell_cluster_conversion = sc_cluster,
#'                              cluster_metric = "complete", dist_metric = "correlation",
#'                              log.transform = FALSE, top.var.gene = 1000)
cluster_dis=function(count_table, cell_cluster_conversion, dist_metric = "correlation", cluster_metric = "complete", top.var.gene = 1000, log.transform = T){
  if (is.null(rownames(count_table))) stop("'count_table' should have gene name as row name")
  if (is.null(colnames(count_table))) stop("'count_table' should have cell name as column name")

  if (!identical(colnames(count_table), rownames(cell_cluster_conversion))) stop("column name of 'count_table' should match the row name of 'cell_cluster_conversion'")

  if (!identical(colnames(cell_cluster_conversion), c("cell_name", "class_label"))) stop("'cell_cluster_conversion' should have column name as 'cell_name' and 'class_label'")

  if (!(dist_metric %in% c("correlation", "euclidean", "manhattan"))) stop("'dist_metric' should be one of 'correlation', 'euclidean', or 'manhattan'")

  cluster_label=base::unique(as.character(cell_cluster_conversion[colnames(count_table),"class_label"]))
  metacellcount=sapply(cluster_label, merge_column, count_table=count_table, cell_cluster_conversion=cell_cluster_conversion)         #we collapse cells in a cluster into one cell

  #remove genes with 0 average count per cluster
  num.0=apply(metacellcount, 1, function(x) return(sum(x==0)))
  metacellcount=metacellcount[which(num.0==0),]

  #log transform data
  if (log.transform){
    logmetacellcount=log10(metacellcount)      #we may have -Inf here if we have 0 in the count_table
  }else{
    logmetacellcount=metacellcount
  }

  #we use the expression of the top top.var.gene genes for distance calculation
  if (top.var.gene>=dim(logmetacellcount)[1]){
    top.var.gene=dim(logmetacellcount)[1]
  }
  var.per.gene=apply(logmetacellcount, 1, stats::var)
  top.logmetacellcount=logmetacellcount[order(var.per.gene, decreasing=T)[1:top.var.gene],]

  data2compute=t(top.logmetacellcount)           #clustering is performed on rows

  #calculate distance
  if (dist_metric=="correlation"){
    dist_matrix=stats::as.dist(1-stats::cor(top.logmetacellcount))
  }
  if (dist_metric=="euclidean" || dist_metric=="manhattan"){
    dist_matrix=stats::dist(data2compute, method=dist_metric)
  }
  #hierarchical cluster
  hc <- stats::hclust(dist_matrix, method = cluster_metric)
  # Create dendrogram
  dend <- stats::as.dendrogram (hc)
  return(list(distance_matrix=dist_matrix, hierarchy=hc, dendrogram=dend))       #dist_matrix here is a matrix with cell type names as row and column names
}


#' Merge expression of cells from the same cell type
#'
#' @description merge_column calculates expression of genes in a cell type by taking the average from cells of the same cell type
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name
#' @param cell_cluster_conversion A data frame with each row representing information of one cell. First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#' @param current_label A character specifying the current cell type name
#'
#' @return A numeric vector with average gene expression for the current cell type
#' @export
#'
merge_column=function(count_table, cell_cluster_conversion, current_label){
  cells=which(cell_cluster_conversion[colnames(count_table),"class_label"] %in% current_label)
  result=rowSums(count_table[,cells])/length(cells)
  return(result)
}


#' Construct a weighted penalty matrix based on cell type hierarchy
#'
#' @description hierarchical_penalty constructs a weighted penalty matrix based on cell type hierarchy.
#'
#' @param weight.matrix A square matrix specifying the pairwise distance between cell types
#' @param cell.type.hierarchy A data frame containing the cell type hierarchy information. Specifically, each row represents one cell type and each column contains the corresponding cell type annotation at a given granularity.
#' @param reference.resolution A character specifying the level of granularity that is used as reference. It must be one of the columns in \code{cell.type.hierarchy}
#' @param current.resolution A character specifying the current level of granularity. It is the level of cell type annotation that is used for gene panel selection. It must be one of the columns in \code{cell.type.hierarchy} and should be a lower level compared to \code{reference.resolution}
#' @param penalty A numeric value specifying the penalty for misclassifications across cell types at the \code{reference.resolution} level. Default is 1.
#'
#' @details During genetic algorithm optimization, a weighted penalty matrix can be provided to give partial credit or extra penalty to classification between certain cell types to reflect custom preference.
#' The matrix is a square matrix with each row and each column representing one cell type. This is the same with the confusion matrix from classification. For each value \eqn{p_ij} off-diagonal in the weighted penalty matrix,
#' if \eqn{p_ij > 1}, an extra penalty is given to misclassifying cells from cell type \eqn{j} to cell type \eqn{i}. If \eqn{p_ij < 1}, a partial credit is given to misclassifying cells from cell type \eqn{j} to cell type \eqn{i}.
#' \eqn{p_ij = 1} means no penalty or partial credit. The weighted penalty matrix is incorporated to the confusion matrix by element-wise multiplication to provide a weighted confusion matrix.
#'
#' Essentially, the weighted penalty matrix can be constructed arbitrarily by user’s preference. \code{hierarchical_penalty} constructs the weighted penalty matrix from cell type hierarchy.
#' Input is a square matrix specifying the pairwise distance between cell types (\code{weight.matrix}).
#' The pairwise distance matrix is then adjusted based on cell type hierarchy (\code{cell.type.hierarchy}) to construct the weighted penalty matrix.
#' Specifically, cell types can be organized in a hierarchical manner. Each level represents a different granularity
#' with higher level representing broad cell types (e.g., neuron) and lower level representing detailed cell types that are further divided from broad cell types, i.e., subpopulations (e.g., excitatory neuron and inhibitory neuron).
#' To adjust the pairwise distance matrix based on cell type hierarchy, a reference level is specified.
#' For cell types below the reference level that are from the same cell type at the reference level, the pairwise distance (between 0 and 1) between them is kept unchanged to reflect partial credit to wrong classifications among them.
#' For cell types below the reference level that are from different cell types at the reference level, the pairwise distance between them is changed to \code{penalty} to reflect extra penalty.
#' This weighted penalty matrix can then be used for hierarchical classification.
#'
#' @return A weighted penalty matrix specifying the partial credit and extra penalty for correct and incorrect classifications between pairs of cell types. This is based on cell type hierarchy information
#' @export
#'
#' @examples
#' data(norm_sc_count)
#' data(sc_cluster)
#' cluster_distance=cluster_dis(count_table = norm_sc_count, cell_cluster_conversion = sc_cluster,
#'                              cluster_metric = "complete", dist_metric = "correlation",
#'                              log.transform = FALSE, top.var.gene = 1000)
#' raw_weight_penalty = as.matrix(cluster_distance$distance_matrix)
#' raw_weight_penalty = raw_weight_penalty/max(raw_weight_penalty)
#' weight_penalty = hierarchical_penalty(weight.matrix = raw_weight_penalty, cell.type.hierarchy = cell_type_hierarchy,
#'                                       reference.resolution = "class", current.resolution = "subclass", penalty = 2)
#' diag(weight_penalty)=1
hierarchical_penalty=function(weight.matrix, cell.type.hierarchy, reference.resolution, current.resolution, penalty = 1){
  if (!is.numeric(penalty)) stop("'penalty' needs to be numeric")

  if (dim(weight.matrix)[1] != dim(weight.matrix)[2]) stop("'weight.matrix' needs to be a square matrix")

  if (!identical(rownames(weight.matrix), colnames(weight.matrix))) stop("row name of 'weight.matrix' should be identical to its column name")

  if (!(reference.resolution %in% colnames(cell.type.hierarchy) && current.resolution %in% colnames(cell.type.hierarchy))) stop("'reference.resolution' and 'current.resolution' need to be column names of cell.type.hierarchy")

  cell.type.list=rownames(weight.matrix)
  for (cell.type in cell.type.list){
    #get the reference cell type that the current cell type belongs to
    reference.cell.type = base::unique(cell.type.hierarchy[[reference.resolution]][which(cell.type.hierarchy[[current.resolution]] == cell.type)])

    #get all the cell types that belong to the same reference cell type (including the current cell type itself)
    same.reference.cell.types = unique(base::subset(cell.type.hierarchy, cell.type.hierarchy[[reference.resolution]] %in% reference.cell.type)[[current.resolution]])

    #get all the cell types that don't belong to the same reference cell type
    different.reference.cell.types = base::setdiff(cell.type.list, same.reference.cell.types)

    #change the weight between the current cell type vs. cell types that don't belong to the same reference cell type
    weight.matrix[cell.type, different.reference.cell.types] = penalty
  }
  return(weight.matrix)
}
