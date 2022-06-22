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

  if (length(base::setdiff(gene_list, rownames(count_matrix)))>0) stop("There are genes in 'gene_list' that are not in 'count_matrix'")

  if (!identical(colnames(cell_cluster_conversion), c("cell_name", "class_label"))) stop("'cell_cluster_conversion' should have column name as 'cell_name' and 'class_label'")

  if (length(cluster_label)>1) stop("length of 'cluster_label' should be 1")
  if (length(base::setdiff(cluster_label, unique(cell_cluster_conversion$class_label)))>0) stop("'cluster_label' is not a cell type in 'cell_cluster_conversion'")

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

  if (length(base::setdiff(gene_list, rownames(count_matrix)))>0) stop("There are genes in 'gene_list' that are not in 'count_matrix'")

  if (!identical(colnames(cell_cluster_conversion), c("cell_name", "class_label"))) stop("'cell_cluster_conversion' should have column name as 'cell_name' and 'class_label'")

  if (length(cluster_label)>1) stop("length of 'cluster_label' should be 1")
  if (length(base::setdiff(cluster_label, unique(cell_cluster_conversion$class_label)))>0) stop("'cluster_label' is not a cell type in 'cell_cluster_conversion'")

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
  if (length(base::setdiff(trans_type, c("log","square_root","logit","none")))>0) stop("'trans_type' should be one of 'log', 'logit', 'square_root', or 'none'")

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
#' weight_penalty = hierarchical_penalty(weight.matrix = raw_weight_penalty,
#'                                       cell.type.hierarchy = cell_type_hierarchy,
#'                                       reference.resolution = "class",
#'                                       current.resolution = "subclass",
#'                                       penalty = 2)
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


#' Adjust variance of scRNA-seq expression data using Pagoda2
#'
#' @description preprocess_normalize takes scRNA-seq expression data as input and adjust the variance to normalize the extent to which genes with different expression magnitudes will contribute to the downstream anlaysis
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name
#' @param n.core Number of cores to use. Default is 1.
#'
#' @return A list with:
#' \item{sparse.matrix}{the adjusted count matrix}
#' \item{pagoda.object}{the pagoda2 object}
#' @export
#'
#' @examples
#' data(sc_count)
#' preprocess_normalize(sc_count, n.core = 2)
preprocess_normalize=function(count_table, n.core = 1){
  sm <- Matrix::Matrix(as.matrix(count_table), sparse = TRUE)
  r <- pagoda2::Pagoda2$new(sm, log.scale = TRUE, n.cores = n.core)
  r$adjustVariance()
  return(list(sparse.matrix=sm, pagoda.object=r))
}


#' Differential gene expression between clusters
#'
#' @description diff_gene_cluster performs differential gene expression between clusters using Pagoda2
#' @param pagoda_object A pagoda2 object
#' @param cell_cluster_conversion A data frame with each row representing information of one cell. First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#' @param n.core Number of cores to use. Default is 1.
#' @param z.threshold A numeric value specifying the minimal absolute Z score (adjusted) to report. Default is 3.
#'
#' @return A list with:
#' \item{diff_result}{differential expression result}
#' \item{pagoda.obj}{the pagoda2 object}
#' @export
#'
#' @examples
#' data(sc_count)
#' data(sc_cluster)
#'
#' #adjust variance
#' sc_count.adjust_variance = preprocess_normalize(sc_count, n.core = 2)
#'
#' #run differential expression
#' diff_expr = diff_gene_cluster(pagoda_object = sc_count.adjust_variance$pagoda.object,
#'                               cell_cluster_conversion = sc_cluster, n.core = 1)
#'
#' #Get differentially expressed genes for a specific cell type, e.g., Pvm:
#' head(diff_expr$diff_result$Pvm)
diff_gene_cluster=function(pagoda_object, cell_cluster_conversion, n.core = 1, z.threshold = 3){
  if (class(pagoda_object)[1] != "Pagoda2") stop("'pagoda_object' needs to be a pagoda2 object")

  r <- pagoda_object

  group=as.character(cell_cluster_conversion[rownames(r$counts),"class_label"])
  group=as.factor((group))
  names(group)=rownames(r$counts)

  r$clusters=group

  result = r$getDifferentialGenes(groups=r$clusters, z.threshold=z.threshold, append.auc=T)
  return(list(pagoda.obj=r, diff_result=result))
}


#' Initialize population for genetic algorithm based on differential gene expression
#'
#' @description initialize_population initializes a population for genetic algorithm optimization using information from differential expression analysis.
#' The idea is that for each gene panel, we include significantlly differentially expressed genes of each cell type (if there is any).
#' Detailed information of how each gene panel is initialized can be found by \code{?initialize_solution}.
#' @param pop.size Size of the population, i.e., the number of gene panels to initialize
#' @param panel.size Size of a gene panel, i.e., the number of genes in a gene panel
#' @param diff_expr_result A list containing the differential expression result. Each slot in the list contains a data frame corresponding to differentially expressed genes in one cell type.
#' Within each data frame, each row is one gene. Columns contain differential expression statistics of a gene, e.g., "AUC", "z-score", "precision", etc.
#' @param diff_metric A character specifying the metric used to select significant differentially expressed genes. Needs to be one of the column names of the data frame in \code{diff_expr_result}, e.g., "AUC".
#' @param diff_metric_cutoff A numeric value representing the significance cutoff for differential expression.
#' @param gene.list A character vector of all candidate genes used to initialize the population
#' @param gene2include A character vector of genes that must be included in each panel of the population. Default is NULL
#'
#' @return A matrix with each row representing one gene panel and each column representing one gene in a gene panel.
#' The genes are encoded by their location in \code{gene.list}
#' @export
#'
#' @examples
#' data(sc_count)
#' data(sc_cluster)
#'
#' gene2include.symbol = sample(rownames(sc_count), 20)
#'
#' #adjust variance
#' sc_count.adjust_variance = preprocess_normalize(sc_count, n.core = 2)
#'
#' #run differential expression
#' diff_expr = diff_gene_cluster(pagoda_object = sc_count.adjust_variance$pagoda.object,
#'                               cell_cluster_conversion = sc_cluster, n.core = 1)
#'
#' initpop.DE=initialize_population(pop.size = 100,
#'                                  panel.size = 200,
#'                                  diff_expr_result = diff_expr$diff_result,
#'                                  diff_metric = "AUC",
#'                                  diff_metric_cutoff = 0.7,
#'                                  gene.list = rownames(sc_count),
#'                                  gene2include = gene2include.symbol)
#'
#' head(initpop.DE)
initialize_population = function(pop.size, panel.size, diff_expr_result, diff_metric, diff_metric_cutoff, gene.list, gene2include=NULL){
  if (!(diff_metric %in% colnames(diff_expr_result[[1]]))) stop(paste("'diff_metric' must be one of:", paste(colnames(diff_expr_result[[1]]), collapse = ", ")))

  if (length(base::setdiff(gene2include, gene.list))>0) stop("genes in 'gene2include' must also be in 'gene.list'")

  clusters = names(diff_expr_result)

  if (is.null(gene2include)){          #if we don't have genes that we must include
    initpop = sapply(seq(1:pop.size), initialize_solution,
                     panel.size = panel.size, diff_expr_result = diff_expr_result, cluster.list = clusters,
                     diff_metric = diff_metric, diff_metric_cutoff = diff_metric_cutoff,
                     gene.list = gene.list, gene2include = gene2include)
  }else{                               #if there are genes that we must include
    new.panel.size=panel.size-length(gene2include)
    initpop = sapply(seq(1:pop.size), initialize_solution,
                     panel.size = new.panel.size, diff_expr_result = diff_expr_result, cluster.list = clusters,
                     diff_metric = diff_metric, diff_metric_cutoff = diff_metric_cutoff,
                     gene.list = gene.list, gene2include = gene2include)
    initpop=base::rbind(initpop, replicate(pop.size, which(gene.list %in% gene2include)))
  }
  return(t(initpop))
}


#' Initialize a gene panel within a population for genetic algorithm based on differential gene expression
#'
#' @description initialize_solution initializes one gene panel within a population using information from differential expression analysis.
#' Specifically, given the panel size and number of cell types, we first calculate on average how many significantly differentially expressed genes we can include from each cell type.
#' Then we select this number of genes from the differentially expressed genes of each cell type (or all the differentially expressed genes if the number of them is smaller than this number).
#' We select the most significant genes first followed by less significant ones.
#' After doing this for all the cell types, for the rest of the panel, we randomly select genes from the rest of the differentially expressed genes.
#' @param sol.num A numeric value, can be any number.
#' @param panel.size Size of a gene panel, i.e., the number of genes in a gene panel
#' @param diff_expr_result A list containing the differential expression result. Each slot in the list contains a data frame corresponding to differentially expressed genes in one cell type.
#' Within each data frame, each row is one gene. Columns contain differential expression statistics of a gene, e.g., "AUC", "z-score", "precision", etc.
#' @param diff_metric A character specifying the metric used to select significant differentially expressed genes. Needs to be one of the column names of the data frame in \code{diff_expr_result}, e.g., "AUC".
#' @param diff_metric_cutoff A numeric value representing the significance cutoff for differential expression.
#' @param cluster.list A character vector of cell type names. It must be the same or a subset of \code{names(diff_expr_result)}
#' @param gene.list A character vector of all candidate genes used to initialize the population
#' @param gene2include A character vector of genes that must be included in each panel of the population. Default is NULL
#'
#' @return A numeric vector representing an initialized solution with \code{panel.size} genes.
#' The genes are encoded by their location in \code{gene.list}
#' @export
#'
initialize_solution = function(sol.num, panel.size, diff_expr_result, diff_metric, diff_metric_cutoff, cluster.list, gene.list, gene2include=NULL){
  if (!(diff_metric %in% colnames(diff_expr_result[[1]]))) stop(paste("'diff_metric' must be one of:", paste(colnames(diff_expr_result[[1]]), collapse = ", ")))

  if (length(base::setdiff(gene2include, gene.list))>0) stop("genes in 'gene2include' must also be in 'gene.list'")

  if (length(base::setdiff(cluster.list, names(diff_expr_result)))>0) stop("cell types in 'cluster.list' must also be in names(diff_expr_result)")

  candidate=c()
  gene2include.per.ct = floor(panel.size/length(cluster.list))    #calculate average number of DE genes we can include given the number of cell types and panel size
  all.DEGs=c()

  #for each cell type, add significant DE genes. We add the most significant ones first
  for(ct in cluster.list){
    DEG.table = subset(diff_expr_result[[ct]], diff_expr_result[[ct]][[diff_metric]]>diff_metric_cutoff)
    #sort DEG table by AUC
    DEG.table = DEG.table[base::order(DEG.table$AUC, decreasing = TRUE),]
    DEGs = DEG.table$Gene

    DEG.pool = base::setdiff(DEGs, gene2include)         #we select from DEGs that are not in gene2include
    all.DEGs = c(all.DEGs, DEG.pool)

    if (length(DEG.pool)>=gene2include.per.ct){
      candidate = c(candidate, sample(DEG.pool, gene2include.per.ct, prob = DEG.table[DEG.pool, diff_metric]))
    }else{
      candidate = c(candidate, DEG.pool)
    }
  }

  all.DEGs=base::unique(all.DEGs)
  candidate=base::unique(candidate)

  #for the rest of the panel, randomly select from the rest DEGs
  candidate = c(candidate, sample(base::setdiff(all.DEGs, candidate), panel.size - length(candidate)))
  return(which(gene.list %in% candidate))
}


#' Initialize population randomly from all candidate genes
#'
#' @param pop.size Size of the population, i.e., the number of gene panels to initialize
#' @param panel.size Size of a gene panel, i.e., the number of genes in a gene panel
#' @param gene.list A character vector of all candidate genes used to initialize the population
#' @param gene2include.id A numeric vector specifying the location of genes in \code{gene.list} that must be included in each panel of the population. Default is NULL
#'
#' @return A matrix with each row representing one gene panel and each column representing one gene in a gene panel.
#' The genes are encoded by their location in \code{gene.list}
#' @export
#'
#' @examples
#' data(sc_count)
#'
#' gene2include.symbol = sample(rownames(sc_count), 20)
#' gene2include.id=which(rownames(sc_count) %in% gene2include.symbol)
#'
#' initpop.random=initialize_population_random(pop.size = 100,
#'                                             panel.size = 200,
#'                                             gene.list = rownames(sc_count),
#'                                             gene2include.id = gene2include.id)
#' head(initpop.random)
initialize_population_random=function(pop.size, panel.size, gene.list, gene2include.id){
  if (length(base::setdiff(gene2include.id, 1:length(gene.list)))>0) stop("genes in 'gene2include.id' must also be in 'gene.list'")

  indices = 1:length(gene.list)
  if (!is.null(gene2include.id)){                  #if we have gene we must include
    pop1 <- t(replicate(pop.size, gene2include.id))
    if (length(gene2include.id)<panel.size){
      pop2 <- t(replicate(pop.size, sample(base::setdiff(indices, gene2include.id), (panel.size-length(gene2include.id)))))
      initpop <- cbind(pop1, pop2)
    }
    if (length(gene2include.id)==panel.size){
      initpop <- pop1
    }
  }else{
    initpop <- t(replicate(pop.size, sample(indices, panel.size)))
  }
  return(initpop)
}


#' Calculate diversity of a population
#'
#' @param x A numeric matrix with each row representing one gene panel and each column represent one gene in a gene panel
#'
#' @return A numeric value representing the diversity
#' @export
#'
#' @examples
#' pop = matrix(sample(1:1000000, 10000), 100, 100)
#' popDiv(pop)
popDiv <- function(x) {
  N <- nrow(x)
  ndiff <- 0
  for (i in 1:(N-1)) {                                   #for each solution
    more <- sapply((i+1):N, function(j) {                #calculate the similarity between this solution and all other solutions
      length(unique(c(x[i,], x[j,]))) - ncol(x)          #calculate similarity. If two solutions are identical, length(unique(c(x[i,], x[j,]))) = ncol(x)
    })
    ndiff <- ndiff+sum(more)
  }
  ndiff/(N*(N-1)/2)
}


#' Subsample scRNA-seq data
#'
#' @description subsample_sc takes scRNA-seq data as input and subsample it by selecting a subset of cells.
#'
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name.
#' @param cell_cluster_conversion A data frame with each row representing information of one cell.
#' First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#' Only needed if \code{sampling_type} is "Subsampling_by_cluster".
#' @param rate A value between 0 and 1 specifying the proportion of cells we want to keep for each cell type during subsampling. 0.8 means we keep 80% of cells for each cell type. Default is 1.
#' @param cluster_size_max Maximum number of cells to keep for each cell type during subsampling. Default is 1000. Only needed if \code{sampling_type} is "Subsampling_by_cluster".
#' @param cluster_size_min Minimum number of cells to keep for each cell type during subsampling. Default is 1. Only needed if \code{sampling_type} is "Subsampling_by_cluster".
#' @param sampling_type A character specifying how the subsample is performed.
#' "Subsampling_by_cluster" means subsampling cells from each cell type separately. "Subsampling" means subsampling all cells together by mixing cells from all cell types.
#' Default is "Subsampling_by_cluster".
#' @param nCV Number of cross validation to perform on the subsampled data. This is to ensure that we have enough cells for cross validation. Default is 1. Only needed if \code{sampling_type} is "Subsampling_by_cluster".
#'
#' @return A subsampled matrix
#' @export
#'
#' @examples
#' data(sc_count)
#' data(sc_cluster)
#'
#' #number of cells per cell type
#' table(sc_cluster[colnames(sc_count), "class_label"])
#'
#' #subsample by cell type
#' subsample_sc_count = subsample_sc(count_table = sc_count,
#'                                   cell_cluster_conversion = sc_cluster,
#'                                   rate = 1,
#'                                   cluster_size_max = 30,
#'                                   cluster_size_min = 5,
#'                                   sampling_type = "Subsampling_by_cluster",
#'                                   nCV = 5)
#'
#' #number of cells per cell type after subsample
#' table(sc_cluster[colnames(subsample_sc_count), "class_label"])
#'
#' #subsample cells from all cell types together
#' dim(sc_count)
#' subsample_sc_count = subsample_sc(count_table = sc_count,
#'                                   rate = 0.5,
#'                                   sampling_type = "Subsampling")
#' dim(subsample_sc_count)
subsample_sc=function(count_table, cell_cluster_conversion = NULL,
                      rate = 1, cluster_size_max = 1000, cluster_size_min = 1,
                      sampling_type = "Subsampling_by_cluster",
                      nCV = 1){
  if (sampling_type=="Subsampling"){
    result = subsampling_all(count_table = count_table, rate = rate)
    return(result)
  }
  if (sampling_type=="Subsampling_by_cluster"){
    result = subsampling_by_cluster(count_table = count_table, cell_cluster_conversion = cell_cluster_conversion, rate = rate, cluster_size_max = cluster_size_max, cluster_size_min = cluster_size_min, nCV = nCV)
    return(result)
  }
}


#' Subsample scRNA-seq data by mixing cells from all cell types
#'
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name.
#' @param rate A value between 0 and 1 specifying the proportion of cells we want to keep for each cell type during subsampling. 0.8 means we keep 80% of cells for each cell type. Default is 1.
#'
#' @return A subsampled matrix
#' @export
#'
#' @examples
#' data(sc_count)
#' dim(sc_count)
#' subsample_sc_count = subsampling_all(sc_count, 0.5)
#' dim(subsample_sc_count)
subsampling_all=function(count_table, rate = 1){        #subsampling from all cells
  if (rate > 1 || rate <0) stop("'rate' must be between 0 and 1")

  num.cells = base::round(dim(count_table)[2]*rate)
  selected.cells = base::sample(seq(1:dim(count_table)[2]), num.cells, replace = FALSE)
  simulated_count_table = count_table[,selected.cells]
  return(simulated_count_table)
}


#' Subsample scRNA-seq data by cell type
#'
#' @param count_table A matrix containing the expression level of each gene in each cell with gene name as row name and cell name as column name.
#' @param cell_cluster_conversion A data frame with each row representing information of one cell.
#' First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#' @param rate A value between 0 and 1 specifying the proportion of cells we want to keep for each cell type during subsampling. 0.8 means we keep 80% of cells for each cell type. Default is 1.
#' @param cluster_size_max Maximum number of cells to keep for each cell type during subsampling. Default is 1000.
#' @param cluster_size_min Minimum number of cells to keep for each cell type during subsampling. Default is 1.
#' @param nCV Number of cross validation to perform on the subsampled data. This is to ensure that we have enough cells for cross validation. Default is 1.
#'
#' @return A subsampled matrix
#' @export
#'
subsampling_by_cluster=function(count_table, cell_cluster_conversion, rate = 1, cluster_size_max = 1000, cluster_size_min = 1, nCV = 1){      #subsampling by cluster (this perserves the same number of cell clusters)
  if (is.null(rownames(count_table))) stop("'count_table' should have gene name as row name")
  if (is.null(colnames(count_table))) stop("'count_table' should have cell name as column name")

  if (!identical(colnames(cell_cluster_conversion), c("cell_name", "class_label"))) stop("'cell_cluster_conversion' should have column name as 'cell_name' and 'class_label'")

  if (length(base::setdiff(colnames(count_table), rownames(cell_cluster_conversion)))>0) stop("There are cells in 'count_table' that are not in 'cell_cluster_conversion'")

  if (rate > 1 || rate <0) stop("'rate' should between 0 and 1")

  unique_cluster_label = unique(cell_cluster_conversion$class_label)              #This works when the cells in cell_cluster_conversion are the same set of cells in count_table
  selected.cells = unlist(lapply(unique_cluster_label, subsample, cell_cluster_conversion = cell_cluster_conversion, rate = rate, cluster_size_max = cluster_size_max, cluster_size_min = cluster_size_min, nCV = nCV))
  simulated_count_table = count_table[, selected.cells]
  return(simulated_count_table)
}


#' Subsample cells from one cell type
#'
#' @param cell_cluster_conversion A data frame with each row representing information of one cell.
#' First column contains the cell name. Second column contains the corresponding cell type name. Row name of the data frame should be the cell name.
#' @param class_label A character specifying the cell type
#' @param rate A value between 0 and 1 specifying the proportion of cells we want to keep for each cell type during subsampling. 0.8 means we keep 80% of cells for each cell type. Default is 1.
#' @param cluster_size_max Maximum number of cells to keep for each cell type during subsampling. Default is 1000.
#' @param cluster_size_min Minimum number of cells to keep for each cell type during subsampling. Default is 1.
#' @param nCV Number of cross validation to perform on the subsampled data. This is to ensure that we have enough cells for cross validation. Default is 1.
#'
#' @return A subsampled matrix
#' @export
#'
subsample=function(cell_cluster_conversion, class_label, rate = 1, cluster_size_max = 1000, cluster_size_min = 1, nCV = 1){
  if (cluster_size_max<=cluster_size_min){
    stop("cluster_size_max should be greater than cluster_size_min")
  }
  if (nCV > cluster_size_min){
    stop("cluster_size_min should be greater or equal to the number of cross validations")
  }
  if (rate > 1){
    stop("rate cannot be greater than 1")
  }
  #if requirements are met, we should have cluster_size_max > cluster_size_min >= nCV. Also we have n<=n.all.cell. As to the relationship between n.all.cell with cluster_size_max, cluster_size_min, and nCV, we don't know.

  #list of cells in the current cluster
  #original_cell_list=as.character(rownames(cell_cluster_conversion)[which(cell_cluster_conversion$class_label == class_label)])       #This works when the cells in cell_cluster_conversion are the same set of cells in count_table
  original_cell_list=which(cell_cluster_conversion$class_label == class_label)
  #cellname=colnames(count_table)
  #original_cell_list=as.character(cellname[which(cell_cluster_conversion[cellname,"class_label"] == class_label)])

  #the number of available cells in the current cluster
  n.all.cell = length(original_cell_list)
  #the number of cells to select after subsampling
  n=ceiling(n.all.cell*rate)

  #############
  #Subsampling#
  #############
  if (n>=nCV){            #if the number of cells after subsampling is greater than the number of cross validation (we need to have at least one cell for each cross validation given a cluster)
    #if the number of cells after subsampling is smaller than cluster_size_min, we select cluster_size_min cells instead of n cells
    if (n<=cluster_size_min){
      if (cluster_size_min > n.all.cell){         #if cluster_size_min is larger than the total number of available cells in the current cluster
        selected_cell_list=c(original_cell_list, original_cell_list[base::sample(seq(1:n.all.cell), cluster_size_min - n.all.cell, replace = TRUE)])     #we sample with replacement (select all available cells first, for the remaining part, sample with replacement)
      }
      if (cluster_size_min <= n.all.cell){        #if cluster_size_min is not larger than the total number of available cells in the current cluster
        selected_cell_list=original_cell_list[base::sample(seq(1:n.all.cell), cluster_size_min, replace = FALSE)]    #we sample without replacement
      }
    }
    #if the number of cells after subsampling is larger than cluster_size_min && smaller than cluster_size_max, we select n cells
    if (n>cluster_size_min && n<=cluster_size_max){
      selected_cell_list=original_cell_list[base::sample(seq(1:n.all.cell), n, replace = FALSE)]    #n is always <= n.all.cell, so we sample without replacement
    }
    #if the number of cells after subsampling is larger than cluster_size_max, we select cluster_size_max cells instead of n cells
    if (n>cluster_size_max){
      selected_cell_list=original_cell_list[base::sample(seq(1:n.all.cell), cluster_size_max, replace = FALSE)]    #since n<= n.all.cell, and in this case cluster_size_max < n, so we have cluster_size_max < n.all.cell, so we sample without replacement
    }
  }else{                  #if the number of cells after subsamping is lower than the number of cross validations, it means it is also smaller than cluster_size_min. So we select cluster_size_min cells instead of n cells
    if (cluster_size_min > n.all.cell){
      selected_cell_list=c(original_cell_list, original_cell_list[base::sample(seq(1:n.all.cell), cluster_size_min - n.all.cell, replace = TRUE)])    #we sample with replacement
    }
    if (cluster_size_min <= n.all.cell){
      selected_cell_list=original_cell_list[base::sample(seq(1:n.all.cell), cluster_size_min, replace = FALSE)]    #we sample without replacement
    }
  }
  return(selected_cell_list)
}


#' Replicate a vector row-wise to generate a matrix
#'
#' @param x A numeric vector.
#' @param n Number of times to repeat \code{x}.
#'
#' @return A numeric matrix.
#' @export
#'
#' @examples
#' x = c(1:10)
#' rep.row(x, 3)
rep.row=function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


#' Replicate a vector column-wise to generate a matrix
#'
#' @param x A numeric vector.
#' @param n Number of times to repeat \code{x}.
#'
#' @return A numeric matrix.
#' @export
#'
#' @examples
#' x = c(1:10)
#' rep.col(x, 3)
rep.col=function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


#' Generate random values from a zero-inflated negative binomial distribution.
#'
#' @param n Number of random values to return.
#' @param size The dispersion paramater used in ‘dnbinom’.
#' @param pie The zero-inflation parameter.
#' @param mu Alternative parametrization of negative binomial distribution via mean. The mean parameter used in ‘dnbinom’.
#'
#' @return A numeric vector containing the randomly generated values.
#' @export
#'
#' @examples
#' zinb_generator(30, size = rep(2, 30), pie = 0.1, mu = rep(5, 30))
zinb_generator=function(n, size, pie, mu){
  counts=1-stats::rbinom(n=n, size=1, prob=pie)
  #for the rest, randomly generate them from a negative binomial distribution (second component of ZINB)
  bool = counts==1
  n.nb=sum(bool)
  counts[which(bool)]=stats::rnbinom(n=n.nb, mu = mu[which(bool)], size = size[which(bool)])
  return(counts)
}

#' Calculate weighted accuracy based on confusion matrix and weighted penalty matrix
#'
#' @param confusion_matrix Confusion matrix
#' @param metric A character specifying the metric to use for evaluating the gene panel's classification performance.
#' Default is "Accuracy", which is the overall accuracy of classification.
#' The other options is "Kappa", which is the Kappa statistics.
#' @param weight_penalty A weighted penalty matrix specifying the partial credit and extra penalty for correct and incorrect classifications between pairs of cell types.
#' It should be a square matrix with cell types as both row and column name.
#'
#' @return A list with elements:
#'   \item{weighted.confusion.matrix}{Weighted confusion matrix.}
#'   \item{weighted.metric}{Weighted fitness value.}
#' @export
#'
weighted_fitness=function(confusion_matrix, metric = "Accuracy", weight_penalty){
  reorder.weight_penalty = weight_penalty[rownames(confusion_matrix), colnames(confusion_matrix)]
  weighted.confusion.matrix = confusion_matrix*reorder.weight_penalty
  if (metric=="Accuracy"){
    weighted.metric = sum(diag(weighted.confusion.matrix))/sum(weighted.confusion.matrix)
  }else{
    weighted.metric=NULL
  }
  return(list(weighted.confusion.matrix=weighted.confusion.matrix, weighted.metric=weighted.metric))
}
