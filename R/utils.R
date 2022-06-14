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


#' A scatterplot with color corresponding to density
#'
#' @description Visualize two dimensional data using scatterplot with density information
#' @param x A numeric vector containing the \code{x} coordinates of points in the plot
#' @param y A numeric vector containing the \code{y} coordinates of points in the plot
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param main an overall title for the plot
#' @param add.diagonal a logical value indicating whether to add \code{x = y} line to the plot
#' @param x_low A numeric value representing the lower limit of \code{x} coordinate
#' @param x_high A numeric value representing the upper limit of \code{x} coordinate
#' @param y_low A numeric value representing the lower limit of \code{y} coordinate
#' @param y_high A numeric value representing the upper limit of \code{y} coordinate
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
#' @param base Base of log transformation, e.g., 10
#'
#' @details For log and logit transofmration, if the value before transformation is 0, we will replace it by a small value before the transformation to avoid negative infinity. The small value is calculated as the minimum of all non-zero values from the input divided by 2
#' @return A numeric vector of transformed values
#' @export
#'
#' @examples
#' x = sample(1:100, size = 10)
#' data_transformation(x, trans_type = "log", base=10)
data_transformation=function(input_data ,trans_type, base=10){
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




