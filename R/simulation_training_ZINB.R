#' Platform effect estimation using Bayesian modelling
#'
#' @description simulation_training_ZINB estimates platform effects between scRNA-seq and spatial transcriptomics technologies using Bayesian modelling. The Bayesian model is fitted using rstan
#' @param sc_count A count matrix containing the expression of each gene in each cell from the scRNA-seq data. Each row represents one gene and each column represents one column, with gene name as row name and cell name as column name.
#' @param spatial_count A count matrix containing the expression of each gene in each cell from the spatial transcriptomics data. Each row represents one gene and each column represents one column, with gene name as row name and cell name as column name.
#' @param overlap_gene A character vector of gene names that show up in both scRNA-seq and spatial transcriptomics data
#' @param unique_cluster_label A character vector of cell type names that show up in both scRNA-seq and spatial transcriptomics data
#' @param sc_cluster A data frame with two columns and cell name as row name. The first column contains the name of cells from scRNA-seq data and the second column contains each cell's corresponding cell type.
#' @param spatial_cluster A data frame with two columns and cell name as row name. The first column contains the name of cells from spatial transcriptomics data and the second column contains each cell's corresponding cell type.
#' @param outputpath A character specifying the path for output
#' @param optimizer Optimizer for Bayesian model fitting.
#' * \code{variational_inference}: The Bayesian model will be fitted using the \code{meanfield} variational inference
#' * \code{sampling}: The Bayesian model will be fitted using the No-U-Turn sampler variant of Hamiltonian Monte Carlo
#' @param mcmc.check a logical value indicating whether to do MCMC diagnostics. Keep it as FALSE if \code{optimizer} is \code{variational_inference}
#' @param saveplot a logical value indicating whether to save plot
#' @param num.iter A positive integer specifying the number of iterations for each chain (including warmup). The default is 2000. It is recommended to use at least 10000 if \code{optimizer} is \code{variational_inference}
#' @param num.chain A positive integer specifying the number of Markov chains. The default is 4.
#' @param num.core The number of cores to use when executing the Markov chains in parallel. The default is to 1 core. However, we recommend setting it to be as many processors as the hardware and RAM allow (up to the number of chains).
#' @param max.treedepth A positive integer specifying the maximum tree depth. The default is 10.
#' @param seed The seed for random number generation, defaulting to 3.
#'
#' @return A list with the following elements:
#' \item{distortion_model}{A stanfit object containing the fitted Bayesian model.}
#' \item{data2fit}{A data frame containing the input data for fitting the Bayesian model:}
#' * \code{sp.count}: number of molecule per gene per cell from the spatial transcriptomis data
#' * \code{sc.prop}: relative expression per gene per cell from the scRNA-seq data
#' * \code{gene}: gene name
#' * \code{cell}: cell name
#' * \code{cluster}: cell type name of cells in the spatial transcriptomics data
#' * \code{cs}: cell size, i.e., total number of molecule per cell, of cells in the spatial transcriptomics data
#'
#' \item{model.summary}{A named list with elements \code{summary} and \code{c_summary}, which contain summaries of a stanfit object.}
#' \item{fixed.effect}{A named list with elements \code{summary} and \code{c_summary}, which contain summaries of specified parameters.}
#' \item{c_i_full}{A matrix containing the estimated gene specific intercept \eqn{c_i} for each gene.}
#' \item{gamma_i_full}{A matrix containing the estimated gene specific coefficient \eqn{\gamma_i} for each gene.}
#' \item{lib.size}{A numeric vector containing the library size, i.e., total number of molecule per cell, of cells in the spatial transcriptomics data.}
#'
#' @export
#'
#' @examples
#' data(sc_count)
#' data(sc_cluster)
#' data(spatial_count)
#' data(spatial_cluster)
#' overlap_gene = intersect(unique(rownames(sc_count)),unique(rownames(spatial_count)))
#' unique_cluster_label=intersect(as.character(unique(sc_cluster$class_label)),
#'                                as.character(unique(spatial_cluster$class_label)))
#' outputpath = "~"
#' simulation.params=simulation_training_ZINB(sc_count = sc_count,
#'                                            spatial_count = spatial_count,
#'                                            overlap_gene = overlap_gene,
#'                                            unique_cluster_label = unique_cluster_label,
#'                                            sc_cluster = sc_cluster,
#'                                            spatial_cluster = spatial_cluster,
#'                                            outputpath = outputpath,
#'                                            optimizer = "variational_inference",
#'                                            mcmc.check = "FALSE",
#'                                            num.iter = 300,
#'                                            num.chain = 4,
#'                                            num.core = 4,
#'                                            max.treedepth = 10,
#'                                            seed = 3,
#'                                            saveplot = FALSE)
#'
#' @seealso \code{\link[rstan]{stan}},  \code{\link[rstan]{vb}} from the \code{rstan} package
simulation_training_ZINB=function(sc_count, spatial_count,
                                  overlap_gene, unique_cluster_label, sc_cluster, spatial_cluster,
                                  outputpath, optimizer, mcmc.check = FALSE, saveplot=FALSE,
                                  num.iter = 2000, num.chain = 4, num.core = 1, max.treedepth = 10, seed = 3){
  current.dir = getwd()

  #################
  #1. prepare data#
  #################
  print("Prepare data")
  gene.list=overlap_gene
  cell.list=colnames(spatial_count)
  num.gene=length(gene.list)
  num.cell=length(cell.list)

  ###spatial data###
  #count from spatial data for each gene in each cell
  sp.count.matrix=as.matrix(spatial_count[gene.list, cell.list])

  #cell size for each cell in spatial data
  cell.size=colSums(spatial_count[, cell.list])
  cell.size.matrix=t(replicate(num.gene, cell.size))
  rownames(cell.size.matrix)=gene.list

  ###scRNA-seq data###
  #get the relative expression for each gene in each cluster from scRNA-seq data (denominator is count from all genes)
  prop_ik_sc_all_gene=sapply(unique_cluster_label, gene_list=gene.list, relative_freq,
                             count_matrix=sc_count,
                             cell_cluster_conversion=sc_cluster)

  #relative expression from scRNA-seq data for each gene in each cell
  sp.cell.type=as.character(spatial_cluster[cell.list, "class_label"])               #cell type of each cell in the spatial data
  sc.prop.matrix=as.matrix(prop_ik_sc_all_gene[gene.list, sp.cell.type])

  #gene name for each gene in each cell
  gene.name.matrix=as.matrix(replicate(num.cell, gene.list))
  rownames(gene.name.matrix)=gene.list

  #cluster label for each gene in each cell
  cluster.label.matrix=as.matrix(t(replicate(num.gene, sp.cell.type)))
  rownames(cluster.label.matrix)=gene.list

  cell.name.matrix=as.matrix(t(replicate(num.gene, cell.list)))
  rownames(cell.name.matrix)=gene.list

  genes.2.keep=gene.list

  sp.count=as.numeric(sp.count.matrix[genes.2.keep,])
  cell.size=as.numeric(cell.size.matrix[genes.2.keep,])
  sc.prop=as.numeric(sc.prop.matrix[genes.2.keep,])
  gene.name=as.character(gene.name.matrix[genes.2.keep,])
  cluster.label=as.character(cluster.label.matrix[genes.2.keep,])
  cell.name=as.character(cell.name.matrix[genes.2.keep,])

  data2fit.all = data.frame(sp.count = sp.count,
                            sc.prop = sc.prop,
                            gene = gene.name, cell = cell.name,
                            cluster = cluster.label, cs = cell.size)

  #preparing data for Stan
  #convert gene.name to numeric value (random effect has to be numeric)
  conv.table = data.frame(gene.name = unique(gene.name),
                          gene.number = 1:length(unique(gene.name)))
  rownames(conv.table) = conv.table$gene.name

  data2fit.stan = list(y = data2fit.all$sp.count,
                       x = data2fit.all$sc.prop,
                       N = dim(data2fit.all)[1],
                       Nrandom = length(unique(data2fit.all$gene)),
                       raneff = conv.table[data2fit.all$gene, "gene.number"],
                       cell = data2fit.all$cell,
                       cluster = data2fit.all$cluster,
                       cs = data2fit.all$cs)

  ############################
  #basic data inspection plot#
  ############################
  #calculate relative expression for scRNA-seq and spatial data
  #scRNA-seq data:
  sub_sc_count=sc_count[gene.list,]
  prop_ik_sc=sapply(unique_cluster_label, gene_list=gene.list, relative_freq,
                    count_matrix=sub_sc_count,
                    cell_cluster_conversion=sc_cluster)
  #spatial data:
  sub_spatial_count=spatial_count[gene.list,]
  prop_ik_spatial=sapply(unique_cluster_label, gene_list=gene.list, relative_freq,
                         count_matrix=sub_spatial_count,
                         cell_cluster_conversion=spatial_cluster)
  if (saveplot){
    print("Generating basic data inspection plot")
    setwd(outputpath)

    grDevices::pdf("Relationship between relative expression in scRNA-seq and relative expression in spatial data (based on overlapping genes).pdf", height=8, width=8)
    densityscatter(x=as.numeric(as.matrix(prop_ik_sc)), y=as.numeric(as.matrix(prop_ik_spatial)),
                   xlab="relative expression in snRNA-seq (original)", ylab="relative expression in spatial data (original)", main = "original spatial data (overlapping genes)")
    densityscatter(x=data_transformation(as.numeric(as.matrix(prop_ik_sc)), "log", base=exp(1)), y=data_transformation(as.numeric(as.matrix(prop_ik_spatial)), "log", base=exp(1)),
                   xlab="relative expression in snRNA-seq (log)", ylab="relative expression in spatial data (log)", main = "original spatial data (overlapping genes)")
    grDevices::dev.off()
    grDevices::pdf("Relationship between relative expression in scRNA-seq and relative expression in spatial data (based on overlapping genes) per gene.pdf", height=8, width=8)
    for (i in 1:dim(prop_ik_sc)[1]){
      densityscatter(x=data_transformation(as.numeric(as.matrix(prop_ik_sc[i,])), "log", base=exp(1)), y=data_transformation(as.numeric(as.matrix(prop_ik_spatial[i,])), "log", base=exp(1)),
                     xlab="relative expression in snRNA-seq (log)", ylab="relative expression in spatial data (log)", main = paste("original spatial data (overlapping genes)", rownames(prop_ik_sc)[i]))
    }
    grDevices::dev.off()

    ###############################
    #check the level of distortion#
    ###############################
    #Estimate distortion based on differential expression
    prop_matrix=cbind(prop_ik_sc,prop_ik_spatial)
    group_label=c(rep("sc",dim(prop_ik_sc)[2]), rep("spatial",dim(prop_ik_spatial)[2]))
    distortion=t(apply(prop_matrix, 1, distortion_test, group_label=group_label))
    colnames(distortion)=c("beta.regular", "beta.deming")

    grDevices::pdf("Distortion based on differential expression.pdf")
    graphics::hist(distortion[,"beta.regular"], xlab = "regression coefficient", main = "regression coefficient (linear regression)",
         breaks = seq(0, ceiling(max(distortion[,"beta.regular"])), by = 1), xlim = range(c(0, ceiling(max(distortion[,"beta.regular"])))))
    graphics::hist(distortion[,"beta.deming"], xlab = "regression coefficient", main = "regression coefficient (Deming regression)",
         breaks = seq(0, ceiling(max(distortion[,"beta.deming"], na.rm=T)), by = 1), xlim = range(c(0, ceiling(max(distortion[,"beta.deming"], na.rm=T)))))

    #hist(log(distortion[,"fc"]), xlab = "Fold change (log)", main = "Fold change (log)")
    graphics::hist(log(distortion[,"beta.regular"]), xlab = "Linear regression coefficient (log)", main = "Linear regression coefficient (log)")
    graphics::hist(log(distortion[,"beta.deming"]), xlab = "Deming regression coefficient (log)", main = "Deming regression coefficient (log)")

    #plot(density(log(distortion[,"fc"])), xlab = "Fold change (log)", main = "Fold change (log)")
    plot(stats::density(log(distortion[,"beta.regular"])), xlab = "Linear regression coefficient (log)", main = "Linear regression coefficient (log)")
    plot(stats::density(log(distortion[,"beta.deming"][!is.na(distortion[,"beta.deming"])])), xlab = "Deming regression coefficient (log)", main = "Deming regression coefficient (log)")
    grDevices::dev.off()
  }


  ##############
  #2. fit model#
  ##############
  print("Start model fitting")

  #specify model fitting parameters
  num.iter = num.iter
  num.chain = num.chain
  num.core = num.core
  max.treedepth = max.treedepth

  #prepare the Bayesian model
  model_code = '
        data {
            int<lower=1> N;                       // sample size (number of entries in the input data for modeling fitting)
            int<lower=1> Nrandom;                  //number of random effect level (for gene specific model, this equals to the number of genes)
            int<lower=0> y[N];                    //counts (spatial counts)
            real<lower=0, upper=1> x[N];          //raw relative proportion
            vector<lower=0>[N] cs;                   //raw cell size for each cell
            int<lower=0, upper=Nrandom> raneff[N];   // index of random effect (mapping between data point (from 1 to N) to random effect (from 1 to Nrandom))
          }

        transformed data {
          vector[N] log_cs = log(cs);
        }

        parameters {
            real alpha_tilde;                        // coef of linear pred for miu
            real beta_tilde;                        // coef of linear pred for theta
            real<lower=0,upper=1> zi;         // zero inflation
            real<lower=0> sigma_alpha;        // standard deviation of prior distribution of alpha
            real<lower=0> sigma_beta;         // standard deviation of prior distribution of beta

            vector<lower=0>[Nrandom] gamma_i;
            real mu_gamma;
            real<lower=0> sigma_gamma;
            vector[Nrandom] c_i;
            real mu_c;
            real<lower=0> sigma_c;
          }

        transformed parameters {
          real alpha;
          real beta;

          alpha = sigma_alpha * alpha_tilde;
          beta = sigma_beta * beta_tilde;
        }

        model {
          vector[N] theta;                  // dispersion of ZINB distribution
          vector[N] miu;                    // mean of ZINB distribution

          vector[N] lambda;
          vector[N] log_lambda;
          for (i in 1:N){
            lambda[i] = inv_logit(gamma_i[raneff[i]]*sqrt(x[i]) + c_i[raneff[i]]);
            log_lambda[i] = log(lambda[i]);
          }

          theta = exp(beta + lambda);

          miu = exp(alpha + log_lambda + log_cs);

          gamma_i ~ lognormal(mu_gamma, sigma_gamma);
          c_i ~ normal(mu_c, sigma_c);

          // priors
          alpha_tilde ~ normal(0,1);
          beta_tilde ~ normal(0,1);
          sigma_alpha ~ cauchy(0,5);
          sigma_beta ~ cauchy(0,5);
          sigma_gamma ~ cauchy(0,5);
          mu_gamma ~ cauchy(0,5);
          mu_c ~ cauchy(0,5);
          sigma_c ~ cauchy(0,5);
          zi ~ beta(1,1);

          // likelihood
          for (i in 1:N){
            if(y[i] == 0){
              target += log_sum_exp(bernoulli_lpmf(1 | zi),               //there is a twist here. zi is the probability of generating 1 from bernoulli_rng but here, we use this probability of generating 1 as the probability of seeing an inflated 0, i.e., y=0. In another word, generating 1 from the bernoulli distribution means have 0 as y
                                    bernoulli_lpmf(0 | zi) +
                                    neg_binomial_2_lpmf(y[i] | miu[i], theta[i]));
            }else{
              target += bernoulli_lpmf(0 | zi) +
                        neg_binomial_2_lpmf(y[i] | miu[i], theta[i]);
            }
          }
        }

        generated quantities {
          vector[N] y_rep;
          vector[N] lambda;
          vector[N] log_lambda;

          vector[N] theta;                  // dispersion of ZINB distribution
          vector[N] miu;                    // mean of ZINB distribution
          vector[N] log_lik;

          // randomly draw count from ZINB distribution
          for(i in 1:N){
            lambda[i] = inv_logit(gamma_i[raneff[i]]*sqrt(x[i]) + c_i[raneff[i]]);
            log_lambda[i] = log(lambda[i]);

            if (bernoulli_rng(zi)==1){
              y_rep[i] = 0;
            }else{
              y_rep[i] = neg_binomial_2_rng(exp(alpha + log_lambda[i] + log_cs[i]), exp(beta + lambda[i]));
            }

            //calculate log_lik
            miu[i] = exp(alpha + log_lambda[i] + log_cs[i]);
            theta[i] = exp(beta + lambda[i]);
            if(y[i] == 0){
              log_lik[i] = log_sum_exp(bernoulli_lpmf(1 | zi),               //there is a twist here. zi is the probability of generating 1 from bernoulli_rng but here, we use this probability of generating 1 as the probability of seeing an inflated 0, i.e., y=0. In another word, generating 1 from the bernoulli distribution means have 0 as y
                                    bernoulli_lpmf(0 | zi) +
                                    neg_binomial_2_lpmf(y[i] | miu[i], theta[i]));
            }else{
              log_lik[i] = bernoulli_lpmf(0 | zi) +
                            neg_binomial_2_lpmf(y[i] | miu[i], theta[i]);
            }
          }
        }
      '

  # model fit
  if (optimizer == "sampling"){
    fit.m <- rstan::stan(model_code = model_code,
                         data=data2fit.stan,
                         iter=num.iter, chains=num.chain, cores = num.core,
                         seed = seed,
                         save_warmup = TRUE,
                         control = list(max_treedepth = max.treedepth))
  }

  if (optimizer == "variational_inference"){
    mod = rstan::stan_model(model_code = model_code)
    fit.m <- rstan::vb(mod, data = data2fit.stan,
                       seed = seed,
                       iter=num.iter, tol_rel_obj = 0.00000001)
  }

  model.summary = rstan::summary(fit.m)

  ######################
  #3. model fit summary#
  ######################
  setwd(outputpath)
  params2check = c("alpha_tilde", "beta_tilde", "zi", "sigma_alpha", "sigma_beta", "mu_gamma", "sigma_gamma", "mu_c", "sigma_c", "alpha", "beta")

  c_i_full = rstan::summary(fit.m, par=c("c_i"))$summary
  rownames(c_i_full) = rownames(conv.table)

  gamma_i_full = rstan::summary(fit.m, par=c("gamma_i"))$summary
  rownames(gamma_i_full) = rownames(conv.table)

  fixed.effect = rstan::summary(fit.m, par=params2check)

  if (saveplot){
    print("Generating model fit summary plot")

    ###basic summary###
    grDevices::pdf("model fit summary.pdf")
    p=rstan::stan_dens(fit.m, pars = params2check)
    print(p)
    grDevices::dev.off()

    ###trace plot###
    grDevices::pdf("trace plot.pdf", height = 16, width = 16)
    p=rstan::traceplot(fit.m, pars = params2check)
    print(p)
    grDevices::dev.off()

    ###Posterior predictive check###
    #In order to use the PPC functions from the bayesplot package we need a vector y of outcome values,
    y <- data2fit.stan$y
    #and a matrix yrep of draws from the posterior predictive distribution,
    yrep_zinb <- rstan::extract(fit.m, pars = "y_rep")[[1]]  #each row is one set of data from posterior predictive distribution

    nsamples = 9

    #1. The first PPC we’ll look at is a comparison of the distribution of y and the distributions of some of the simulated datasets (rows) in the yrep matrix.
    bayesplot::color_scheme_set("brightblue")
    grDevices::pdf("Posterior predictive check - ppc_dens_overlay.pdf")
    p1=bayesplot::ppc_dens_overlay(y, yrep_zinb[1:nsamples, ])
    #In the plot above, the dark line is the distribution of the observed outcomes y and each of the 50 lighter lines is the kernel density estimate of one of the replications of y from the posterior predictive distribution (i.e., one of the rows in yrep).
    #To see the discrepancy at the lower values of more clearly we can use the xlim function from ggplot2 to restrict the range of the x-axis:
    p2=bayesplot::ppc_dens_overlay(y, yrep_zinb[1:nsamples, ]) + ggplot2::xlim(0, 10)
    p3=bayesplot::ppc_dens_overlay(y, yrep_zinb[1:nsamples, ]) + ggplot2::xlim(0, 2.5)
    print(p1)
    print(p2)
    print(p3)
    grDevices::dev.off()

    grDevices::png("Posterior predictive check - qqplot.png", height = 1000, width = 1000)
    graphics::par(mfrow=c(3,3))
    for (i in 1:nsamples){
      stats::qqplot(y, yrep_zinb[i,], xlab="observed y", ylab=paste("simulated y", i), main="QQ plot")
      graphics::abline(0,1)
    }
    grDevices::dev.off()

    #2. check the distribution of test statistics, e.g., estimated 0s (are we generating enough 0? Too many 0 or too few?)
    #First we define a function that takes a vector as input and returns the proportion of zeros:
    prop_zero <- function(x) mean(x == 0)
    prop_zero(y) # check proportion of zeros in y
    #The stat argument to ppc_stat accepts a function or the name of a function for computing a test statistic from a vector of data. In our case we can specify stat = "prop_zero"
    grDevices::pdf("Posterior predictive check - ppc_stat.pdf")
    p=bayesplot::ppc_stat(y, yrep_zinb, stat = "prop_zero", binwidth = 0.005) + ggplot2::ggtitle("Proportion of zero")
    print(p)
    grDevices::dev.off()
    #The dark line is at the value T(y), i.e. the value of the test statistic computed from the observed y, in this case prop_zero(y).
    #The lighter area on the left is actually a histogram of the proportion of zeros in in the yrep simulations

    #3. for each simulation, check the percentage of 0 count per gene between observed data and simulated data
    N=length(genes.2.keep)
    M=dim(sp.count.matrix)[2]
    grDevices::pdf("Percentage of 0 per gene between original data and simulated data.pdf")
    for (i in 1:nsamples){
      data2fit.all$pred=yrep_zinb[i,]
      #Do we have probe failure in this simulation? Let's first calculate the percentage of 0 counts per gene
      num.0 = num.0.simu = rep(0, N)
      for (j in 1:N){           #for each gene, we calculate the number of 0s across all cells
        data2fit.per.gene = subset(data2fit.all, data2fit.all$gene==genes.2.keep[j])
        num.0[j]=sum(data2fit.per.gene$sp.count==0)
        num.0.simu[j]=sum(data2fit.per.gene$pred==0)
      }
      perc.0 = num.0/M
      perc.0.simu = num.0.simu/M
      densityscatter(perc.0, perc.0.simu, xlab="percentage of 0 count per gene in observed data", ylab="percentage of 0 count per gene in simulated data", main = paste("simulation", i), add.diagonal = T)
    }
    grDevices::dev.off()


    ###agreement between scRNA-seq data and spatial data check###
    grDevices::pdf("Relationship between relative expression in scRNA-seq and simulated spatial data.pdf", height=8, width=8)
    for (s in 1:nsamples){
      data2fit.all$simu.sp.count = yrep_zinb[s,]               #we choose one of the simulations as the simulated spatial count

      #change the data from long format back to wide format
      temp=data2fit.all[,c("simu.sp.count", "gene", "cell")]
      simu_spatial_count = stats::reshape(temp, idvar = "gene", timevar = "cell", direction = "wide")
      rownames(simu_spatial_count)=simu_spatial_count$gene
      simu_spatial_count=simu_spatial_count[,-1]
      colnames(simu_spatial_count)=cell.list

      # simu_spatial_count[is.na(simu_spatial_count)]=0
      #for cells from cell types with relative proportion=0, their count here will be NA because these entries are filtered out from data2fit.all
      # saveRDS(simu_spatial_count, paste0("simu_spatial_count_", s,".rds"))

      #get the average expression for each gene in each cluster from spatial data
      # ave_count_ik_sp_gene2keep_simu=sapply(unique_cluster_label, gene_list=genes.2.keep, average_count,
      #                                       count_matrix=simu_spatial_count,
      #                                       cell_cluster_conversion=spatial_cluster)
      # #get the relative expression for each gene in each cluster from spatial data
      prop_ik_sp_gene2keep_simu=sapply(unique_cluster_label, gene_list=genes.2.keep, relative_freq,
                                       count_matrix=simu_spatial_count,
                                       cell_cluster_conversion=spatial_cluster)

      #(Plot) Relative expression in scRNA-seq data vs. relative expression expression in simulated spatial data
      sc=as.numeric(as.matrix(prop_ik_sc_all_gene[genes.2.keep,]))
      # sp=as.numeric(as.matrix(ave_count_ik_sp_gene2keep_simu))

      # densityscatter(x=data_transformation(sc, "log", base=exp(1)), y=data_transformation(sp, "log", base=exp(1)),
      #                xlab="relative expression in snRNA-seq (log)", ylab="average expression in spatial data (log)", main = paste("simulated spatial data", s))
      densityscatter(x=data_transformation(as.numeric(as.matrix(prop_ik_sc[genes.2.keep,])), "log", base=exp(1)), y=data_transformation(as.numeric(as.matrix(prop_ik_sp_gene2keep_simu)), "log", base=exp(1)),
                     xlab="relative expression in snRNA-seq (log)", ylab="relative expression in spatial data (log)", main = paste("simulated spatial data", s, "(overlapping genes)"))
    }
    grDevices::dev.off()
  }

  ###MCMC summary plot###
  if (mcmc.check){
    print("Start MCMC check")

    #specify parameters to check
    params.list=params.interest=params2check
    transformation = list("zi" = "logit", "sigma_alpha" = "log", "sigma_beta" = "log",  "sigma_gamma" = "log", "sigma_c" = "log")

    #Posterior uncertainty intervals
    grDevices::pdf("MCMC draws - Posterior uncertainty intervals.pdf")
    p=bayesplot::mcmc_intervals(fit.m, pars = params.list)
    print(p)
    grDevices::dev.off()

    #Bivariate plots
    grDevices::png("MCMC draws - Bivariate plots.png", height = 2000, width = 2000)
    p1=bayesplot::mcmc_pairs(fit.m, pars = params.list,
                             off_diag_args = list(size = 1.5))
    p2=bayesplot::mcmc_pairs(fit.m, pars = params.list,
                             transform = transformation,
                             off_diag_args = list(size = 1.5))
    print(p1)
    print(p2)
    grDevices::dev.off()

    #autocorrelation
    grDevices::pdf("MCMC diagnostics - mmcmc_acf (autocorrelation).pdf", height = 16, width = 16)
    p=bayesplot::mcmc_acf(fit.m, pars = params.interest, lags = 100)
    print(p)
    grDevices::dev.off()

    if (optimizer == "sampling"){
      #Diagnostics for the No-U-Turn Sampler
      lp <- bayesplot::log_posterior(fit.m)
      np <- bayesplot::nuts_params(fit.m)

      #Energy and Bayesian fraction of missing information
      setwd(outputpath)
      grDevices::pdf("MCMC diagnostics - mcmc_nuts_energy.pdf")
      bayesplot::color_scheme_set("red")
      p=bayesplot::mcmc_nuts_energy(np)
      print(p)
      grDevices::dev.off()

      #General MCMC diagnostics
      rhats <- bayesplot::rhat(fit.m)
      grDevices::pdf("MCMC diagnostics - mcmc_rhat.pdf")
      bayesplot::color_scheme_set("brightblue") # see help("color_scheme_set")
      p1=bayesplot::mcmc_rhat(rhats)
      p2=bayesplot::mcmc_rhat_hist(rhats)
      print(p1)
      print(p2)
      grDevices::dev.off()

      #Effective sample size
      ratios <- bayesplot::neff_ratio(fit.m)      #get the neff value
      grDevices::pdf("MCMC diagnostics - mcmc_neff (Effective sample size).pdf")
      p=bayesplot::mcmc_neff(ratios, size = 2)
      print(p)
      grDevices::dev.off()
    }
  }

  ###########
  #4. output#
  ###########
  print("Return result")
  lib.size=colSums(spatial_count[, cell.list])

  simulation.model=list(distortion_model = fit.m,
                        data2fit = data2fit.all,
                        model.summary = model.summary,
                        fixed.effect = fixed.effect,
                        c_i_full = c_i_full,
                        gamma_i_full = gamma_i_full,
                        lib.size = lib.size)
  setwd(current.dir)
  return(simulation.model)
}
