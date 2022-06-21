#' Genetic algorithm for gene panel selection
#'
#' @description gpsFISH_optimize implements a genetic algorithm for gene panel selection. The function searches for a subset
#' of a fixed size, \code{k} gene, from the \code{n} candidate gene, such that user-supplied function \code{OF} is minimized at that subset.
#' The selection step is done by tournament selection based on ranks, and elitism may be used to retain the best solutions from one generation to the next.
#' The crossover step is done by uniform crossover.
#' Population objective function values can be evaluated in parallel.
#'
#' @param n Total number of genes for gene panel selection.
#' @param k Size of a gene panel, i.e., the number of genes in a gene panel.
#' @param OF Fitness function, a function that evaluate the fitness of each gene panel. Default is \code{fitness}.
#' Additional arguments can be passed to \code{OF} through \code{\dots}. See \link[gpsFISH]{fitness} for more details.
#' @param popsize Size of the population, i.e., the number of gene panels in a population. Default is 200.
#' @param keepbest The \code{keepbest} least fit offspring of each generation are replaced by the \code{keepbest} most fit members of the previous generation. Used to implement elitism. Default is calculated based on \code{popsize}.
#' @param ngen The number of generations to run the optimization. Default is 50.
#' @param tourneysize The number of individuals involved in each tournament during selection. Default is calculated based on \code{popsize}.
#' @param mutprob Mutation probability of each gene in a gene panel. This probability can be set indirectly through \code{mutfrac}. Default is 0.01.
#' @param mutfrac The average fraction of offspring that will experience at least one mutation. Only used if mutprob is not supplied.
#' @param initpop A numeric matrix specifying the initial population with each row representing one gene panel and each column representing one gene in a gene panel.
#' The genes are encoded by their location in the column name of \code{full_count_table}.
#' If not provided (default is NULL), it will be randomly initialized from all the candidate genes.
#' The final populations from one gene panel selection can be passed as the initial population of the next selection.
#' Possibly useful if using this function in an adaptive, iterative, or parallel scheme.
#' @param verbose An integer controlling the display of progress during optimization.
#' If \code{verbose} is greater than 0, then the iteration number and best fitness value are displayed at the console every generation.
#' Otherwise nothing is displayed. Default is zero (no display).
#' @param cluster The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one, and parallelization requires at least two cores.
#' @param save.intermediate A logical value specifying whether intermediate populations are outputted to the current working directory. Default is FALSE.
#' @param gene2include.id A numeric vector specifying the location of genes in the column name of \code{full_count_table} that must be included in each panel of the population. Default is NULL.
#' @param gene.weight A data frame specifying the weight for each gene. Default is NULL.
#' Each row is a gene. The first column contains the gene name. The second column contains the gene weight.
#' Row name of the data frame is gene name.
#' A gene with higher weight will be (1) more likely to be selected during crossover, (2) less likely to be mutated if it is already in the population, (3) and more likely to be introduced into the population through mutation if it is not in the current population.
#' @param earlyterm If \code{earlyterm} is provided, then the optimization will stop when it reaches \code{ngen} generations, or the maximum fitness improvement (fitness difference between two consecutive generations) in the past \code{earlyterm} generations are smaller or equal to \code{converge.cutoff}, whichever comes first. Default is \code{ngen}.
#' @param converge.cutoff A numeric value specifying the cutoff of maximum fitness improvement. Default is 0.
#' @param ... Additional arguments passed to \code{OF}.
#'
#' @return A list with the following elements:
#'   \item{old}{A list holding information about the search progress. Its elements are:}
#'   \item{old$best}{A matrix containing the best gene panel for each generation.}
#'   \item{old$obj}{A numeric vector containing the fitness values corresponding to the gene panels
#'     in \code{old$best}.}
#'   \item{old$avg}{A numeric vector containing the average population fitness value for each generation.}
#'   \item{pop}{A matrix containing the final population, row-sorted
#'     in order of increasing objective function. Each row is an index vector representing
#'     one gene panel.}
#'   \item{obj}{The fitness values corresponding to each row of \code{pop}.}
#'   \item{bestgeneration}{Number of the best generation.}
#'   \item{bestsol}{A vector of length \code{k} containing the best solution. This is the optimized gene panel.}
#'   \item{bestobj}{The fitness value for the best solution.}
#'   \item{diversity}{A numeric vector containing the diversity of the population in each generation.}
#'   \item{ave_confusionMatrix}{Average confusion matrix over all gene panels for each generation.}
#'   \item{best_confusionMatrix}{Confusion matrix of the best gene panel for each generation.}
#'   \item{norm_ave_confusionMatrix}{Average normalized confusion matrix over all gene panels for each generation.}
#'   \item{norm_best_confusionMatrix}{Normalized confusion matrix of the best gene panel for each generation.}
#'   \item{stats_byclass}{Average classification statistics over all gene panels for each generation.}
#'   \item{best_stats_byclass}{Classification statistics of the best gene panel for each generation.}
#'   \item{best_pred_prob}{Prediction probability of each cell in each cell type in each cross validation of the best gene panel for each generation.}
#'
#'@details
#' \itemize{
#' \item The selection strategy we used is tournaments. Specifically, randomly selected candidate gene panels face each other 1 vs. 1.
#'       The one with a higher fitness value will be used as parent. In addition, candidate gene panels with higher fitness values will be more likely to be selected in the tournaments.
#' \item After having the parent gene panels, uniform crossover is performed to generate the offspring gene panels.
#'       Duplicated genes after uniform crossover will be replaced by randomly sampled genes in the parent candidate gene panels but not in the offspring gene panel.
#' }
#'
#' @export
#'
#' @examples
gpsFISH_optimize = function (n, k, OF = fitness, popsize = 200, keepbest = floor(popsize/10),
                             ngen = 50, tourneysize = max(ceiling(popsize/10), 2), mutprob = 0.01,
                             mutfrac = NULL, initpop = NULL, verbose = 0, cluster = NULL, save.intermediate = F,
                             gene2include.id = NULL,
                             gene.weight = NULL,
                             earlyterm = ngen,
                             converge.cutoff = 0,
                             ...)
{
  ################
  #Argument check#
  ################
  if (!is.null(mutfrac)) {
    if (base::missing(mutprob)) {
      stopifnot(mutfrac >= 0, mutfrac <= 1)
      mutprob = 1 - (1 - mutfrac)^(1/k)
    }
    else {
      warning("Both mutprob and mutfrac supplied. Ignoring mutfrac.")
    }
  }
  stopifnot(n%%1 == 0, n > 0, n >= k, k%%1 == 0, k > 0, popsize%%1 ==
              0, popsize >= 2, popsize > keepbest, popsize > tourneysize,
            keepbest%%1 == 0, keepbest >= 0, ngen%%1 == 0, ngen >=
              1, tourneysize%%1 == 0, tourneysize >= 2, mutprob >=
              0, mutprob <= 1, all(dim(initpop) == c(popsize, k)),
            verbose%%1 == 0)

  ###############
  #Create object#
  ###############
  indices = 1:n
  if (keepbest > 0) {
    elitespots = 1:keepbest
    newspots = (keepbest + 1):popsize
  }
  fitness.old = vector(mode = "numeric", length = popsize)
  fitness.new = vector(mode = "numeric", length = popsize)
  offspring = matrix(0, nrow = popsize, ncol = k)
  old = list()
  old$best = matrix(0, nrow = ngen + 1, ncol = k)                    #store the best solution (initial population + every generation)
  old$obj = vector(mode = "numeric", length = ngen + 1)              #store the best fitness value (small fitness -> better result)
  old$avg = vector(mode = "numeric", length = ngen + 1)              #average fitness value of the current population
  Tourneys = matrix(0, nrow = popsize, ncol = 2*tourneysize)
  Tourneys1 = matrix(0, nrow = popsize, ncol = tourneysize)
  Tourneys2 = matrix(0, nrow = popsize, ncol = tourneysize)

  #diversity of population#
  diversity = vector(mode = "numeric", length = ngen + 1)            #diversity of the population for each generation

  #confusion matrix of population#
  ave_confusionMatrix = norm_ave_confusionMatrix = vector(mode = "list", length = ngen + 1)     #average over all (normalized) confusion matrix of all chromosomes in a population for each generation
  best_confusionMatrix = norm_best_confusionMatrix = vector(mode = "list", length = ngen + 1)     #(normalized) confusion matrix of the best chromosome in a population for each generation

  #statistics by class from the classifier#
  stats_byclass = vector(mode = "list", length = ngen + 1)                    #average over all classification statistics by class of all chromosomes in a population for each generation
  best_stats_byclass = vector(mode = "list", length = ngen + 1)               #classification statistics by class of the best chromosome in a population for each generation

  #prediction probability and by class AUC#
  best_pred_prob = vector(mode = "list", length = ngen + 1)                    #prediction probability of the best chromosome in a population for each generation
  best_AUC_by_class = vector(mode = "list", length = ngen + 1)                 #classification AUC by class of the best chromosome in a population for each generation

  #calculate feature importance from classifier#
  best_feature_imp = vector(mode = "list", length = ngen + 1)                 #feature importance of the best chromosome in a population for each generation

  ###################################
  #parameters for parallel computing#
  ###################################
  if (cluster>1){
    useparallel = TRUE
  }else{
    useparallel = FALSE
  }

  ####################
  #initial population#
  ####################
  if (is.null(initpop)){                             #randomly generate initial population if it is not provided
    if (!is.null(gene2include.id)){                  #if we have gene we must include
      pop1 = t(replicate(popsize, gene2include.id))
      if (length(gene2include.id)<k){
        if (length(base::setdiff(indices, gene2include.id))==1){
          initpop.candidate = base::setdiff(indices, gene2include.id)
        }else{
          initpop.candidate = sample(base::setdiff(indices, gene2include.id), (k-length(gene2include.id)))
        }
        pop2 = t(replicate(popsize, initpop.candidate))
        pop = cbind(pop1, pop2)
      }
      if (length(gene2include.id)==k){
        pop = pop1
      }
    }else{
      pop = t(replicate(popsize, sample(indices, k)))
    }
  }else{
    pop = initpop
  }

  ###################################
  #calculate diversity of population#
  ###################################
  diversity[1] = popDiv(pop)                                         #calculate diversity of the initial population

  ###########################
  #generate fitness function#
  ###########################
  ###generate fitness function for single-thread case###
  if (!useparallel) {
    getfitness = function(P) {
      f = vector(mode = "numeric", length = popsize)        #fitness value for each chromosome
      cfm = vector(mode = "list", length = popsize)         #confusion matrix for each chromosome
      norm_cfm = vector(mode = "list", length = popsize)         #normalized confusion matrix for each chromosome
      stats_by_class = vector(mode = "list", length = popsize)   #classification statistics by class for each chromosome
      pred_prob = vector(mode = "list", length = popsize)
      AUC_by_class = vector(mode = "list", length = popsize)
      variable_importance = vector(mode = "list", length = popsize)      #variable importance for each chromosome
      for (i in 1:popsize){
        fitness_result = OF(P[i, ], ...)
        f[i] = fitness_result$fitness_value
        cfm[[i]] = fitness_result$confusionMatrix                        #confusion matrix of each chromosome (each one is the average confusion matrix over all cross validations)
        norm_cfm[[i]] = fitness_result$norm.confusionMatrix
        stats_by_class[[i]] = fitness_result$stats_by_class
        pred_prob[[i]] = fitness_result$pred_prob
        AUC_by_class[[i]] = fitness_result$AUC_by_class
        variable_importance[[i]] = fitness_result$var_imp
      }
      best.pos=which(f==min(f))[1]
      ave_cfm=Reduce("+", cfm) / length(cfm)                 #average over all confusion matrices from all chromosomes
      best_cfm=cfm[[best.pos]]                    #small fitness value means high accuracy
      norm_ave_cfm=Reduce("+", norm_cfm) / length(norm_cfm)
      best_norm_ave_cfm=norm_cfm[[best.pos]]
      ave_stats_by_class=Reduce("+", stats_by_class) / length(stats_by_class)
      best_stats_by_class=stats_by_class[[best.pos]]
      best_predprob=pred_prob[[best.pos]]
      best_AUCbyclass=AUC_by_class[[best.pos]]
      best_featureimp=variable_importance[[best.pos]]
      return(list(f=f, ave_cfm=ave_cfm, best_cfm=best_cfm, norm_ave_cfm=norm_ave_cfm, best_norm_ave_cfm=best_norm_ave_cfm,
                  ave_stats_by_class=ave_stats_by_class, best_stats_by_class=best_stats_by_class,
                  best_predprob=best_predprob, best_AUCbyclass=best_AUCbyclass, best_featureimp=best_featureimp,
                  best.pos=best.pos))
    }
  }
  ###generate fitness function for multi-thread case###
  if (useparallel) {
    list(...)
    #multi-thread & no shared memory
    fcn = function(v) OF(v, ...)                                #fcn here is the same to the fitness function (v is a subset of P)
    submit_fcn = function(index, pop.matrix){
      fcn(pop.matrix[index,])
    }
    getfitness = function(P, cluster) {
      parallel::mclapply(1:dim(P)[1], submit_fcn, pop.matrix = P, mc.cores = cluster)
      #         mc.preschedule = FALSE, affinity.list = rep(1:cluster,length.out=dim(P)[1]))
    }
  }

  if (save.intermediate){
    saveRDS(pop, "initial_pop.rds")
  }

  #####################################
  #Actual calculation of fitness value#
  #####################################
  ####for single-thread case or multi-thread, calculate fitness on the initial population###
  if (useparallel){
    #multi-thread
    pop.fitness = getfitness(pop, cluster = cluster)                  #for multi-thread without shared memory, the getfitness function is the same with OF

    fitness.old = sapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$fitness_value)                   #get fitness value of each chromosome
    best.pos=which(fitness.old==min(fitness.old))[1]

    ave_cr_confusionMatrix = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$confusionMatrix)      #confusion matrix for each chromosome
    ave_confusionMatrix[[1]] = Reduce("+", ave_cr_confusionMatrix)/length(ave_cr_confusionMatrix)             #average over all confusion matrices
    best_confusionMatrix[[1]] = ave_cr_confusionMatrix[[best.pos]]             #confusion matrix with best fitness (small fitness means high accuracy)

    norm_ave_cr_confusionMatrix = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$norm.confusionMatrix)      #normalized confusion matrix for each chromosome
    norm_ave_confusionMatrix[[1]] = Reduce("+", norm_ave_cr_confusionMatrix)/length(norm_ave_cr_confusionMatrix)             #average over all normalized confusion matrices
    norm_best_confusionMatrix[[1]] = norm_ave_cr_confusionMatrix[[best.pos]]

    stats_cr_byclass = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$stats_by_class)      #classification stats by class for each chromosome
    stats_byclass[[1]] = Reduce("+", stats_cr_byclass)/length(stats_cr_byclass)             #average over all classification stats by class
    best_stats_byclass[[1]] = stats_cr_byclass[[best.pos]]

    pred_prob_population = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$pred_prob)         #prediction probability for each chromosome
    best_pred_prob[[1]] = pred_prob_population[[best.pos]]                                               #prediction probability of the best chromosome
    AUC_by_class_population = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$AUC_by_class)      #AUC by class for each chromosome
    best_AUC_by_class[[1]] = AUC_by_class_population[[best.pos]]                                            #AUC by class of the best chromosome
    feature_importance_population  = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$var_imp)   #feature importnace for each chromosome
    best_feature_imp[[1]]=feature_importance_population[[best.pos]]                                        #feature importance of the best chromosome
  }
  if (!useparallel){
    #single thread
    pop.fitness = getfitness(pop)

    fitness.old = pop.fitness$f

    ave_confusionMatrix[[1]] = pop.fitness$ave_cfm
    best_confusionMatrix[[1]] = pop.fitness$best_cfm
    norm_ave_confusionMatrix[[1]]=pop.fitness$norm_ave_cfm
    norm_best_confusionMatrix[[1]]=pop.fitness$best_norm_ave_cfm
    stats_byclass[[1]]=pop.fitness$ave_stats_by_class
    best_stats_byclass[[1]]=pop.fitness$best_stats_by_class
    best_pred_prob[[1]]=pop.fitness$best_predprob
    best_AUC_by_class[[1]]=pop.fitness$best_AUCbyclass
    best_feature_imp[[1]]=pop.fitness$best_featureimp
    best.pos=pop.fitness$best.pos
  }


  #old$best[1, ] = sort(pop[rank(fitness.old, ties.method = "random") == 1, ])      #store the best solution (this needs to change! We need to choose the one like this:which(fitness.old==min(fitness.old))[1]. Otherwise the best solution may not corresponds to other information when there are ties)
  old$best[1, ] = sort(pop[best.pos, ])                                             #store the best solution
  old$obj[1] = min(fitness.old)                                                     #store the best fitness value (small fitness -> better result)
  old$avg[1] = mean(fitness.old)                                                    #average fitness value of the current population

  ###output information of the current generation###
  if (verbose > 0) {
    base::cat("Initial population best OF value = ", old$obj[1],
        ". Population diversity = ",  diversity[1], "\n")
  }

  #####################
  #for each generation#
  #####################
  for (gen in 1:ngen) {
    ###generate offspring based on the performance on the initial population###
    ###Selection###
    Tourneys[, ] = t(replicate(popsize, sample(1:popsize,
                                                2*tourneysize)))
    Tourneys1 = Tourneys[, 1:tourneysize]
    Tourneys2 = Tourneys[, (tourneysize + 1): (2*tourneysize)]

    # Tourneys1[, ] = t(replicate(popsize, sample(1:popsize,                 #for each solution/chromosome in popsize, randomly pick tourneysize of solutions as potential parent1 (for example, popsize = 3 and tourneysize = 2, then each row of Tourneys1 here is a random combination of 2 values from 1,2,3 and Tourneys1 has popsize number of rows with different values)
    #                                              tourneysize)))
    # Tourneys2[, ] = t(replicate(popsize, sample(1:popsize,                 #for each solution/chromosome in popsize, randomly pick tourneysize of solutions as potential parent1
    #                                              tourneysize)))

    pickfun = function(v) sample(v, 1, prob = base::rank(-fitness.old[v]))       #solutions/chromosomes with low fitness value -> higher probability of being picked
    Parents1 = apply(Tourneys1, 1, pickfun)                                #for each row in Tourneys1 (tourneysize solutions/chromosomes), pick the solution/chromosome with best fitness from tourneysize of solutions as parent 1
    Parents2 = apply(Tourneys2, 1, pickfun)                                #for each row in Tourneys2 (tourneysize solutions/chromosomes, pick the solution/chromosome with best fitness from tourneysize of solutions as parent 2

    ###Crossover###
    for (i in 1:popsize) {                                                  #generate offspring through crossover
      p1 = pop[Parents1[i], ]
      p2 = pop[Parents2[i], ]
      combo = unique(c(p1, p2))             #for each solution/chromosome, we have a pair of parents. Here we take all potential genes in the two parents. If we have genes must include, these genes will show up in combo for sure because every solution has these genes

      ###Crossover by uniform crossover###
      if (is.null(gene.weight)){
        loc = sample(0:1, k, replace = T)      #for each location, randomly generate 0 and 1 values.
        offspring[i, ] = p1*(1-loc)+p2*loc                    #0 means selecting from p1 and 1 means selecting from p2
        duplicated.loc = base::duplicated(offspring[i, ])         #there may be duplications
        if (sum(duplicated.loc)>0){
          crossover.candidate = base::setdiff(combo, offspring[i, ])
          if (length(crossover.candidate)==1){
            offspring[i, duplicated.loc] = crossover.candidate
          }else{
            offspring[i, duplicated.loc] = sample(crossover.candidate, sum(duplicated.loc))    #we replace those duplications with genes in the two parents but not in the current offspring
          }
        }
      }else{
        #get gene weight at each position for p1 and p2
        weight.matrix = matrix(NA, nrow = 2, ncol = k)
        weight.matrix[1, ] = gene.weight[p1, "weight"]     #weight for genes in p1
        weight.matrix[2, ] = gene.weight[p2, "weight"]     #weight for genes in p2

        #for each position, select 0 (representing gene from p1) or 1 (representing gene from p2) based on the two genes' weight
        loc = apply(weight.matrix, 2, function(x) sample(0:1, 1, prob = x))

        offspring[i, ] = p1*(1-loc)+p2*loc                    #0 means selecting from p1 and 1 means selecting from p2

        duplicated.loc = base::duplicated(offspring[i, ])         #there may be duplications
        if (sum(duplicated.loc)>0){                          #if there are duplications
          replacement.pool = base::setdiff(combo, offspring[i, ])      #we replace those duplications with genes in the two parents but not in the current offspring
          if (length(replacement.pool)==1){
            offspring[i, duplicated.loc] = replacement.pool
          }else{
            offspring[i, duplicated.loc] = sample(replacement.pool, sum(duplicated.loc), prob = gene.weight[replacement.pool, "weight"])    #genes with higher weight will be more likely to be selected
          }
        }
      }

      ###Crossover by random sampling from unique values###
      # if (is.null(gene.weight)){
      #   offspring[i, ] = sample(combo, k)                                     #then we carried out crossover (reproduction) by combining the unique elements of both parents and keeping k of them, chosen at random.
      # }else{
      #   offspring[i, ] = sample(combo, k, prob = gene.weight[combo, "weight"]) #then we carried out crossover (reproduction) by combining the unique elements of both parents and keeping k of them, chosen given each genes weight.
      # }
      # if (length(unique(offspring[i, ]))<50){
      #   saveRDS(list(p1=p1, p2=p2, loc=loc, offspring.panel = offspring[i, ]), "problem_in_crossover.rds")
      # }
    }

    ###Mutation###
    chosen = matrix(as.logical(rbinom(popsize * k, 1, mutprob)),           #randonly pick genes in each solution to mutate based on the mutation probability
                     nrow = popsize)
    nchosen = apply(chosen, 1, sum)                                        #total number of genes to mutate for each solution

    for (i in 1:popsize) {                                                  #for each solution
      if (nchosen[i] > 0) {                                                 #if we have at least one gene to mutate
        candidates = data.table::indices[-offspring[i, ]]                               #genes not in the current solution/chromosome
        if (is.null(gene.weight)){
          if (length(candidates)==1){
            toadd = candidates
          }else{
            toadd = sample(candidates, nchosen[i])                            #randomly select genes from genes not in the current solution
          }
        }else{
          if (length(candidates)==1){
            toadd = candidates
          }else{
            toadd = sample(candidates, nchosen[i], prob = gene.weight[candidates, "weight"])               #select genes from genes not in the current solution given their weight
          }
        }
        offspring[i, chosen[i, ]] = toadd                                   #add the selected genes to the new solution to replace old ones
      }

      #add genes must include (if we have them) back to the offspring population
      if (!is.null(gene2include.id)){
        #offspring[i,] = sample(unique(c(offspring[i,], gene2include.id)), k, replace = FALSE)
        already.in = base::intersect(offspring[i,], gene2include.id)          #must include genes that are already in the current solution
        if (length(already.in)<length(gene2include.id)){                #if some/all of the must include genes are not in the solution
          if (length(already.in)>0){                                    #if a subset of must include genes are already in the solution
            pos2change.pool = which(!(offspring[i,] %in% already.in))
            if (length(pos2change.pool)==1){                            #if there is only one position to change, sample(x, 1) will not give you x, instead, it will give you a random number between 1 and x. To avoid this, we separate length(pos2change.pool)=1 and >1
              pos2change = pos2change.pool
            }else{
              if (is.null(gene.weight)){
                pos2change = sample(pos2change.pool, length(gene2include.id)-length(already.in))     #we only update the rest (if the panel size is 30. We have 10 genes we must include and 3 of them are already included. Then this step will select 7 positions from the 27 positions in the panel and change them to the 7 genes we must include but haven't included yet)
              }else{
                pos2change = sample(pos2change.pool, length(gene2include.id)-length(already.in), prob = 1-gene.weight[offspring[i,pos2change.pool], "weight"])      #genes with higher weight will be less likely to be replaced by must include values (if we have more positions we can change then the genes we need to include)
              }
            }
            offspring[i, pos2change] = base::setdiff(gene2include.id, already.in)
          }
          if (length(already.in)==0){                                    #if none of the must include genes are in the solution
            if (is.null(gene.weight)){
              offspring[i, sample(1:k, length(gene2include.id))]=gene2include.id
            }else{
              pos2change = sample(1:k, length(gene2include.id), prob = 1-gene.weight[offspring[i, 1:k], "weight"])             #genes with higher weight will be less likely to be replaced by must include values (if we have more positions we can change then the genes we need to include)
              offspring[i, pos2change]=gene2include.id
            }
          }
        }
      }
    }

    ###calculate fitness for offspring###
    if (useparallel){
      #multi-thread
      pop.fitness = getfitness(offspring, cluster = cluster)                  #for multi-thread without shared memory, the getfitness function is the same with OF

      fitness.new = sapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$fitness_value)
      best.pos=which(fitness.new==min(fitness.new))[1]

      ave_cr_confusionMatrix = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$confusionMatrix)
      ave_confusionMatrix[[gen + 1]] = Reduce("+", ave_cr_confusionMatrix)/length(ave_cr_confusionMatrix)
      best_confusionMatrix[[gen + 1]] = ave_cr_confusionMatrix[[best.pos]]

      norm_ave_cr_confusionMatrix = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$norm.confusionMatrix)      #normalized confusion matrix for each chromosome
      norm_ave_confusionMatrix[[gen + 1]] = Reduce("+", norm_ave_cr_confusionMatrix)/length(norm_ave_cr_confusionMatrix)             #average over all normalized confusion matrices
      norm_best_confusionMatrix[[gen + 1]] = norm_ave_cr_confusionMatrix[[best.pos]]

      stats_cr_byclass = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$stats_by_class)      #classification stats by class for each chromosome
      stats_byclass[[gen + 1]] = Reduce("+", stats_cr_byclass)/length(stats_cr_byclass)             #average over all classification stats by class
      best_stats_byclass[[gen + 1]] = stats_cr_byclass[[best.pos]]

      pred_prob_population = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$pred_prob)         #prediction probability for each chromosome
      best_pred_prob[[gen + 1]] = pred_prob_population[[best.pos]]
      AUC_by_class_population = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$AUC_by_class)      #AUC by class for each chromosome
      best_AUC_by_class[[gen + 1]] = AUC_by_class_population[[best.pos]]
      feature_importance_population  = lapply(1:length(pop.fitness), function(x) pop.fitness[[x]]$var_imp)   #feature importnace for each chromosome
      best_feature_imp[[gen + 1]]=feature_importance_population[[best.pos]]                                        #feature importance of the best chromosome
    }
    if (!useparallel){
      #single-thread
      pop.fitness = getfitness(offspring)

      fitness.new = pop.fitness$f

      ave_confusionMatrix[[gen + 1]] = pop.fitness$ave_cfm
      best_confusionMatrix[[gen + 1]] = pop.fitness$best_cfm
      norm_ave_confusionMatrix[[gen + 1]]=pop.fitness$norm_ave_cfm
      norm_best_confusionMatrix[[gen + 1]]=pop.fitness$best_norm_ave_cfm
      stats_byclass[[gen + 1]]=pop.fitness$ave_stats_by_class
      best_stats_byclass[[gen + 1]]=pop.fitness$best_stats_by_class
      best_pred_prob[[gen + 1]]=pop.fitness$best_predprob
      best_AUC_by_class[[gen + 1]]=pop.fitness$best_AUCbyclass
      best_feature_imp[[gen + 1]]=pop.fitness$best_featureimp
      best.pos=pop.fitness$best.pos
    }


    ###update population###
    #when keepbest!=0, the best solution of each generation will be the best from both the previous generation and the current generation. Therefore, the fitness value is guaranteed to improve
    if (keepbest == 0) {                                                  #update population when keepbest==0
      pop = offspring                                                    #here we simply replace the old population with the new population
      fitness.old = fitness.new                                          #fitness value is also replaced by the new fitness values
    }else {                                                                            #update population when keepbest!=0 (instead of replacing everything with the new population, here we keep the best of the old population and replace others with the new population)
      old.keep = base::rank(fitness.old, ties.method = "first") <= keepbest                #we keep the top keepbest solutions in the old population
      new.keep = base::rank(fitness.new, ties.method = "first") <= popsize - keepbest      #for the rest of the solutions, we replace them with the top (popsize-keepbest) solutions from the offspring population
      pop[elitespots, ] = pop[old.keep, ]                                 #update population (add the top solutions in the old population back to the population)
      pop[newspots, ] = offspring[new.keep, ]                             #add the top solutions from the offspring population to the rest of the population
      fitness.old[elitespots] = fitness.old[old.keep]                     #update fitness values of solutions
      fitness.old[newspots] = fitness.new[new.keep]
    }

    ###save the result for the current generation###
    if (keepbest == 0){                                 #if we completely replace the old population with the offspring population (now pop=offspring)
      old$best[gen + 1, ] = sort(pop[best.pos, ])      #then the best.pos for the offspring population will be the best solution
    }else{
      #if the new population is a mixture of the previous generation and current generation (offspring population)
      #then the best solution is either the best solution in the previous generation or the best solution in the current generation (i.e., offspring[best.pos, ])
      if (old$obj[gen] < min(fitness.new)){                     #if the best solution in the previous round is better than the one in the current round
        old$best[gen + 1, ] = old$best[gen, ]                   #then the best solution in the mixed population will be the best solution of the previous round

        #we need to change all the other result to the result from the previous generation also
        ave_confusionMatrix[[gen + 1]] = ave_confusionMatrix[[gen]]
        best_confusionMatrix[[gen + 1]] = best_confusionMatrix[[gen]]
        norm_ave_confusionMatrix[[gen + 1]] = norm_ave_confusionMatrix[[gen]]
        norm_best_confusionMatrix[[gen + 1]] = norm_best_confusionMatrix[[gen]]
        stats_byclass[[gen + 1]] = stats_byclass[[gen]]
        best_stats_byclass[[gen + 1]] = best_stats_byclass[[gen]]
        best_pred_prob[[gen + 1]] = best_pred_prob[[gen]]
        best_AUC_by_class[[gen + 1]] = best_AUC_by_class[[gen]]
        best_feature_imp[[gen + 1]] = best_feature_imp[[gen]]
      }
      if (old$obj[gen] >= min(fitness.new)){                    #if the best solution in the generation round is better
        old$best[gen + 1, ] = sort(offspring[best.pos, ])      #then the best solution in the mixed population will be the best.pos solution in the offspring population (of note, pop != offspring in this case so we need to use offspring[best.pos,] instead of pop[best.pos,])
        #we don't need to modify the rest of the result because they are already based on the offspring population
      }
    }

    if (verbose > 0){
      saveRDS(pop, paste0("pop",gen,".rds"))
    }

    old$obj[gen + 1] = min(fitness.old)                  #this doesn't change because even if keepbest!=0, the best solution will be the minimum of the best solution in two rounds and its fitness value is the minimum so this information matches with old$best[gen + 1, ]
    old$avg[gen + 1] = mean(fitness.old)

    ###################################
    #calculate diversity of population#
    ###################################
    diversity[gen + 1] = popDiv(pop)                                    #calculate diversity of the new population

    ###output information of the current generation###
    if (verbose > 0) {
      cat("Finished iteration ", gen,
          ". Best OF value = ", old$obj[gen + 1],
          ". Population diversity = ",  diversity[gen + 1], "\n")
    }

    #####################################################################################
    #Early termination if best fitness doesn't change after a given number of iterations#
    #####################################################################################
    if (gen >= earlyterm){
      fitness.short.history = old$obj[(gen-earlyterm + 1): (gen + 1)]      #get the fitness value in the past earlyterm iterations (past 3 iteration means fitness value of 4 iterations (including the current one))
      improvement = abs(diff(fitness.short.history))                       #get the fitness difference between each iteration and its previous iteration
      if (max(improvement) <= converge.cutoff){                                #if the maximum improvement is smaller than the cutoff, we stop the optimization
        break
      }
    }
  }

  ############################
  #output for all generations#
  ############################
  out = list()

  old$best = old$best[1:(gen+1),]       #if early terminate the program, then we only keep the result that is already finished
  old$obj = old$obj[1:(gen+1)]
  old$avg = old$avg[1:(gen+1)]

  out$old = old                     #this contains result for each generation

  ord = order(fitness.old)
  out$pop = matrix(0, nrow = popsize, ncol = k)
  for (i in 1:popsize) {
    # out$pop[i, ] = sort(pop[ord[i], ])
    out$pop[i, ] = pop[ord[i], ]                                   #This is the last population, with solutions in increasing order (from best to worst)
  }
  out$obj = fitness.old[ord]                                      #This is the fitness value for all solutions in the last population in increasing order (from best to worst)

  alltimebest = which(old$obj == min(old$obj, na.rm = TRUE))      #among all the generations, which one have the best fitness (we may have situations where the best solution of multiple generations have the same fitness value)
  #alltimebest = which(old$obj == min(old$obj[1:(gen+1)], na.rm = TRUE))
  alltimebest = tail(alltimebest, 1)                              #we choose the last generate as the best when there are ties
  out$bestgeneration = alltimebest
  out$bestsol = out$old$best[alltimebest, ]                       #the best solution from the best generation (usually the last generation)
  out$bestobj = out$old$obj[alltimebest]                          #fitness value of the best solution from the best generation

  #add diversity of population to output#
  out$diversity = diversity[1:(gen+1)]

  #add confusion matrix of population to output#
  out$ave_confusionMatrix = ave_confusionMatrix[1:(gen+1)]
  out$best_confusionMatrix = best_confusionMatrix[1:(gen+1)]
  out$norm_ave_confusionMatrix = norm_ave_confusionMatrix[1:(gen+1)]
  out$norm_best_confusionMatrix = norm_best_confusionMatrix[1:(gen+1)]

  #add stats by class to output#
  out$stats_byclass=stats_byclass[1:(gen+1)]
  out$best_stats_byclass=best_stats_byclass[1:(gen+1)]

  #add prediction probability and by class AUC#
  out$best_pred_prob=best_pred_prob[1:(gen+1)]
  out$best_AUC_by_class=best_AUC_by_class[1:(gen+1)]

  #add feature importance#
  out$best_feature_importance=best_feature_imp[1:(gen+1)]

  class(out) = "gpsFISH"
  out
}

