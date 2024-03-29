% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{initialize_population_random}
\alias{initialize_population_random}
\title{Initialize population randomly from all candidate genes}
\usage{
initialize_population_random(pop.size, panel.size, gene.list, gene2include.id)
}
\arguments{
\item{pop.size}{Size of the population, i.e., the number of gene panels to initialize}

\item{panel.size}{Size of a gene panel, i.e., the number of genes in a gene panel}

\item{gene.list}{A character vector of all candidate genes used to initialize the population}

\item{gene2include.id}{A numeric vector specifying the location of genes in \code{gene.list} that must be included in each panel of the population. Default is NULL}
}
\value{
A matrix with each row representing one gene panel and each column representing one gene in a gene panel.
The genes are encoded by their location in \code{gene.list}
}
\description{
Initialize population randomly from all candidate genes
}
\examples{
data(sc_count)

gene2include.symbol = sample(rownames(sc_count), 20)
gene2include.id=which(rownames(sc_count) \%in\% gene2include.symbol)

initpop.random=initialize_population_random(pop.size = 100,
                                            panel.size = 200,
                                            gene.list = rownames(sc_count),
                                            gene2include.id = gene2include.id)
head(initpop.random)
}
