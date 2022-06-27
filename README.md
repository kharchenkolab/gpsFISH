# gpsFISH
- [Introduction](#introduction)
- [Installation](#installation)
- [Tutorials](#tutorials)
  * [Platform effect estimation](#platform-effect-estimation)
  * [Gene panel selection](#gene-panel-selection)
- [Citation](#citation)


## Introduction
Methods for spatial transcriptomic analysis based on highly multiplexed single-molecule fluorescence in situ hybridization (hm-smFISH) hold particular promise in analysis of complex tissues. Most such methods, however, measure only a panel of genes, which need to be selected in advance to recognize cell types or processes relevant to the tissue being studied. 

Here, we develop GPS-FISH, a computational method that selects a panel of genes to optimize detection of previously annotated cell types or cell type hierarchies, based on scRNA-seq data. Systematic difference of transcript detection rates is observed between different technologies, which distort the resulting transcriptional profile estimates. Such distortions lead to reduced ability to distinguish cell types in the resulting hm-smFISH data. GPS-FISH was designed to model and adjust for the platform effects, yielding more informative gene panels and better cell type classification. In addition, GPS-FISH can incorporate the hierarchical structure of cell types to account for complex cell type relationships. Finally, we observed a high redundancy of gene panels and take advantage of this redundancy to provide flexibility for various custom preferences. Overall, GPS-FISH outperforms other gene selection methods by taking platform effects across technologies into consideration and provides flexible options for cell type hierarchy and complex custom preferences. 

<img src="inst/workflow.jpg" align="center" height="600">

## Installation
To use the latest version of gpsFISH from GitHub, install with the following:

``` r
devtools::install_github("kharchenkolab/gpsFISH")
```

## Tutorials

Please see the following tutorials for detailed examples of how to use gpsFISH: 

### Platform effect estimation:
* [HTML version](https://htmlpreview.github.io/?https://github.com/kharchenkolab/gpsFISH/blob/main/doc/simulation_training.html)
* [Markdown version](doc/simulation_training.Rmd)

### Gene panel selection:
* [HTML version](https://htmlpreview.github.io/?https://github.com/kharchenkolab/gpsFISH/blob/main/doc/gene_panel_selection.html)
* [Markdown version](doc/gene_panel_selection.Rmd)


## Citation

If you find gpsFISH useful for your publication, please cite:
