
[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/gpsFISH.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/gpsFISH)

# gpsFISH
- [Introduction](#introduction)
- [Installation](#installation)
- [Tutorials](#tutorials)
  * [Platform effect estimation](#platform-effect-estimation)
  * [Gene panel selection](#gene-panel-selection)
- [Citation](#citation)


## Introduction
Methods for spatial transcriptomic analysis based on highly multiplexed single-molecule fluorescence in situ hybridization or sequencing (hm-smFISH) hold particular promise in analysis of complex tissues. Most such methods, however, measure only a limited panel of transcripts, which need to be selected in advance to inform on the cell types or processes being studied. 

Here we describe gpsFISH, a computational method to perform such selection in a way that optimizes detection of previously annotated cell types or cell type hierarchies. We show systematic difference of transcript detection rates between different technologies, which distorts the resulting transcriptional profile estimates. Such distortions reduce the resulting ability to distinguish cell types in the hm-smFISH data relative to what would be expected from scRNA-seq data. gpsFISH was designed to model and adjust for the platform effects, yielding more informative gene panels and better cell type classification. In addition, gpsFISH can incorporate the hierarchical structure of cell types to better account for complex cell type relationships. Finally, we observed a high redundancy of gene panels, and extended gpsFISH to incorporate custom gene preferences. Overall, gpsFISH outperforms other gene selection methods, improving support for hierarchical cell classification as well as flexible options to accommodate diverse design requirements.

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
* [Markdown version](doc/simulation_training.md)

### Gene panel selection:
* [HTML version](https://htmlpreview.github.io/?https://github.com/kharchenkolab/gpsFISH/blob/main/doc/gene_panel_selection.html)
* [Markdown version](doc/gene_panel_selection.md)


## Citation

If you find gpsFISH useful for your publication, please cite:
