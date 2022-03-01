# pcaExplorer - Interactive exploration of Principal Components of Samples and Genes in RNA-seq data

<a href="https://doi.org/10.1186/s12859-019-2879-1"><img src="https://img.shields.io/badge/doi-pcaExplorer-blue.svg"><a>

## Software status

[![R build status](https://github.com/federicomarini/pcaExplorer/workflows/R-CMD-check/badge.svg)](https://github.com/federicomarini/pcaExplorer/actions)

| Platforms |  OS  | R CMD check |
|:----------------:|:----------------:|:----------------:|
| Bioc ([_devel_](http://bioconductor.org/packages/devel/bioc/html/pcaExplorer.html)) | Multiple | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/pcaExplorer.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/pcaExplorer) |
| Bioc ([_release_](http://bioconductor.org/packages/release/bioc/html/pcaExplorer.html)) | Multiple | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/pcaExplorer.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/pcaExplorer) |

[![codecov.io](https://codecov.io/github/federicomarini/pcaExplorer/coverage.svg?branch=master)](https://codecov.io/github/federicomarini/pcaExplorer?branch=master)

`pcaExplorer` is a Bioconductor package containing a Shiny application for
analyzing expression data in different conditions and experimental factors. 

It is a general-purpose interactive companion tool for RNA-seq analysis, which 
guides the user in exploring the Principal Components of the data under inspection.

`pcaExplorer` provides tools and functionality to detect outlier samples, genes
that show particular patterns, and additionally provides a functional interpretation of 
the principal components for further quality assessment and hypothesis generation
on the input data. 

Moreover, a novel visualization approach is presented to simultaneously assess 
the effect of more than one experimental factor on the expression levels.

Thanks to its interactive/reactive design, it is designed to become a practical
companion to any RNA-seq dataset analysis, making exploratory data analysis 
accessible also to the bench biologist, while providing additional insight also
for the experienced data analyst.

## Installation

`pcaExplorer` can be easily installed using `BiocManager::install()`:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("pcaExplorer")
```

or, optionally, 

``` r
BiocManager::install("federicomarini/pcaExplorer")
# or alternatively...
devtools::install_github("federicomarini/pcaExplorer")
```

## Quick start

This command loads the `pcaExplorer` package

``` r
library("pcaExplorer")
```

The `pcaExplorer` app can be launched in different modes:

- `pcaExplorer(dds = dds, dst = dst)`, where `dds` is a `DESeqDataSet` object and `dst` is a `DESeqTransform`
object, which were created during an existing session for the analysis of an RNA-seq
dataset with the `DESeq2` package

- `pcaExplorer(dds = dds)`, where `dds` is a `DESeqDataSet` object. The `dst` object is automatically 
computed upon launch.

- `pcaExplorer(countmatrix = countmatrix, coldata = coldata)`, where `countmatrix` is a count matrix, generated
after assigning reads to features such as genes via tools such as `HTSeq-count` or `featureCounts`, and `coldata`
is a data frame containing the experimental covariates of the experiments, such as condition, tissue, cell line,
run batch and so on.

- `pcaExplorer()`, and then subsequently uploading the count matrix and the covariates data frame through the 
user interface. These files need to be formatted as tab separated files, which is a common format for storing
such count values.

Additional parameters and objects that can be provided to the main `pcaExplorer` function are:

- `pca2go`, which is an object created by the `pca2go` function, which scans the genes with high loadings in 
each principal component and each direction, and looks for functions (such as GO Biological Processes) that 
are enriched above the background. The offline `pca2go` function is based on the routines and algorithms of 
the `topGO` package, but as an alternative, this object can be computed live during the execution of the app
exploiting the `goana` function, provided by the `limma` package. Although this likely provides more general
(and probably less informative) functions, it is a good compromise for obtaining a further data interpretation.

- `annotation`, a data frame object, with `row.names` as gene identifiers (e.g. ENSEMBL ids) identical to the 
row names of the count matrix or `dds` object, and an extra column `gene_name`, containing e.g. HGNC-based 
gene symbols. This can be used for making information extraction easier, as ENSEMBL ids (a usual choice when
assigning reads to features) do not provide an immediate readout for which gene they refer to. This can be
either passed as a parameter when launching the app, or also uploaded as a tab separated text file.

## Contact

For additional details regarding the functions of **pcaExplorer**, please consult the documentation or 
write an email to marinif@uni-mainz.de. 

## Code of Conduct

Please note that the pcaExplorer project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

### Bug reports/Issues/New features

Please use https://github.com/federicomarini/pcaExplorer/issues for reporting bugs, issues or for 
suggesting new features to be implemented. 
