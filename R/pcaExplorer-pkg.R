#' pcaExplorer: analyzing time-lapse microscopy imaging, from detection to tracking
#'
#' pcaExplorer provides functionality for interactive visualization of RNA-seq datasets
#' based on Principal Components Analysis. The methods provided allow for quick information
#' extraction and effective data exploration. A Shiny application encapsulates the whole analysis.
#'
#' pcaExplorer provides functionality for interactive visualization of RNA-seq datasets
#' based on Principal Components Analysis. The methods provided allow for quick information
#' extraction and effective data exploration. A Shiny application encapsulates the whole analysis.
#'
#' @import DESeq2
#' @import SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom genefilter rowVars
#' @importFrom heatmaply heatmaply
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom scales brewer_pal hue_pal
#' @importFrom NMF aheatmap
#' @import plyr
#' @importFrom limma goana topGO
#' @importFrom AnnotationDbi select Term mapIds
#' @importMethodsFrom GOstats hyperGTest summary
#' @import GO.db
#' @import shiny
#' @import shinydashboard
#' @importFrom shinyBS bsTooltip bsCollapse bsCollapsePanel
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom DT datatable
#' @importFrom shinyAce aceAutocomplete aceEditor getAceModes getAceThemes
#'   updateAceEditor
#' @import threejs
#' @import biomaRt
#' @importFrom pheatmap pheatmap
#' @importFrom base64enc dataURI
#' @importFrom tidyr gather
#' @import knitr
#' @import rmarkdown
#' @importFrom grDevices dev.off pdf rainbow colorRamp rgb
#' @import methods
#'
#' @author
#' Federico Marini \email{marinif@@uni-mainz.de}, 2016
#'
#' Maintainer: Federico Marini \email{marinif@@uni-mainz.de}
#' @name pcaExplorer-pkg
#' @docType package
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription("pcaExplorer", fields = "Version")
  msg <- paste0("Welcome to pcaExplorer v", pkgVersion, "\n\n")
  citation <- paste0("If you use pcaExplorer in your work, please cite:\n\n",
                     "pcaExplorer: an R/Bioconductor package for interacting with RNA-seq principal components\n",
                     "Federico Marini, Harald Binder\n",
                     "BMC Bioinformatics, 2019 - https://doi.org/10.1186/s12859-019-2879-1\n")
  packageStartupMessage(paste0(msg, citation))
}
