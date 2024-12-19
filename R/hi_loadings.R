#' Extract genes with highest loadings
#'
#' @param pcaobj A `prcomp` object
#' @param whichpc An integer number, corresponding to the principal component of
#' interest
#' @param topN Integer, number of genes with top and bottom loadings
#' @param exprTable A `matrix` object, e.g. the counts of a [DESeq2::DESeqDataSet()].
#' If not NULL, returns the counts matrix for the selected genes
#' @param annotation A `data.frame` object, with row.names as gene identifiers (e.g. ENSEMBL ids)
#' and a column, `gene_name`, containing e.g. HGNC-based gene symbols
#' @param title The title of the plot
#'
#' @return A ggplot2 object, or a `matrix`, if `exprTable` is not null
#'
#' @examples
#' dds <- makeExampleDESeqDataSet_multifac(betaSD = 3, betaSD_tissue = 1)
#' rlt <- DESeq2::rlogTransformation(dds)
#' pcaobj <- prcomp(t(SummarizedExperiment::assay(rlt)))
#' hi_loadings(pcaobj, topN = 20)
#' hi_loadings(pcaobj, topN = 10, exprTable = dds)
#' hi_loadings(pcaobj, topN = 10, exprTable = counts(dds))
#'
#' @export
hi_loadings <- function(pcaobj, whichpc = 1, topN = 10, exprTable = NULL,
                        annotation = NULL, title="Top/bottom loadings") {
  if (whichpc < 0)
    stop("Use a positive integer value for the principal component to select")
  if (whichpc > nrow(pcaobj$x))
    stop("You can not explore a principal component that is not in the data")

  geneloadings_sorted <- sort(pcaobj$rotation[, whichpc])
  geneloadings_extreme <- c(tail(geneloadings_sorted, topN), head(geneloadings_sorted, topN))

  if (!is.null(exprTable)) {
    tab <- exprTable[names(geneloadings_extreme), ]
    if (!is.null(annotation))
      rownames(tab) <- annotation$gene_name[match(rownames(tab), rownames(annotation))]
    return(tab)
  }

  if (!is.null(annotation))
    names(geneloadings_extreme) <- annotation$gene_name[match(names(geneloadings_extreme), rownames(annotation))]

  barplot(geneloadings_extreme, las = 2, col = c(rep("steelblue", topN), rep("coral", topN)),
          main = paste0(title, "PC", whichpc))
  # mydf <- data.frame(loadings=geneloadings_extreme,
  #                    geneID=names(geneloadings_extreme),
  #                    mycol = c(rep("steelblue",topN),rep("coral",topN)))
  # mydf$geneID <- factor(mydf$geneID, levels = mydf$geneID)
  # p <- ggplot(mydf,aes_string(x="geneID",y="loadings")) + geom_col(aes_string(fill = "mycol")) + theme_bw() +
  #   theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + guides(fill = FALSE) +
  #   ggtitle(paste0(title, " - PC", whichpc))
  # p
}
