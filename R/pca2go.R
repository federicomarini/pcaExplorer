#' Functional interpretation of the principal components
#'
#' Extracts the genes with the highest loadings for each principal component, and
#' performs functional enrichment analysis on them using routines and algorithms from
#' the \code{topGO} package
#'
#' @param se A \code{\link{DESeqTransform}} object, with data in \code{assay(se)},
#' produced for example by either \code{\link{rlog}} or
#' \code{\link{varianceStabilizingTransformation}}
#' @param pca_ngenes Number of genes to use for the PCA
#' @param annotation A \code{data.frame} object, with row.names as gene identifiers (e.g. ENSEMBL ids)
#' and a column, \code{gene_name}, containing e.g. HGNC-based gene symbols
#' @param inputType Input format type of the gene identifiers. Will be used by the routines of \code{topGO}
#' @param organism Character abbreviation for the species, using \code{org.XX.eg.db} for annotation
#' @param ensToGeneSymbol Logical, whether to expect ENSEMBL gene identifiers, to convert to gene symbols
#' with the \code{annotation} provided
#' @param loadings_ngenes Number of genes to extract the loadings (in each direction)
#' @param background_genes Which genes to consider as background.
#' @param scale Logical, defaults to FALSE, scale values for the PCA
#' @param return_ranked_gene_loadings Logical, defaults to FALSE. If TRUE, simply returns
#' a list containing the top ranked genes with hi loadings in each PC and in each direction
#' @param annopkg String containing the name of the organism annotation package. Can be used to 
#' override the \code{organism} parameter, e.g. in case of alternative identifiers used
#' in the annotation package (Arabidopsis with TAIR)
#' @param ... Further parameters to be passed to the topGO routine
#'
#' @return A nested list object containing for each principal component the terms enriched
#' in each direction. This object is to be thought in combination with the displaying feature
#' of the main \code{\link{pcaExplorer}} function
#'
#' @examples
#' library(airway)
#' library(DESeq2)
#' data(airway)
#' airway
#' dds_airway <- DESeqDataSet(airway, design= ~ cell + dex)
#' \dontrun{
#' rld_airway <- rlogTransformation(dds_airway)
#' # constructing the annotation object
#' anno_df <- data.frame(gene_id = rownames(dds_airway),
#'                       stringsAsFactors = FALSE)
#' library("AnnotationDbi")
#' library("org.Hs.eg.db")
#' anno_df$gene_name <- mapIds(org.Hs.eg.db,
#'                             keys = anno_df$gene_id,
#'                             column = "SYMBOL",
#'                             keytype = "ENSEMBL",
#'                             multiVals = "first")
#' rownames(anno_df) <- anno_df$gene_id
#' bg_ids <- rownames(dds_airway)[rowSums(counts(dds_airway)) > 0]
#' library(topGO)
#' pca2go_airway <- pca2go(rld_airway,
#'                         annotation = anno_df,
#'                         organism = "Hs",
#'                         ensToGeneSymbol = TRUE,
#'                         background_genes = bg_ids)
#' }
#'
#' @export
pca2go <- function(se,
                   pca_ngenes = 10000,
                   annotation = NULL,
                   inputType = "geneSymbol",
                   organism = "Mm",
                   ensToGeneSymbol = FALSE,
                   loadings_ngenes = 500,
                   background_genes = NULL,
                   scale = FALSE,
                   return_ranked_gene_loadings = FALSE,
                   annopkg = NULL,
                   ... # further parameters to be passed to the topgo routine
                   ) {

  # to force e.g. in case other identifiers are used, e.g. in Arabidopsis with TAIR
  if (is.null(annopkg))
    annopkg <- paste0("org.", organism, ".eg.db")
  
  if (!require(annopkg, character.only = TRUE)) {
    stop("The package", annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
  }
  exprsData <- assay(se)

  if (is.null(background_genes)) {
    if (is(se, "DESeqDataSet")) {
      BGids <- rownames(se)[rowSums(counts(se)) > 0]
    } else if (is(se, "DESeqTransform")) {
      BGids <- rownames(se)[rowSums(assay(se)) != 0]
    } else {
      BGids <- rownames(se)
    }
  } else {
    BGids <- background_genes
  }

  exprsData <- assay(se)[order(rowVars(assay(se)), decreasing = TRUE), ]
  exprsData <- exprsData[1:pca_ngenes, ]

  ## convert ensembl to gene symbol ids
  if (ensToGeneSymbol & !(is.null(annotation))) {
    rownames(exprsData) <- annotation$gene_name[match(rownames(exprsData), rownames(annotation))]
    BGids <- annotation$gene_name[match(BGids, rownames(annotation))]
  }

  rv <- rowVars(exprsData)
  dropped <- sum(rv == 0)
  if (dropped > 0)
    message(paste("Dropped", dropped, "genes with 0 variance"))

  exprsData <- exprsData[rv > 0, ]

  message("After subsetting/filtering for invariant genes, working on a ", nrow(exprsData), "x", ncol(exprsData), " expression matrix\n")

  p <- prcomp(t(exprsData), scale = scale, center = TRUE)

  message("Ranking genes by the loadings ...")
  probesPC1pos <- rankedGeneLoadings(p, pc = 1, decreasing = TRUE)[1:loadings_ngenes]
  probesPC1neg <- rankedGeneLoadings(p, pc = 1, decreasing = FALSE)[1:loadings_ngenes]
  probesPC2pos <- rankedGeneLoadings(p, pc = 2, decreasing = TRUE)[1:loadings_ngenes]
  probesPC2neg <- rankedGeneLoadings(p, pc = 2, decreasing = FALSE)[1:loadings_ngenes]
  probesPC3pos <- rankedGeneLoadings(p, pc = 3, decreasing = TRUE)[1:loadings_ngenes]
  probesPC3neg <- rankedGeneLoadings(p, pc = 3, decreasing = FALSE)[1:loadings_ngenes]
  probesPC4pos <- rankedGeneLoadings(p, pc = 4, decreasing = TRUE)[1:loadings_ngenes]
  probesPC4neg <- rankedGeneLoadings(p, pc = 4, decreasing = FALSE)[1:loadings_ngenes]
  
  if (return_ranked_gene_loadings)
    return(
      list(PC1 = list(posLoad = probesPC1pos, negLoad = probesPC1neg),
           PC2 = list(posLoad = probesPC2pos, negLoad = probesPC2neg),
           PC3 = list(posLoad = probesPC3pos, negLoad = probesPC3neg),
           PC4 = list(posLoad = probesPC4pos, negLoad = probesPC4neg)
      )
    )
  
  message("Ranking genes by the loadings ... done!")
  message("Extracting functional categories enriched in the gene subsets ...")
  topGOpc1pos <- topGOtable(probesPC1pos, BGids, annot = annFUN.org, mapping = annopkg, ...)
  topGOpc1neg <- topGOtable(probesPC1neg, BGids, annot = annFUN.org, mapping = annopkg, ...)
  topGOpc2pos <- topGOtable(probesPC2pos, BGids, annot = annFUN.org, mapping = annopkg, ...)
  topGOpc2neg <- topGOtable(probesPC2neg, BGids, annot = annFUN.org, mapping = annopkg, ...)
  topGOpc3pos <- topGOtable(probesPC3pos, BGids, annot = annFUN.org, mapping = annopkg, ...)
  topGOpc3neg <- topGOtable(probesPC3neg, BGids, annot = annFUN.org, mapping = annopkg, ...)
  topGOpc4pos <- topGOtable(probesPC4pos, BGids, annot = annFUN.org, mapping = annopkg, ...)
  topGOpc4neg <- topGOtable(probesPC4neg, BGids, annot = annFUN.org, mapping = annopkg, ...)

  goEnrichs <- list(PC1 = list(posLoad = topGOpc1pos, negLoad = topGOpc1neg),
                    PC2 = list(posLoad = topGOpc2pos, negLoad = topGOpc2neg),
                    PC3 = list(posLoad = topGOpc3pos, negLoad = topGOpc3neg),
                    PC4 = list(posLoad = topGOpc4pos, negLoad = topGOpc4neg)
  )
  message("Extracting functional categories enriched in the gene subsets ... done!")
  attr(goEnrichs, "n_genesforpca") <- pca_ngenes
  return(goEnrichs)
}


rankedGeneLoadings <- function(x, pc = 1, decreasing = TRUE) {
  # works directly on the prcomp object
  return(rownames(x$rotation)[order(x$rotation[, pc], decreasing = decreasing)])
}


#' Extract functional terms enriched in the DE genes, based on topGO
#'
#' A wrapper for extracting functional GO terms enriched in the DE genes, based on
#' the algorithm and the implementation in the topGO package
#' 
#' Allowed values assumed by the \code{topGO_method2} parameter are one of the
#' following: \code{elim}, \code{weight}, \code{weight01}, \code{lea}, 
#' \code{parentchild}. For more details on this, please refer to the original 
#' documentation of the \code{topGO} package itself
#'
#' @param DEgenes A vector of (differentially expressed) genes
#' @param BGgenes A vector of background genes, e.g. all (expressed) genes in the assays
#' @param ontology Which Gene Ontology domain to analyze: \code{BP} (Biological Process), \code{MF} (Molecular Function), or \code{CC} (Cellular Component)
#' @param annot Which function to use for annotating genes to GO terms. Defaults to \code{annFUN.org}
#' @param mapping Which \code{org.XX.eg.db} to use for annotation - select according to the species
#' @param geneID Which format the genes are provided. Defaults to \code{symbol}, could also be
#' \code{entrez} or \code{ENSEMBL}
#' @param topTablerows How many rows to report before any filtering
#' @param fullNamesInRows Logical, whether to display or not the full names for the GO terms
#' @param addGeneToTerms Logical, whether to add a column with all genes annotated to each GO term
#' @param plotGraph Logical, if TRUE additionally plots a graph on the identified GO terms
#' @param plotNodes Number of nodes to plot
#' @param writeOutput Logical, if TRUE additionally writes out the result to a file
#' @param outputFile Name of the file the result should be written into
#' @param topGO_method2 Character, specifying which of the methods implemented by \code{topGO} should be used, in addition to the \code{classic} algorithm. Defaults to \code{elim}
#' @param do_padj Logical, whether to perform the adjustment on the p-values from the specific
#' topGO method, based on the FDR correction. Defaults to FALSE, since the assumption of 
#' independent hypotheses is somewhat violated by the intrinsic DAG-structure of the Gene
#' Ontology Terms
#'
#' @import topGO
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#'
#' @examples
#' library(airway)
#' library(DESeq2)
#' data(airway)
#' airway
#' dds_airway <- DESeqDataSet(airway, design= ~ cell + dex)
#' # Example, performing extraction of enriched functional categories in
#' # detected significantly expressed genes
#' \dontrun{
#' dds_airway <- DESeq(dds_airway)
#' res_airway <- results(dds_airway)
#' library("AnnotationDbi")
#' library("org.Hs.eg.db")
#' res_airway$symbol <- mapIds(org.Hs.eg.db,
#'                             keys = row.names(res_airway),
#'                             column = "SYMBOL",
#'                             keytype = "ENSEMBL",
#'                             multiVals = "first")
#' res_airway$entrez <- mapIds(org.Hs.eg.db,
#'                             keys = row.names(res_airway),
#'                             column = "ENTREZID",
#'                             keytype = "ENSEMBL",
#'                             multiVals = "first")
#' resOrdered <- as.data.frame(res_airway[order(res_airway$padj),])
#' de_df <- resOrdered[resOrdered$padj < .05 & !is.na(resOrdered$padj),]
#' de_symbols <- de_df$symbol
#' bg_ids <- rownames(dds_airway)[rowSums(counts(dds_airway)) > 0]
#' bg_symbols <- mapIds(org.Hs.eg.db,
#'                      keys = bg_ids,
#'                      column = "SYMBOL",
#'                      keytype = "ENSEMBL",
#'                      multiVals = "first")
#' library(topGO)
#' topgoDE_airway <- topGOtable(de_symbols, bg_symbols,
#'                              ontology = "BP",
#'                              mapping = "org.Hs.eg.db",
#'                              geneID = "symbol")
#' }
#'
#' @export
topGOtable <- function(DEgenes,                  # Differentially expressed genes
                       BGgenes,                 # background genes
                       ontology = "BP",            # could use also "MF"
                       annot = annFUN.org,       # parameters for creating topGO object
                       mapping = "org.Mm.eg.db",
                       geneID = "symbol",       # could also beID = "entrez"
                       topTablerows = 200,
                       fullNamesInRows = TRUE,
                       addGeneToTerms = TRUE,
                       plotGraph = FALSE, 
                       plotNodes = 10,
                       writeOutput = FALSE, 
                       outputFile = "",
                       topGO_method2 = "elim",
                       do_padj = FALSE) {
  # checking the additional topGO_method2
  topgo_methods <- c("elim", "weight", "weight01", "lea", "parentchild")
  if (!(topGO_method2 %in% topgo_methods))
     stop("Please provide one of the following topGO methods in addition to the classic method:\n",
          paste(topgo_methods, collapse = ", "))
  # creating the vectors
  DEgenes_input <- factor(as.integer(BGgenes %in% DEgenes))
  names(DEgenes_input) <- BGgenes
  # instantiating the topGOdata object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = DEgenes_input,
                nodeSize = 10,
                annot = annot,
                mapping = mapping,
                ID = geneID)
  # performing the test(s)
  result_method2 <- runTest(GOdata, algorithm = topGO_method2, statistic = "fisher")
  resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  sTab <- GenTable(GOdata,
                   p.value_method2 = result_method2,
                   p.value_classic = resultClassic,
                   orderBy = "p.value_method2",
                   ranksOf = "p.value_classic",
                   topNodes = length(score(resultClassic)))
  
  names(sTab)[which(names(sTab) == "p.value_method2")] <- paste0("p.value_", topGO_method2)
  
  sTab[["p.value_classic"]] <- as.numeric(sTab[["p.value_classic"]])
  sTab[[paste0("p.value_", topGO_method2)]] <- as.numeric(sTab[[paste0("p.value_", topGO_method2)]])
  
  # if FDR, then apply it here
  if (do_padj)
    sTab[[paste0("padj_BY_", topGO_method2)]] <- p.adjust(sTab[[paste0("p.value_", topGO_method2)]], method = "BY")
  
  # subset to specified number of rows
  topTablerows <- min(nrow(sTab), topTablerows)
  sTab <- sTab[seq_len(topTablerows), ]
  
  if (fullNamesInRows) {
    sTab$Term <- sapply(sTab$GO.ID, function(go) {
      Term(GOTERM[[go]])
    })
  }

  if (addGeneToTerms) {
    # adapted from an elegant one liner found here: https://support.bioconductor.org/p/65856/
    SignificantGenes <- sigGenes(GOdata)
    sTab$genes <- sapply(sTab$GO.ID, function(x) {
      genes <- genesInTerm(GOdata, x)
      tmp <- genes[[1]][genes[[1]] %in% SignificantGenes]
    })
    # coerce the list to a comma separated vector
    sTab$genes <- unlist(lapply(sTab$genes, function(arg) paste(arg, collapse = ",")))
  }

  # write all entries of the table
  if (writeOutput) write.table(sTab, file = outputFile, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  if (plotGraph) showSigOfNodes(GOdata, topGO::score(result_method2), firstSigNodes = plotNodes, useInfo = "all")
  #   if(outputToLatex) sTabSig <- xtable(apply(sTabSig[1:15,], 2, as.character)) # take a smaller subset

  # and returns the significant ones # or all, like here
  return(sTab)
}

#' Functional interpretation of the principal components, based on simple
#' overrepresentation analysis
#'
#' Extracts the genes with the highest loadings for each principal component, and
#' performs functional enrichment analysis on them using the simple and quick routine
#' provided by the \code{limma} package
#'
#' @param se A \code{\link{DESeqTransform}} object, with data in \code{assay(se)},
#' produced for example by either \code{\link{rlog}} or
#' \code{\link{varianceStabilizingTransformation}}
#' @param pca_ngenes Number of genes to use for the PCA
#' @param inputType Input format type of the gene identifiers. Deafults to \code{ENSEMBL}, that then will
#' be converted to ENTREZ ids. Can assume values such as \code{ENTREZID},\code{GENENAME} or \code{SYMBOL},
#' like it is normally used with the \code{select} function of \code{AnnotationDbi}
#' @param organism Character abbreviation for the species, using \code{org.XX.eg.db} for annotation
#' @param loadings_ngenes Number of genes to extract the loadings (in each direction)
#' @param background_genes Which genes to consider as background.
#' @param scale Logical, defaults to FALSE, scale values for the PCA
#' @param ... Further parameters to be passed to the goana routine
#'
#' @return A nested list object containing for each principal component the terms enriched
#' in each direction. This object is to be thought in combination with the displaying feature
#' of the main \code{\link{pcaExplorer}} function
#'
#' @examples
#' library(airway)
#' library(DESeq2)
#' library(limma)
#' data(airway)
#' airway
#' dds_airway <- DESeqDataSet(airway, design = ~ cell + dex)
#' \dontrun{
#' rld_airway <- rlogTransformation(dds_airway)
#' goquick_airway <- limmaquickpca2go(rld_airway,
#'                                    pca_ngenes = 10000,
#'                                    inputType = "ENSEMBL",
#'                                    organism = "Hs")
#' }
#'
#' @export
limmaquickpca2go <- function(se,
                             pca_ngenes = 10000,
                             inputType = "ENSEMBL",
                             organism = "Mm",
                             loadings_ngenes = 500,
                             background_genes = NULL,
                             scale = FALSE,
                             ... # further parameters to be passed to the topgo routine
                             ) {
  annopkg <- paste0("org.", organism, ".eg.db")
  if (!require(annopkg, character.only = TRUE)) {
    stop("The package", annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
  }
  exprsData <- assay(se)

  if (is.null(background_genes)) {
    if (is(se, "DESeqDataSet")) {
      BGids <- rownames(se)[rowSums(counts(se)) > 0]
    } else if (is(se, "DESeqTransform")) {
      BGids <- rownames(se)[rowSums(assay(se)) != 0]
    } else {
      BGids <- rownames(se)
    }
  } else {
    BGids <- background_genes
  }

  exprsData <- assay(se)[order(rowVars(assay(se)), decreasing = TRUE), ]
  exprsData <- exprsData[1:pca_ngenes, ]

  rv <- rowVars(exprsData)
  dropped <- sum(rv == 0)
  if (dropped > 0)
    message(paste("Dropped", dropped, "genes with 0 variance"))

  exprsData <- exprsData[rv > 0, ]

  message("After subsetting/filtering for invariant genes, working on a ", nrow(exprsData), "x", ncol(exprsData), " expression matrix\n")

  p <- prcomp(t(exprsData), scale = scale, center = TRUE)

  message("Ranking genes by the loadings ...")
  probesPC1pos <- rankedGeneLoadings(p, pc = 1, decreasing = TRUE)[1:loadings_ngenes]
  probesPC1neg <- rankedGeneLoadings(p, pc = 1, decreasing = FALSE)[1:loadings_ngenes]
  probesPC2pos <- rankedGeneLoadings(p, pc = 2, decreasing = TRUE)[1:loadings_ngenes]
  probesPC2neg <- rankedGeneLoadings(p, pc = 2, decreasing = FALSE)[1:loadings_ngenes]
  probesPC3pos <- rankedGeneLoadings(p, pc = 3, decreasing = TRUE)[1:loadings_ngenes]
  probesPC3neg <- rankedGeneLoadings(p, pc = 3, decreasing = FALSE)[1:loadings_ngenes]
  probesPC4pos <- rankedGeneLoadings(p, pc = 4, decreasing = TRUE)[1:loadings_ngenes]
  probesPC4neg <- rankedGeneLoadings(p, pc = 4, decreasing = FALSE)[1:loadings_ngenes]

  ## convert ensembl to entrez ids
  probesPC1pos_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = probesPC1pos, column = "ENTREZID", keytype = inputType)
  probesPC1neg_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = probesPC1neg, column = "ENTREZID", keytype = inputType)
  probesPC2pos_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = probesPC2pos, column = "ENTREZID", keytype = inputType)
  probesPC2neg_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = probesPC2neg, column = "ENTREZID", keytype = inputType)
  probesPC3pos_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = probesPC3pos, column = "ENTREZID", keytype = inputType)
  probesPC3neg_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = probesPC3neg, column = "ENTREZID", keytype = inputType)
  probesPC4pos_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = probesPC4pos, column = "ENTREZID", keytype = inputType)
  probesPC4neg_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = probesPC4neg, column = "ENTREZID", keytype = inputType)
  bg_ENTREZ <- AnnotationDbi::mapIds(eval(parse(text = annopkg)), keys = BGids, column = "ENTREZID", keytype = inputType)
  message("Ranking genes by the loadings ... done!")

  message("Extracting functional categories enriched in the gene subsets ...")
  quickGOpc1pos <- limma::topGO(limma::goana(probesPC1pos_ENTREZ, bg_ENTREZ, species = organism, ...), ontology = "BP", number = 200); message("1")
  quickGOpc1neg <- limma::topGO(limma::goana(probesPC1neg_ENTREZ, bg_ENTREZ, species = organism, ...), ontology = "BP", number = 200); message("2")
  quickGOpc2pos <- limma::topGO(limma::goana(probesPC2pos_ENTREZ, bg_ENTREZ, species = organism, ...), ontology = "BP", number = 200); message("3")
  quickGOpc2neg <- limma::topGO(limma::goana(probesPC2neg_ENTREZ, bg_ENTREZ, species = organism, ...), ontology = "BP", number = 200); message("4")
  quickGOpc3pos <- limma::topGO(limma::goana(probesPC3pos_ENTREZ, bg_ENTREZ, species = organism, ...), ontology = "BP", number = 200); message("5")
  quickGOpc3neg <- limma::topGO(limma::goana(probesPC3neg_ENTREZ, bg_ENTREZ, species = organism, ...), ontology = "BP", number = 200); message("6")
  quickGOpc4pos <- limma::topGO(limma::goana(probesPC4pos_ENTREZ, bg_ENTREZ, species = organism, ...), ontology = "BP", number = 200); message("7")
  quickGOpc4neg <- limma::topGO(limma::goana(probesPC4neg_ENTREZ, bg_ENTREZ, species = organism, ...), ontology = "BP", number = 200); message("8")

  quickGOpc1pos <- quickGOpc1pos[order(quickGOpc1pos$P.DE), ]
  quickGOpc1neg <- quickGOpc1neg[order(quickGOpc1neg$P.DE), ]
  quickGOpc2pos <- quickGOpc2pos[order(quickGOpc2pos$P.DE), ]
  quickGOpc2neg <- quickGOpc2neg[order(quickGOpc2neg$P.DE), ]
  quickGOpc3pos <- quickGOpc3pos[order(quickGOpc3pos$P.DE), ]
  quickGOpc3neg <- quickGOpc3neg[order(quickGOpc3neg$P.DE), ]
  quickGOpc4pos <- quickGOpc4pos[order(quickGOpc4pos$P.DE), ]
  quickGOpc4neg <- quickGOpc4neg[order(quickGOpc4neg$P.DE), ]

  goEnrichs <- list(PC1 = list(posLoad = quickGOpc1pos, negLoad = quickGOpc1neg),
                    PC2 = list(posLoad = quickGOpc2pos, negLoad = quickGOpc2neg),
                    PC3 = list(posLoad = quickGOpc3pos, negLoad = quickGOpc3neg),
                    PC4 = list(posLoad = quickGOpc4pos, negLoad = quickGOpc4neg)
  )
  message("Extracting functional categories enriched in the gene subsets ... done!")
  attr(goEnrichs, "n_genesforpca") <- pca_ngenes
  return(goEnrichs)
}
