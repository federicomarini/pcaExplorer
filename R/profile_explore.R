#' Extract and plot the expression profile of genes
#'
#' @param se A [DESeq2::DESeqDataSet()] object, or a
#' [DESeq2::DESeqTransform()] object.
#' @param genelist An array of characters, including the names of the genes of
#' interest of which the profile is to be plotted
#' @param intgroup A factor, needs to be in the `colnames` of `colData(se)`
#' @param plotZ Logical, whether to plot the scaled expression values. Defaults to
#' `FALSE`
#'
#' @return A plot of the expression profile for the genes
#' @export
#'
#' @examples
#' dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3, betaSD_tissue = 1)
#' rlt <- DESeq2::rlogTransformation(dds)
#' geneprofiler(rlt, paste0("gene", sample(1:1000, 20)))
#' geneprofiler(rlt, paste0("gene", sample(1:1000, 20)), plotZ = TRUE)
geneprofiler <- function(se, genelist = NULL, intgroup = "condition", plotZ = FALSE) {
  if (is.null(genelist))
    stop("Provide at least one gene to the genelist parameter")
  # check that at least one gene is found
  genelist <- unique(genelist)
  message("You provided ", length(genelist), " unique identifiers")
  inthedata <- genelist %in% rownames(se)
  if (sum(inthedata) == 0)
    stop("None of the provided genes were found in the experiment data")
  message(sum(inthedata), " out of ", length(genelist), " provided genes were found in the data")

  mydata <- as.data.frame(t(assay(se)[genelist, ]))

  # resort the order of the rows according to the groups that are selected
  mygroups <- interaction(as.data.frame(colData(se)[intgroup]))
  mydata <- mydata[order(mygroups), ]

  if (plotZ) {
    # remove 0 variance genes
    rv <- rowVars(t(mydata))
    mydata <- mydata[, rv > 0]

    mydata <- scale(mydata, center = TRUE, scale=TRUE)
    # was...
    # mydata <- NMF:::scale_mat(mydata,"col")
  }
  mylabels <- colnames(se)[order(mygroups)]
  mycols <- scales::hue_pal()(length(levels(mygroups)))[sort(mygroups)]

  par(mar=c(7.1, 4.1, 2.1, 2.1))
  plot(mydata[, 1], type = "l", xaxt = "n", las = 2, ylim = range(mydata), xlab = "", ylab = ifelse(plotZ, "scaled expression value", "expression value"))
  Map(function(x, y, z)
    axis(1, at = x, col.axis = y, labels = z, lwd = 0, las = 2),
    1:nrow(mydata),
    mycols,
    mylabels
  )
  axis(1, at = 1:nrow(mydata), labels = FALSE)

  for (i in 2:(ncol(mydata) - 1)){
      lines(mydata[, i], type = "l", xaxt = "n", las = 2, col = i)
  }
  ## TODO: if desired, plot only the avg pro group -> maybe as boxplot?
}
