#' Plot distribution of expression values
#'
#' @param rld A [DESeq2::DESeqTransform()] object.
#' @param plot_type Character, choose one of `boxplot`, `violin` or
#' `density`. Defaults to `density`
#'
#' @return A plot with the distribution of the expression values
#' @export
#'
#' @examples
#' dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3, betaSD_tissue = 1)
#' rlt <- DESeq2::rlogTransformation(dds)
#' distro_expr(rlt)
distro_expr <- function(rld, plot_type="density") {
  allrld <- tidyr::gather(as.data.frame(assay(rld)))
  names(allrld) <- c("Sample", "rlogExpression")

  if (plot_type == "boxplot") {
    p <- ggplot(allrld, aes_string(x = "Sample", y = "rlogExpression")) + geom_boxplot(aes_string(col = "Sample", fill = "Sample"), alpha = 0.5)
  }

  if (plot_type == "violin") {
    p <- ggplot(allrld, aes_string(x = "Sample", y = "rlogExpression")) + geom_violin(aes_string(col = "Sample", fill = "Sample"), alpha = 0.5)
  }

  if (plot_type == "density") {
    p <- ggplot(allrld, aes_string(x = "rlogExpression")) + geom_density(aes_string(color = "Sample"), alpha = 0.1)
  }
  p <- p + theme_bw()
  return(p)
}
