#' Sample PCA plot for transformed data
#'
#' Plots the results of PCA on a 2-dimensional space
#'
#' @param x A \code{\link{DESeqTransform}} object, with data in \code{assay(x)},
#' produced for example by either \code{\link{rlog}} or
#' \code{\link{varianceStabilizingTransformation}}
#' @param intgroup Interesting groups: a character vector of
#' names in \code{colData(x)} to use for grouping
#' @param ntop Number of top genes to use for principal components,
#' selected by highest row variance
#' @param returnData logical, if TRUE returns a data.frame for further use, containing the
#' selected principal components and intgroup covariates for custom plotting
#' @param title The plot title
#' @param pcX The principal component to display on the x axis
#' @param pcY The principal component to display on the y axis
#' @param text_labels Logical, whether to display the labels with the sample identifiers
#' @param point_size Integer, the size of the points for the samples
#' @param ellipse Logical, whether to display the confidence ellipse for the selected groups
#' @param ellipse.prob Numeric, a value in the interval [0;1)
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @examples
#' dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3, betaSD_tissue = 1)
#' rlt <- DESeq2::rlogTransformation(dds)
#' pcaplot(rlt, ntop = 200)
#'
#' @export
pcaplot <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE,title = NULL,
                    pcX = 1, pcY = 2, text_labels = TRUE, point_size = 3,
                    ellipse = TRUE, ellipse.prob = 0.95) # customized principal components
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))

  percentVar <- pca$sdev^2/sum(pca$sdev^2)

  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(x)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC1 = pca$x[, pcX], PC2 = pca$x[, pcY], group = group,
                  intgroup.df, names = colnames(x))
  colnames(d)[1] <- paste0("PC", pcX)
  colnames(d)[2] <- paste0("PC", pcY)

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }

  # clever way of positioning the labels - worked good, then no need with ggrepel
  d$hjust <- ifelse((sign(d[, paste0("PC", pcX)]) == 1), 0.9, 0.1)# (1 + varname.adjust * sign(PC1))/2)

  g <- ggplot(data = d, aes_string(x = paste0("PC", pcX), y = paste0("PC", pcY), color = "group")) +
    geom_point(size = point_size) +
    xlab(paste0("PC", pcX, ": ", round(percentVar[pcX] * 100, digits = 2), "% variance")) +
    ylab(paste0("PC", pcY ,": ", round(percentVar[pcY] * 100, digits = 2), "% variance"))

  ## plot confidence ellipse
  # credit to vince vu, author of ggbiplot
  if (ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))

    ell <- ddply(d, "group", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x[[paste0("PC", pcX)]], x[[paste0("PC", pcY)]]))
      mu <- c(mean(x[[paste0("PC", pcX)]]), mean(x[[paste0("PC", pcY)]]))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$group[1])
    })
    # names(ell)[1:2] <- c('xvar', 'yvar')
    if (nrow(ell) > 0) {
      g <- g + geom_path(data = ell, aes_string(x = "X1", y = "X2", color = "groups", group = "groups"))
    }
  }

  if (text_labels)
    g <- g + geom_label_repel(mapping = aes_string(label = "names", fill = "group"),
                              color = "white", show.legend = TRUE)

  plot_title <- paste0("PCA plot - top ", ntop, " variable genes")
  if (!is.null(title)) {
    g <- g + ggtitle(title)
  } else {
    g <- g + ggtitle(plot_title)
  }
  g <- g + theme_bw()
  # as in http://www.huber.embl.de/msmb/Chap-Graphics.html
  # "well-made PCA plots usually have a width that’s larger than the height"
  g <- g + coord_fixed()
  g
}


#' Scree plot of the PCA on the samples
#'
#' Produces a scree plot for investigating the proportion of explained variance, or
#' alternatively the cumulative value
#'
#' @param obj A \code{prcomp} object
#' @param type Display absolute proportions or cumulative proportion. Possible values:
#' "pev" or "cev"
#' @param pc_nr How many principal components to display max
#' @param title Title of the plot
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @examples
#' dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3, betaSD_tissue = 1)
#' rlt <- DESeq2::rlogTransformation(dds)
#' pcaobj <- prcomp(t(SummarizedExperiment::assay(rlt)))
#' pcascree(pcaobj, type = "pev")
#' pcascree(pcaobj, type = "cev", title = "Cumulative explained proportion of variance - Test dataset")
#'
#' @export
pcascree <- function(obj, type = c("pev", "cev"),pc_nr=NULL,title=NULL)
{
  type <- match.arg(type)
  d <- obj$sdev^2
  yvar <- switch(type, pev = d/sum(d), cev = cumsum(d)/sum(d))
  yvar.lab <- switch(type, pev = "proportion of explained variance",
                     cev = "cumulative proportion of explained variance")
  # df <- data.frame(PC = 1:length(d), yvar = yvar)

  if (!is.null(pc_nr)) {
    colsize <- pc_nr
    yvar <- yvar[1:pc_nr]
  } else {
    colsize <- length(d)
    yvar <- yvar[1:length(d)]
  }

  pc_df <- data.frame(PC_count = 1:colsize, var = yvar)

  if(type=="pev"){
    p <- ggplot(pc_df, aes_string(x = "PC_count", y = "var")) + geom_bar(stat = "identity")
    p <- p + scale_x_continuous(breaks = 1:length(d))
    p <- p + ylab(yvar.lab) + xlab("principal components")
    # p
  } else {
    p <- ggplot(pc_df, aes_string(x = "PC_count", y = "var")) + geom_point() +
      geom_path() + scale_x_continuous(breaks = 1:length(d))
    p <- p + ylab(yvar.lab) + xlab("principal components") + ylim(0,max(pc_df$var))
    # p
  }
  if(!is.null(title)) p <- p + ggtitle(title)
  p <- p + theme_bw()
  p
}


#' Sample PCA plot for transformed data
#'
#' Plots the results of PCA on a 3-dimensional space, interactively
#'
#' @param x A \code{\link{DESeqTransform}} object, with data in \code{assay(x)},
#' produced for example by either \code{\link{rlog}} or
#' \code{\link{varianceStabilizingTransformation}}
#' @param intgroup Interesting groups: a character vector of
#' names in \code{colData(x)} to use for grouping
#' @param ntop Number of top genes to use for principal components,
#' selected by highest row variance
#' @param returnData logical, if TRUE returns a data.frame for further use, containing the
#' selected principal components and intgroup covariates for custom plotting
#' @param title The plot title
#' @param pcX The principal component to display on the x axis
#' @param pcY The principal component to display on the y axis
#' @param pcZ The principal component to display on the z axis
#' @param text_labels Logical, whether to display the labels with the sample identifiers
#' @param point_size Integer, the size of the points for the samples
#'
#' @return A html-based visualization of the 3d PCA plot
#' @export
#'
#' @examples
#' dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3, betaSD_tissue = 1)
#' rlt <- DESeq2::rlogTransformation(dds)
#' pcaplot3d(rlt, ntop = 200)
pcaplot3d <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE,title=NULL,
                     pcX = 1, pcY = 2, pcZ = 3, text_labels=TRUE,point_size=3)
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))

  percentVar <- pca$sdev^2/sum(pca$sdev^2)

  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(x)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC1 = pca$x[, pcX], PC2 = pca$x[, pcY], PC3 = pca$x[,pcZ],
                  group = group,
                  intgroup.df, names = colnames(x))
  colnames(d)[1] <- paste0("PC", pcX, ": ", round(percentVar[pcX] * 100,digits = 2), "% variance")
  colnames(d)[2] <- paste0("PC", pcY, ": ", round(percentVar[pcY] * 100,digits = 2), "% variance")
  colnames(d)[3] <- paste0("PC", pcZ, ": ", round(percentVar[pcZ] * 100,digits = 2), "% variance")

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }

  nrgroups <- length(levels(d$group))
  cols <- hue_pal()(nrgroups)[d$group]

  scatterplot3js(as.matrix(d[, 1:3]),
                 color = cols,
                 # renderer = "canvas",
                 size = 1.3,
                 labels = rownames(d), label.margin ="50px 50px 50px 50px")
}
