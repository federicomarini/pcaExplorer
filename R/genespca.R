#' Principal components analysis on the genes
#'
#' Computes and plots the principal components of the genes, eventually displaying
#' the samples as in a typical biplot visualization.
#'
#' The implementation of this function is based on the beautiful `ggbiplot`
#' package developed by Vince Vu, available at https://github.com/vqv/ggbiplot.
#' The adaptation and additional parameters are tailored to display typical genomics data
#' such as the transformed counts of RNA-seq experiments
#'
#' @param x A [DESeq2::DESeqTransform()] object, with data in `assay(x)`,
#' produced for example by either [DESeq2::rlog()] or
#' [DESeq2::varianceStabilizingTransformation()]
#' @param ntop Number of top genes to use for principal components,
#' selected by highest row variance
#' @param choices Vector of two numeric values, to select on which principal components to plot
#' @param arrowColors Vector of character, either as long as the number of the samples, or one single value
#' @param groupNames Factor containing the groupings for the input data. Is efficiently chosen
#' as the (interaction of more) factors in the colData for the object provided
#' @param biplot Logical, whether to additionally draw the samples labels as in a biplot representation
#' @param scale Covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1,
#' the inner product between the variables approximates the covariance and the
#' distance between the points approximates the Mahalanobis distance.
#' @param pc.biplot Logical, for compatibility with biplot.princomp()
#' @param obs.scale Scale factor to apply to observations
#' @param var.scale Scale factor to apply to variables
#' @param groups Optional factor variable indicating the groups that the observations
#' belong to. If provided the points will be colored according to groups
#' @param ellipse Logical, draw a normal data ellipse for each group
#' @param ellipse.prob Size of the ellipse in Normal probability
#' @param labels optional Vector of labels for the observations
#' @param labels.size Size of the text used for the labels
#' @param alpha Alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param var.axes Logical, draw arrows for the variables?
#' @param circle Logical, draw a correlation circle? (only applies when prcomp
#' was called with scale = TRUE and when var.scale = 1)
#' @param circle.prob Size of the correlation circle in Normal probability
#' @param varname.size Size of the text for variable names
#' @param varname.adjust  Adjustment factor the placement of the variable names,
#' '>= 1' means farther from the arrow
#' @param varname.abbrev  Logical, whether or not to abbreviate the variable names
#' @param returnData Logical, if TRUE returns a data.frame for further use, containing the
#' selected principal components for custom plotting
#' @param coordEqual Logical, default FALSE, for allowing brushing. If TRUE, plot using
#' equal scale cartesian coordinates
#' @param scaleArrow Multiplicative factor, usually >=1, only for visualization purposes,
#' to allow for distinguishing where the variables are plotted
#' @param useRownamesAsLabels Logical, if TRUE uses the row names as labels for plotting
#' @param point_size Size of the points to be plotted for the observations (genes)
#' @param annotation A `data.frame` object, with row.names as gene identifiers (e.g. ENSEMBL ids)
#' and a column, `gene_name`, containing e.g. HGNC-based gene symbols
#'
#' @return An object created by `ggplot`, which can be assigned and further customized.
#'
#' @examples
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3, betaSD_tissue = 1)
#' rlt <- rlogTransformation(dds)
#' groups <- colData(dds)$condition
#' groups <- factor(groups, levels = unique(groups))
#' cols <- scales::hue_pal()(2)[groups]
#' genespca(rlt, ntop=100, arrowColors = cols, groupNames = groups)
#'
#' groups_multi <- interaction(as.data.frame(colData(rlt)[, c("condition", "tissue")]))
#' groups_multi <- factor(groups_multi, levels = unique(groups_multi))
#' cols_multi <- scales::hue_pal()(length(levels(groups_multi)))[factor(groups_multi)]
#' genespca(rlt, ntop = 100, arrowColors = cols_multi, groupNames = groups_multi)
#'
#' @export
genespca <- function(x,
                     ntop,
                     choices = c(1, 2),
                     arrowColors = "steelblue",
                     groupNames="group",
                     biplot = TRUE,
                     scale = 1, pc.biplot = TRUE,
                     obs.scale = 1 - scale, var.scale = scale, groups = NULL,
                     ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3,
                     alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69,
                     varname.size = 4, varname.adjust = 1.5, varname.abbrev = FALSE,
                     returnData = FALSE, coordEqual = FALSE, scaleArrow = 1,
                     useRownamesAsLabels = TRUE, point_size = 2, annotation = NULL) {

  stopifnot(length(choices) == 2)
  if (length(arrowColors) != 1 & length(arrowColors) != ncol(x))
    stop("Please provide either one color or a vector as long as the number of samples")

  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp((assay(x)[select, ]))

  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  if (!biplot) {
    nobs.factor <- sqrt(nrow(pca$x) - 1)
    devs <- pca$sdev
    pcast <- pca
    pcast$x <- sweep(pca$x, 2, 1 / (devs * nobs.factor), FUN = "*") * nobs.factor
    d <- data.frame(PC1 = pcast$x[, choices[1]],
                    PC2 = pcast$x[, choices[2]],
                    names = rownames((assay(x)[select, ])))

    if (returnData) {
      attr(d, "percentVar") <- percentVar
      return(d)
    }

    ggplot(data = d, aes_string(x = "PC1", y = "PC2")) +
      geom_point(size = 3) +
      xlab(paste0("PC", choices[1], ": ", round(percentVar[choices[1]] * 100), "% variance")) +
      ylab(paste0("PC", choices[2], ": ", round(percentVar[choices[2]] * 100), "% variance")) +
      # geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F) +
      ggtitle("title") + theme_bw()
  } else {
    if (inherits(pca, "prcomp")) {
      nobs.factor <- sqrt(nrow(pca$x) - 1)
      d <- pca$sdev
      u <- sweep(pca$x, 2, 1 / (d * nobs.factor), FUN = "*")
      v <- pca$rotation
    }

    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale,
                                FUN = "*"))
    v <- sweep(v, 2, d^var.scale, FUN = "*")
    df.v <- as.data.frame(v[, choices])
    names(df.u) <- c("xvar", "yvar")
    names(df.v) <- names(df.u)
    if (pc.biplot) {
      df.u <- df.u * nobs.factor
    }

    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2)) ^ (1/4)
    v.scale <- rowSums(v^2)
    df.v <- r * df.v / sqrt(max(v.scale))
    if (obs.scale == 0) {
      u.axis.labs <- paste("standardized PC", choices, sep = "")
    } else {
      u.axis.labs <- paste("PC", choices, sep = "")
    }
    u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)",
                                              100 * pca$sdev[choices]^2 / sum(pca$sdev^2)))
    if (!is.null(labels)) {
      df.u$labels <- labels
    }
    if (!is.null(groups)) {
      df.u$groups <- groups
    }

    # additionally...
    df.u$ids <- rownames(df.u)
    if(!is.null(annotation)) {
      df.u$geneNames <- annotation$gene_name[match(df.u$ids, rownames(annotation))]
    } else {
      df.u$geneNames <- df.u$ids
    }
    if (varname.abbrev) {
      df.v$varname <- abbreviate(rownames(v))
    } else {
      df.v$varname <- rownames(v)
    }
    df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
    df.v$hjust <- with(df.v, (1 - varname.adjust * sign(xvar)) / 2)

    if (returnData) {
      return(df.u)
    }



    g <- ggplot(data = df.u, aes_string(x = "xvar", y = "yvar")) + xlab(u.axis.labs[1]) +
      ylab(u.axis.labs[2]) # + coord_equal() # REMOVED OTHERWISE BRUSH DOES NOT WORK PROPERLY
    if (coordEqual) g <- g + coord_equal()

    if (!is.null(df.u$labels)) {
      if (!is.null(df.u$groups)) {
        g <- g + geom_text(aes(label = labels, color = groups),
                           size = labels.size)
      } else {
        g <- g + geom_text(aes(label = labels), size = labels.size)
      }
    } else {
      if (!is.null(df.u$groups)) {
        g <- g + geom_point(aes(color = groups), size = point_size, alpha = alpha)
      } else {
        g <- g + geom_point(size = point_size, alpha = alpha)
      }
    }

    if (useRownamesAsLabels) {
      g <- g + geom_text(aes_string(label = "geneNames"), size = labels.size, hjust = 0.25, vjust = -0.75)
    }

    if (!is.null(df.u$groups) && ellipse) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- cbind(cos(theta), sin(theta))
      ell <- ddply(df.u, "groups", function(x) {
        if (nrow(x) <= 2) {
          return(NULL)
        }
        sigma <- var(cbind(x$xvar, x$yvar))
        mu <- c(mean(x$xvar), mean(x$yvar))
        ed <- sqrt(qchisq(ellipse.prob, df = 2))
        data.frame(sweep(circle %*% chol(sigma) * ed, 2,
                         mu, FUN = "+"), groups = x$groups[1])
      })
      names(ell)[1:2] <- c("xvar", "yvar")
      g <- g + geom_path(data = ell, aes(color = groups, group = groups))
    }
    # moved down to have the arrows drawn on top of the points and not vice versa
    if (var.axes) {
      if (circle) {
        theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi,
                                                  length = 50))
        circle <- data.frame(xvar = r * cos(theta), yvar = r *
                               sin(theta))
        g <- g + geom_path(data = circle, color = "steelblue",
                           size = 1/2, alpha = 1/3)
      }
      df.v$scaleArrow <- scaleArrow # quick fix for mapping scaling of the arrows
      arrowColors <-  factor(arrowColors, levels = unique(arrowColors))
      df.v$arrowColors <- factor(arrowColors, levels = unique(arrowColors))
      df.v$groupNames <- factor(groupNames, levels = unique(groupNames))
      df.v$sca_x <- df.v$xvar * scaleArrow
      df.v$sca_y <- df.v$yvar * scaleArrow
      df.v$sta_x <- 0
      df.v$sta_y <- 0
      g <- g + geom_segment(data = df.v, aes_string(x = "sta_x", y = "sta_y", xend = "sca_x", yend = "sca_y", color = "arrowColors"),
                            arrow = arrow(length = unit(1/2, "picas"))) +
        scale_color_manual(values = levels(arrowColors), name = "Group", labels = levels(groupNames))
    }

    if (var.axes) {
      g <- g + geom_text(data = df.v, aes_string(label = "varname",
                                                 x = "sca_x", y = "sca_y", # angle = angle,
                                                 hjust = "hjust"),
                         color = arrowColors, size = varname.size)
    }
    g <- g + theme_bw()
    return(g)
  }
}
