#' Get an annotation data frame from biomaRt
#'
#' @param dds A \code{\link{DESeqDataSet}} object
#' @param biomart_dataset A biomaRt dataset to use. To see the list, type
#' \code{mart = useMart('ensembl')}, followed by \code{listDatasets(mart)}.
#' @param idtype Character, the ID type of the genes as in the row names of
#' \code{dds}, to be used for the call to \code{\link{getBM}}
#'
#' @return A data frame for ready use in \code{pcaExplorer}, retrieved from biomaRt.
#' @export
#'
#' @examples
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'                                              colData = colData(airway),
#'                                              design = ~dex+cell)
#' \dontrun{
#' get_annotation(dds_airway, "hsapiens_gene_ensembl", "ensembl_gene_id")
#' }
get_annotation <- function(dds, biomart_dataset, idtype) {
  if (is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")

  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  host = "www.ensembl.org",
                  dataset = biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)

  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id", "gene_name", "description")

  return(anns2)
}


#' Get an annotation data frame from org db packages
#'
#' @param dds A \code{\link{DESeqDataSet}} object
#' @param orgdb_species Character string, named as the \code{org.XX.eg.db}
#' package which should be available in Bioconductor
#' @param idtype Character, the ID type of the genes as in the row names of
#' \code{dds}, to be used for the call to \code{\link{mapIds}}
#' @param key_for_genenames Character, corresponding to the column name for the 
#' key in the orgDb package containing the official gene name (often called 
#' gene symbol). 
#' This parameter defaults to "SYMBOL", but can be adjusted in case the key is not
#' found in the annotation package (e.g. for \code{org.Sc.sgd.db}).
#'
#' @return A data frame for ready use in \code{pcaExplorer}, retrieved from the
#' org db packages
#' @export
#'
#' @examples
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'                                              colData = colData(airway),
#'                                              design = ~dex+cell)
#' anno_df <- get_annotation_orgdb(dds_airway, "org.Hs.eg.db", "ENSEMBL")
#' head(anno_df)
get_annotation_orgdb <- function(dds, orgdb_species, idtype, key_for_genenames = "SYMBOL") {
  if (is.null(orgdb_species))
    stop("Select a species to generate the corresponding annotation")

  orgdbpkgs <- data.frame(
    pkg = c("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
            "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Hs.ipi.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pf.plasmo.db",
            "org.Pt.eg.db", "org.Rn.eg.db", "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db", "org.Xl.eg.db"),
    descr = c("Anopheles", "Arabidopsis", "Bovine", "Worm", "Canine", "Fly", "Zebrafish", "E coli strain K12", "E coli strain Sakai", "Chicken",
              "Human", "org.Hs.ipi.db", "Mouse", "Rhesus", "Malaria", "Chimp", "Rat", "Yeast", "Streptomyces coelicolor", "Pig", "Toxoplasma gondii",
              "Xenopus"),
    stringsAsFactors = FALSE
  )

  if (!(orgdb_species %in% orgdbpkgs$pkg)) {
    message("The orgDB package is most likely not existent in Bioconductor")
    message("It should be one of", orgdbpkgs$pkg)
  }

  if (!require(orgdb_species, character.only = TRUE))
    stop("The package ", orgdb_species, " is not installed/available. Try installing it with BiocManager::install('", orgdb_species, "')")

  if (!(idtype %in% keytypes(eval(parse(text = orgdb_species))))) {
    stop("The key you provided is not listed as key for the annotation package. Please try one of ",
         paste(keytypes(eval(parse(text = orgdb_species))), collapse = ","))
  }

  if (!(key_for_genenames %in% keytypes(eval(parse(text = orgdb_species))))) {
    stop("The key specified for containing gene names is not included in the annotation package. Please try one of ",
         paste(keytypes(eval(parse(text = orgdb_species))), collapse = ","))
  }
  
  
  pkg <- eval(parse(text = orgdb_species))

  if (idtype == "SYMBOL")
    warning("You probably do not need to convert symbol to symbol") # the performance would somehow be affected

  anns_vec <- mapIds(pkg, keys = rownames(dds), column = key_for_genenames,
                     keytype = idtype)

  anns <- data.frame(
    gene_id = rownames(dds),
    gene_name = anns_vec,
    stringsAsFactors = FALSE,
    row.names = rownames(dds)
  )
  return(anns)
}
