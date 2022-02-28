# pcaExplorer 2.22.0

## Other notes

* `get_annotation_orgdb()` gains an additional argument, `key_for_genenames`, which defaults to "SYMBOL". This should not change the behavior of the function, if not specified, but accommodates for the use of annotation packages where the information has been encoded differently (e.g. org.Sc.sgd.db where the info is contained in the "ORF" column)

# pcaExplorer 2.20.0

## Other notes

* the tables in the PCA2GO tab panel can be compacted only if they are computed via the `pca2go` function (offline) - at runtime, `limmaquickpca2go` is used and no compaction is required
* if an annotation is provided with a column `gene_id`, these values are actually overwriting the rownames (makes the object more robust with respect to its provenance)

# pcaExplorer 2.16.0

## Other notes

* Replaced dependency from `d3heatmap` with the functionality of `heatmaply`

# pcaExplorer 2.12.0

## Bug fixes

* Fixed an error in the initialization of the app due to a new behavior introduced by `shinyAce` in version >= 0.4.0
* `topGOtable` does not generate rows with NAs if providing a too high number for the categories to report

## Other notes

* The type of the columns in the data.frame returned by `topGOtable` are now correctly referring to the type they contain - e.g. the p values are now stored as numeric values
* Citation now refers to the published manuscript - https://doi.org/10.1186/s12859-019-2879-1

# pcaExplorer 2.10.0

## New features

* Added extra parameters to `topGOtable` to offer more control on the method used, and the option to correct the p-values for multiple testing (via the Benjamini-Yekutieli procedure)
* `pca2go` has now an option to return (early) a list with vectors of genes with high loadings
* Better preview of the uploaded data with modal dialog windows triggered by buttons which appear once corresponding inputs are available
* Improved notification system: in addition to the progress bar, info that all input is correctly specified, or suggest optional annotation loading
* Added flexibility to select separator for each of the uploadable files
* The pairwise correlation plots can now use logarithmic scale in the axes, use smaller subsets of the data for quicker inspection, and resizes the correlation info proportionally to its intensity
* The sample to sample heatmap supports additionally manhattan and correlation-based distances
* There is a new vignette with a detailed quick start, "Up and running with pcaExplorer", specifying how the workflow with `pcaExplorer` can look like, demonstrated on the `airway` dataset
* In the Instructions panel, we added buttons to access the fully rendered documentation, either local or online if e.g. deployed to a server. Related to this, `pcaExplorer` has a new parameter, `runLocal`, to control this behavior
* An additional parameter, `annopkg`, has been added to `pca2go()` to override the behavior with the `organism` parameter (this is useful when the name of the annotation package is not conform to the classical `org.Xx.eg.db`, e.g. for Arabidopsis Thaliana); a detailed use case has been added in the main vignette

## Other notes

* The computing of the required objects now requires the explicit action on the dedicated button, and the tooltip informs the user on what steps are taken (including normalization)
* An information box has been added to provide detailed information on the required input formats
* Added notification to specify how to install the airway package for demonstration purposes if not already available
* Added startup message upon loading the package
* The content in the Instructions tab is now contained in collapsible elements
* The file formats accepted by `pcaExplorer` are now specified both in the vignette text, as well as in the app at runtime
* The content of the Instructions tab is now more compact, containing the rendered "Up and running with pcaExplorer" vignette. The full vignettes can be accessed via buttons in the same panel
* Added instructions to install phantomJS via the `webshot` package - would raise an error when previewing the report

# pcaExplorer 2.8.0

## New features

* Added a `NEWS.md` file to track changes to the package
* PCA plots now are correctly generated with fixed coordinates
* Introduced use of conditionalPanels for better handling of errors in the app tabs
* Added possibility to use different transformations, also reflected in the change of one of the main arguments (previously `rlt`, now `dst`, i.e. `DESeqTransform`): rlog, vst, shifted log, ... The transformation type is tracked in the reactive values. 
* More modular loading of data, by splitting generation of `dds` and `dst`
* `pca2go` is now also picking values from the input widgets

## Other notes

* Built project website via pkgdown, with customized reference structure
* Correctly adding the resources to shinyBS, loaded via `.onLoad`, and also better placement for bstooltips
* Editor options start collapsed in the Report Editor tab
* Vignette and template report are updated to reflect the new parameter names
* Uniformed style for ggplot2 plots
* Better tooltip placement in the main page
* Replaced `print` calls with more appropriate `message`s
* Displaying user returned messages in long (plotting) operations

## Bug fixes

* Fixed behavior of rendering inline the content of the report - did not work properly for server deployed instances

# pcaExplorer 2.6.0

## New features

* Automatically computing size factors where required
* Added progress indication when compiling the report

## Bug fixes

* Fixed after changes in threejs package
* Edited dropdown menu to remove unused green badge
* Menus start expanded on the side, again
* `theme_bw` applied when needed, corrected previous behavior

## Other notes

* Updated citation infos
* Slight difference in handling validate/need errors

# pcaExplorer 2.2.0

## New features

* Added Demo data, loadable via demo button

## Bug fixes

* Plots work now without cutting out points when zooming in

## Other notes

* Saved reactive values are now exported to dedicate environments (instead of assigning to global)

# pcaExplorer 1.99.0

## Other notes

* Reflecting the major feature added, will trigger a major version number bump. Welcome soon, pcaExplorer 2.0.0!

# pcaExplorer 1.1.5

## New features

* Automated report generation - template available + editor in the app tab for advance user customization
* Support for state saving, in the global environment as well as with binary data
* All plots generated can be now exported with the dedicated button 
* Added confidence ellipse for PCA plot
* Added 3d pca plot
* Added functions to automatically retrieve the annotation in format ready to use for the app
* Added profile explorer function, for plotting behaviour across all samples for subset of genes
* Added distribution plots
* Added pairwise correlation plot
* Added table to enhance readability of the gene finder plot, also by annotating sample names

## Bug fixes

* Minor typos fixed in the tabs
* Added option row.names to read.delim for allowing row names when uploading the data

## Other notes

* Added extra info in the about section
* Instructions and vignette rewritten to reflect new design of the app

# pcaExplorer 1.1.3

## Bug fixes

* Remove y axis limits to gene boxplots
* Fixed: correct labels and colors assignements for genespca

# pcaExplorer 1.0.0

## Other notes

* Released in Bioconductor 3.3

# pcaExplorer 0.99.1

## Other notes

* Changed format of the NEWS file

# pcaExplorer 0.99.0

## Other notes

* Ready for submission to Bioconductor

# pcaExplorer 0.9.0

## Other notes

* Added TravisCI integration for both branches
* Added appveyor integration - plus badges in the README.md
* Code cleanup
* Added screenshots for the vignette
* Removed some lengthy tests

# pcaExplorer 0.8.0

## New features

* Selection of identifier type available in pca2go

## Bug fixes

* Couple of layout fixes

## Other notes

* MIT license
* Added TravisCI integration
* Added codecov integration
* Enhanced documentation

# pcaExplorer 0.7.0

## New features

* Vignette full draft done

# pcaExplorer 0.6.4

## Other notes

* Updated NEWS file

# pcaExplorer 0.6.3

## New features

* About and Instructions done by now
* Added some missing details on the documentations

# pcaExplorer 0.6.2

## Other notes

* Corrected wordings for (cor)relations of principal components with covariates
* Added a couple of checks if correct objects are provided

# pcaExplorer 0.6.1

## New features

* Added function to remove selected samples suspected to be deemed as outliers, in order to see the effect of clustering on the good ones

# pcaExplorer 0.6.0

## Other notes

* Documentation completed
* Examples fully working, cleaned up further a little more.

# pcaExplorer 0.5.0

## Other notes

* Further steps in direction of R CMD check

# pcaExplorer 0.4.0

## New features

* Added pca2go live functionality

# pcaExplorer 0.3.0

## New features

* Added color palette to choose, and dependent on the samples and factors available/selected

# pcaExplorer 0.2.0

## New features

* Multifactorial exploration completed and adaptable to each dataset

# pcaExplorer 0.1.0

## New features

* Restyling and (re)packaging mostly completed
