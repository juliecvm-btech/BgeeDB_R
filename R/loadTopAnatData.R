#' @title Retrieve data from Bgee to perform GO-like enrichment of anatomical terms, mapped to genes by expression patterns.
#'
#' @description This function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology.
#'
#' @details The expression calls come from Bgee (\url{http://bgee.org}), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) from multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference atlas of normal gene expression. Anatomical structures are identified using IDs from the Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}). The mapping from genes to anatomical structures includes only the evidence of expression in these specific structures, and not the expression in their substructures (i.e., expression data are not propagated). The retrieval of propagated expression data might be implemented in the future, but meanwhile, it can be obtained using specialized packages such as topGO, see the \code{topAnat.R} function.
#'
#' @param myBgeeObject An output object from Bgee$new().
#'
#' @param callType A character of indicating the type of expression calls to be used for enrichment. Only calls for significant detection of expression are implemented so far ("presence"). Differential expression calls, based on differential expression analysis, might be implemented in the future.
#'
#' @param stage A character indicating the targeted developmental stages for the analysis. Developmental stages can be chosen from the developmental stage ontology used in Bgee (available at \url{https://github.com/obophenotype/developmental-stage-ontologies}). If a stage is specified, the expression pattern mapped to this stage and all children developmental stages (substages) will be retrieved. Default is NULL, meaning that expression patterns of genes are retrieved regardless of the developmental stage displaying expression; this is equivalent to specifying stage="UBERON:0000104" (life cycle, the root of the stage ontology).
#' For information, the most useful stages (going no deeper than level 3 of the ontology) include:
#' \itemize{
#'   \item{UBERON:0000068 (embryo stage)}
#'   \itemize{
#'     \item{UBERON:0000106 (zygote stage)}
#'     \item{UBERON:0000107 (cleavage stage)}
#'     \item{UBERON:0000108 (blastula stage)}
#'     \item{UBERON:0000109 (gastrula stage)}
#'     \item{UBERON:0000110 (neurula stage)}
#'     \item{UBERON:0000111 (organogenesis stage)}
#'     \item{UBERON:0007220 (late embryonic stage)}
#'     \item{UBERON:0004707 (pharyngula stage)}
#'   }
#'   \item{UBERON:0000092 (post-embryonic stage)}
#'   \itemize{
#'     \item{UBERON:0000069 (larval stage)}
#'     \item{UBERON:0000070 (pupal stage)}
#'     \item{UBERON:0000066 (fully formed stage)}
#'   }
#' }
#'
#' @param confidence A character indicating if only high quality present calls should be
#' retrieved. For Bgee releases prior to 14, options are "all" (default) or "high_quality".
#' For Bgee release 14 and above, options are "silver" (default) and "gold".
#'
#' @param timeout local timeout used when the function is run. It allows to modify the
#' timeout used by default in R to download files. If not provided the timeout will be
#' fixed at 1800 secondes. When the functions exits (either naturally or as the result
#' of an error) the download timeout go back to the original one used in the R session.
#'
#' @return A list of 4 elements:
#' \itemize{
#'   \item{A \code{gene2anatomy} list, mapping genes to anatomical structures based on expression calls.}
#'   \item{A \code{organ.names} data frame, with the name corresponding to UBERON IDs.}
#'   \item{A \code{organ.relationships} list, giving the relationships between anatomical structures in the UBERON ontology (based on parent-child "is_a" and "part_of" relationships).}
#'   \item{The Bgee class object thta was used to retrieve the data.}
#' }
#'
#' @author Julien Roux, Julien Wollbrett
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Danio_rerio", dataType = "rna_seq")
#'   myTopAnatData <- loadTopAnatData(bgee)
#' }
#'
#' @import utils digest
#' @export
loadTopAnatData <- function(myBgeeObject, callType="presence", confidence=NULL, stage=NULL, timeout = 1800){
  OLD_WEBSERVICE_VERSION = '13.2'

  ## check that fields of Bgee object are not empty
  if (length(myBgeeObject$speciesId) == 0 | length(myBgeeObject$topAnatUrl) == 0 | length(myBgeeObject$dataType) == 0 | length(myBgeeObject$pathToData) == 0 | length(myBgeeObject$sendStats) == 0){
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object is valid.")
  }
  if ( callType != "presence" ){
    stop("ERROR: no other call types than \"presence\" expression calls can be retrieved for now.")
  }
  if ( compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0 ){
    if ( is.null(confidence) ){
      confidence = "silver"
    }
    if ( (confidence != "gold") && (confidence != "silver") ){
      stop(paste0("ERROR: the data confidence parameter specified is not among the allowed values. For Bgee ", myBgeeObject$release, " allowed values are \"silver\" or \"gold\".\nBy default \"silver\" quality is selected."))
    }
  } else {
    if(is.null(confidence)){
      confidence = "all"
    }
    if ( (confidence != "all") && (confidence != "high_quality") ){
      stop(paste0("ERROR: the data confidence parameter specified is not among the allowed values. For Bgee ", myBgeeObject$release, " allowed values are \"all\" or \"high_quality\".\nBy default \"all\" quality is selected."))
    }
  }

  ## Download (or retrieve from cache) the three files needed for TopAnat analysis
  filePaths <- downloadTopAnatFiles(myBgeeObject = myBgeeObject, callType = callType, confidence = confidence, stage = stage, timeout = timeout, OLD_WEBSERVICE_VERSION = OLD_WEBSERVICE_VERSION)

  ## Read the downloaded files and build the objects needed for TopAnat analysis
  preparedData <- prepareTopAnatData(organRelationshipsFilePath = filePaths$organRelationshipsFilePath,
                                      organNamesFilePath = filePaths$organNamesFilePath,
                                      gene2anatomyFilePath = filePaths$gene2anatomyFilePath)
  preparedData$bgee.object <- myBgeeObject
  return(preparedData)
}

