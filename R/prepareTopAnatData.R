#' @title Parse downloaded TopAnat TSV files into the objects needed for analysis
#'
#' @description Reads the three TSV files downloaded by \code{downloadTopAnatFiles()} and
#' builds the list of objects needed to run a TopAnat analysis: a gene-to-anatomy mapping, the
#' anatomical ontology relationships (as a DAG, rooted at a single "BGEE:0" root), and the
#' anatomical entity names. Called by \code{loadTopAnatData()}.
#'
#' @param organRelationshipsFilePath Local path to the TSV file listing the relationships
#' between anatomical structures in the UBERON ontology (columns SOURCE_ID, TARGET_ID).
#'
#' @param organNamesFilePath Local path to the TSV file listing anatomical structure names
#' (columns ID, NAME).
#'
#' @param gene2anatomyFilePath Local path to the TSV file listing gene-to-anatomical-structure
#' expression calls (columns GENE_ID, ANAT_ENTITY_ID).
#'
#' @return A list of 3 elements:
#' \itemize{
#'   \item{A \code{gene2anatomy} list, mapping genes to anatomical structures based on expression calls.}
#'   \item{An \code{organ.relationships} list, giving the relationships between anatomical structures in the UBERON ontology (based on parent-child "is_a" and "part_of" relationships).}
#'   \item{An \code{organ.names} data frame, with the name corresponding to UBERON IDs.}
#' }
#'
#' @author Julien Roux, Julien Wollbrett
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Danio_rerio", dataType = "rna_seq")
#'   filePaths <- downloadTopAnatFiles(myBgeeObject = bgee, callType = "presence",
#'                                      confidence = "silver", stage = NULL, timeout = 1800,
#'                                      OLD_WEBSERVICE_VERSION = "13.2")
#'   myTopAnatData <- prepareTopAnatData(organRelationshipsFilePath = filePaths$organRelationshipsFilePath,
#'                                        organNamesFilePath = filePaths$organNamesFilePath,
#'                                        gene2anatomyFilePath = filePaths$gene2anatomyFilePath)
#' }
#'
#' @export
prepareTopAnatData <- function(organRelationshipsFilePath, organNamesFilePath, gene2anatomyFilePath){
  ## Process the data and build the final list to return
  cat("\nParsing the results.............................................\n")

  ## Relationships between organs
  if (file.exists(organRelationshipsFilePath)){
    if (file.info(organRelationshipsFilePath)$size != 0) {
      tab <- read.table(organRelationshipsFilePath, header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
      organRelationships <- tapply(as.character(tab$TARGET_ID), as.character(tab$SOURCE_ID), unique)
    } else {
      stop(paste0("File ", basename(organRelationshipsFilePath), " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
  } else {
    stop(paste0("File ", basename(organRelationshipsFilePath), " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  ## Organ names
  if (file.exists(organNamesFilePath)){
    if (file.info(organNamesFilePath)$size != 0) {
      organNames <- read.table(organNamesFilePath, header=TRUE, sep="\t", comment.char="", blank.lines.skip=TRUE, as.is=TRUE, quote = "")
    } else {
      stop(paste0("File ", basename(organNamesFilePath), " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
  } else {
    stop(paste0("File ", basename(organNamesFilePath), " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  ## Mapping of genes to tissues
  if (file.exists(gene2anatomyFilePath)){
    if (file.info(gene2anatomyFilePath)$size != 0) {
      tab <- read.table(gene2anatomyFilePath, header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
      if(length(tab$GENE_ID) != 0){
        gene2anatomy <- tapply(as.character(tab$ANAT_ENTITY_ID), as.character(tab$GENE_ID), unique)
      } else {
        stop("There was no mapping of genes to anatomical structures found. Probably the parameters are too stringent, or this data type is absent in this species. See listBgeeSpecies() for data types availability.")
      }
    } else {
      stop(paste0("File ", basename(gene2anatomyFilePath), " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters. It is also possible that the parameters are too stringent and returned no data, please try to relax them."))
    }
  } else {
    stop(paste0("File ", basename(gene2anatomyFilePath), " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters. It is also possible that the parameters are too stringent and returned no data, please try to relax them."))
  }

  cat("\nAdding BGEE:0 as unique root of all terms of the ontology.......\n")
  ## There can be multiple roots among all the terms downloaded. We need to add one unique root for topGO to work: BGEE:0
  ## Add all organs from organNames that are not source (child / names of the list) in organsRelationship to the organsRelationship list (with value / target / parent = BGEE:0)
  ## Some organs are not present in the organRelationships file because they have no relations to other organs (not linked to a target).
  ## That's why we use all organs present in organNames to find missingParents
  missingParents <- setdiff(organNames[, 1], names(organRelationships))

  ## Add new values
  organRelationships <- c(organRelationships, as.list(rep("BGEE:0", times=length(missingParents))))
  ## Add new keys
  names(organRelationships)[(length(organRelationships) - length(missingParents) + 1):length(organRelationships)] = as.character(missingParents)
  ## Add BGEE:0	/ root to organNames
  organNames <- rbind(organNames, c("BGEE:0", "root"))

  ## Check if some organ names are missing, and add them if necessary
  missingNames <- setdiff(unique(unique(unlist(organRelationships, use.names = FALSE)), names(organRelationships)), organNames$ID)
  if (length(missingNames) > 0){
    cat(paste0("\nWARNING: some organs names appear to be missing. There might be some problem with the ontology data.\n"))
    organNames <- rbind(organNames, setNames(data.frame(missingNames, "?"), names(organNames)))
  }

  cat("\nDone.\n")
  return(list(gene2anatomy = gene2anatomy, organ.relationships = organRelationships, organ.names = organNames))
}
