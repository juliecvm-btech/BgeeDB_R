#' @title Download the TSV files needed for a TopAnat analysis
#'
#' @description Downloads (or reuses cached copies of) the three TSV files needed to run a
#' TopAnat analysis: anatomical entity relationships, anatomical entity names, and the mapping
#' of genes to anatomical entities. Called by \code{loadTopAnatData()}.
#'
#' @details If a file for the requested species/parameters is already present in
#' \code{myBgeeObject$pathToData}, it is reused and not redownloaded.
#'
#' @param myBgeeObject An output object from Bgee$new().
#'
#' @param callType A character indicating the type of expression calls to be used for
#' enrichment. Only calls for significant detection of expression are implemented so far
#' ("presence").
#'
#' @param confidence A character indicating if only high quality present calls should be
#' retrieved. For Bgee releases prior to 14, options are "all" or "high_quality". For Bgee
#' release 14 and above, options are "silver" or "gold".
#'
#' @param stage A character indicating the targeted developmental stage(s) for the analysis, or
#' NULL to retrieve expression patterns regardless of developmental stage.
#'
#' @param timeout local timeout used when the function is run, in seconds. See
#' \code{loadTopAnatData()} for details.
#'
#' @param OLD_WEBSERVICE_VERSION A character giving the last Bgee webservice version (e.g.
#' "13.2") using the older query API. Used to decide which URL format and confidence values to
#' use.
#'
#' @return A list of 3 elements, giving the local file paths of the (possibly newly downloaded)
#' TSV files:
#' \itemize{
#'   \item{\code{organRelationshipsFilePath}}
#'   \item{\code{organNamesFilePath}}
#'   \item{\code{gene2anatomyFilePath}}
#' }
#'
#' @author Julien Roux, Julien Wollbrett
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Danio_rerio", dataType = "rna_seq")
#'   filePaths <- downloadTopAnatFiles(myBgeeObject = bgee, callType = "presence",
#'                                      confidence = "silver", stage = NULL, timeout = 1800,
#'                                      OLD_WEBSERVICE_VERSION = "13.2")
#' }
#'
#' @export
downloadTopAnatFiles <- function(myBgeeObject, callType, confidence, stage, timeout, OLD_WEBSERVICE_VERSION){
  ## Set the timeout option to timeout value to let some time to the server to send data (default is 1800 sec.)
  op <- options(timeout = timeout)
  ## on exit change back options to initial values
  on.exit(options(op))

  ## First query: organ relationships
  organRelationshipsFileName <- paste0("topAnat_AnatEntitiesRelationships_", myBgeeObject$speciesId, ".tsv")
  ## Check if file is already in cache
  if (file.exists(file.path(myBgeeObject$pathToData, organRelationshipsFileName))){
    cat(paste0("\nNOTE: an organ relationships file was found in the download directory ", myBgeeObject$pathToData,
        ". Data will not be redownloaded.\n"))
  } else {
    cat("\nBuilding URLs to retrieve organ relationships from Bgee.........\n")
    myUrl <- myBgeeObject$topAnatUrl
    if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
      myUrl <- paste0(myUrl, "?page=r_package&action=get_anat_entity_relations&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=SOURCE_ID&attr_list=TARGET_ID&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    } else {
      myUrl <- paste0(myUrl, "?page=dao&action=org.bgee.model.dao.api.ontologycommon.RelationDAO.getAnatEntityRelations&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=SOURCE_ID&attr_list=TARGET_ID&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }
    ## Query webservice
    cat(paste0("   URL successfully built (", myUrl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    success <- bgee_download_file(url = myUrl, destfile = paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName, ".tmp"))

    if (success == 0){
      ## Read 5 last lines of file: should be empty indicating success of data transmission
      ## We cannot use a system call to UNIX command since some user might be on Windows
      tmp <- tail(read.table(paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
      if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
        ## The file transfer was successful, we rename the temporary file
        file.rename(paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName, ".tmp"), paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName))
      } else {
        ## delete the temporary file
        file.remove(paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName, ".tmp"))
        stop(paste0("File ", organRelationshipsFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
      }
      cat(paste0("   Got results from Bgee webservice. Files are written in \"", myBgeeObject$pathToData, "\"\n"))
    } else {
      serverAnswer = try(getURL(myUrl))
      if (class(serverAnswer) == "try-error"){
        stop("ERROR: the query to the server was not successful. Is your internet connection working?\n")
      } else {
        stop(paste0("ERROR: the query to the server was not successful. The server returned the following answer:\n", serverAnswer))
      }
    }
  }

  ## Second query: organ names
  organNamesFileName <- paste0("topAnat_AnatEntitiesNames_", myBgeeObject$speciesId, ".tsv");
  ## Check if file is already in cache
  if (file.exists(file.path(myBgeeObject$pathToData, organNamesFileName))){
    cat(paste0("\nNOTE: an organ names file was found in the download directory ", myBgeeObject$pathToData,
               ". Data will not be redownloaded.\n"))

  } else {
    cat("\nBuilding URLs to retrieve organ names from Bgee.................\n")
    myUrl <- myBgeeObject$topAnatUrl
    if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
      myUrl <- paste0(myUrl, "?page=r_package&action=get_anat_entities&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=ID&attr_list=NAME&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }else {
      myUrl <- paste0(myUrl, "?page=dao&action=org.bgee.model.dao.api.anatdev.AnatEntityDAO.getAnatEntities&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=ID&attr_list=NAME&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }

    ## Query webservice
    cat(paste0("   URL successfully built (", myUrl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    success <- bgee_download_file(url = myUrl, destfile = paste0(myBgeeObject$pathToData, "/", organNamesFileName, ".tmp"))

    if (success == 0){
      ## Read 5 last lines of file: should be empty indicating success of data transmission
      ## We cannot use a system call to UNIX command since some user might be on Windows
      tmp <- tail(read.table(paste0(myBgeeObject$pathToData, "/", organNamesFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE, quote = ""), n=5)
      if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
        ## The file transfer was successful, we rename the temporary file
        file.rename(paste0(myBgeeObject$pathToData, "/", organNamesFileName, ".tmp"), paste0(myBgeeObject$pathToData, "/", organNamesFileName))
      } else {
        ## delete the temporary file
        file.remove(paste0(myBgeeObject$pathToData, "/", organNamesFileName, ".tmp"))
        stop(paste0("File ", organNamesFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
      }
     cat(paste0("   Got results from Bgee webservice. Files are written in \"", myBgeeObject$pathToData, "\"\n"))
    } else {
      serverAnswer = try(getURL(myUrl))
      if (class(serverAnswer) == "try-error"){
        stop("ERROR: the query to the server was not successful. Is your internet connection working?\n")
      } else {
        stop(paste0("ERROR: the query to the server was not successful. The server returned the following answer:\n", serverAnswer))
      }
    }
  }

  ## Third query: gene to organs mapping

  # The Java API does not distinguish between full length and droplet based single cell. A topAnat analysis can be run
  # on all single cell data but not on data coming from a subset of single cell technologies.
  # Bgee objects have been designed to allow the download of expression files for any datatype and then distinguish between
  # full length and droplet based single cell.
  # In order to solve that mismatch we update the type used to run topAnat analysis. If either sc_full_length or sc_droplet_based
  # datatype is selected, we run a topAnat analysis including all single cell technologies (full length AND droplet based)

  # First write a warning if only one single cell technology is selected
  if ("sc_full_length" %in% myBgeeObject$dataType & ! "sc_droplet_based" %in% myBgeeObject$dataType |
    "sc_droplet_based" %in% myBgeeObject$dataType & ! "sc_full_length" %in% myBgeeObject$dataType) {
    message("WARNING: TopAnat can not be run on one single cell technology. Both full length and droplet based single cell data will",
      " be queried for this topAnat analysis. If you do not want to query single cell data please remove \"sc_full_length\" or \"sc_droplet_based\"",
      " from the list of datatypes of your Bgee object.")
  }

  # Then update the list of datatypes used to run topAnat
  topAnat_dataType <- myBgeeObject$dataType
  if ("sc_full_length" %in% topAnat_dataType | "sc_droplet_based" %in% topAnat_dataType) {
    topAnat_dataType <- topAnat_dataType[! topAnat_dataType %in% c("sc_full_length", "sc_droplet_based")]
    topAnat_dataType <- append(topAnat_dataType, "sc_rna_seq")
  }

  gene2anatomyFileName <- paste0("topAnat_GeneToAnatEntities_", myBgeeObject$speciesId, "_", toupper(callType))
  ## If a stage is specified, add it to file name
  if ( !is.null(stage) ){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, "_", gsub(":", "_", stage))
  }
  ## If all data types specified, no need to add anything to file name. Otherwise, specify data types in file name
  if ( sum(topAnat_dataType %in% c("rna_seq","affymetrix","est","in_situ", "sc_rna_seq")) < 5 ){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, "_", toupper(paste(sort(topAnat_dataType), collapse="_")))
  }
  ## If high quality data needed, specify in file name. Otherwise not specified
  if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, "_", toupper(confidence))
  } else {
    if ( confidence == "high_quality" ){
      gene2anatomyFileName <- paste0(gene2anatomyFileName, "_HIGH")
    }
  }
  gene2anatomyFileName <- paste0(gene2anatomyFileName, ".tsv")

  ## Check if file is already in cache
  if (file.exists(file.path(myBgeeObject$pathToData, gene2anatomyFileName))){
    cat(paste0("\nNOTE: a gene to organs mapping file was found in the download directory ", myBgeeObject$pathToData,
               ". Data will not be redownloaded.\n"))

  } else {
    cat("\nBuilding URLs to retrieve mapping of gene to organs from Bgee...\n")
    myUrl <- myBgeeObject$topAnatUrl
    if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
      myUrl <- paste0(myBgeeObject$topAnatUrl, "?page=r_package&action=get_expression_calls&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=GENE_ID&attr_list=ANAT_ENTITY_ID&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }else {
      myUrl <- paste0(myUrl, "?page=dao&action=org.bgee.model.dao.api.expressiondata.ExpressionCallDAO.getExpressionCalls&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=GENE_ID&attr_list=ANAT_ENTITY_ID&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }

    ## Add data type to file name: only if not all data types asked
    if ( sum(topAnat_dataType %in% c("rna_seq","sc_rna_seq","affymetrix","est","in_situ")) < 5 ){
      for (type in toupper(sort(topAnat_dataType))){
        myUrl <- paste0(myUrl, "&data_type=", type)
      }
    }
    ## Add data quality
    if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
      myUrl <- paste0(myUrl, "&data_qual=", toupper(confidence))
    } else {
      if(confidence == "high_quality"){
        myUrl <- paste0(myUrl, "&data_qual=HIGH")
      }
    }

    if ( !is.null(stage) ){
      myUrl <- paste0(myUrl, "&stage_id=", stage)
    }

    ## Query webservice
    cat(paste0("   URL successfully built (", myUrl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    ## this download correspond to a file that can either be 1) generated on the fly if it is the first time this combination of conditions
    ## is queried or 2) stored on our server otherwise.
    ## If the file has to be generated it can take a lot of time (sometimes more than one hour if single cell data exist). In order to solve
    ## download errors due to apache connection stopping while the file is still generating and then throw an error 502, we decided to check
    ## the error retrieved by that download and restart the download if a specific error was retrieved.
    catchedError <- 1
    while (catchedError) {
      startTime <- Sys.time()
      tryCatch(
        {
          bgee_download_file(url = myUrl, destfile = paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName, ".tmp"))
          catchedError <- 0
        },
        error = function(x) {
          #If it took a lot of time to throw an error there is a high probability that the file is currently generated.
          # In that case we just wait one minute and try to download the file again
          # by default the timeout is 60 seconds. In order to be safe we check that the error took more than 50
          # 40 sec to be thrown
          stopTime <- Sys.time()
          if (as.numeric(difftime(time1 = stopTime, time2 = startTime, units = "sec")) <= 40) {
            stop(paste0("ERROR: the query to the server was not successful. The server returned the following answer:\n", x))
          }
        },
        warning = function(x) {}
      )
    }

    tmp <- tail(read.table(paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
    if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
      ## The file transfer was successful, we rename the temporary file
      file.rename(paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName, ".tmp"), paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName))
    } else {
      ## delete the temporary file
      file.remove(paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName, ".tmp"))
      stop(paste0("File ", gene2anatomyFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
    cat(paste0("   Got results from Bgee webservice. Files are written in \"", myBgeeObject$pathToData, "\"\n"))
  }

  return(list(
    organRelationshipsFilePath = file.path(myBgeeObject$pathToData, organRelationshipsFileName),
    organNamesFilePath = file.path(myBgeeObject$pathToData, organNamesFileName),
    gene2anatomyFilePath = file.path(myBgeeObject$pathToData, gene2anatomyFileName)
  ))
}
