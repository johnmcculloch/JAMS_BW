#' find_accession_from_keyterm(expvec = NULL, analyses = NULL, SEobj = NULL, keyterms = NULL)
#'
#' Returns the accession number(s) of a feature(s) matching a general expression present in a SummarizedExperiment or metagenomeSeq object.
#' @export

find_accession_from_keyterm <- function(expvec = NULL, analyses = NULL, SEobj = NULL, keyterms = NULL){

    Features_of_Interest_df <- NULL

    SEobj_list <- list()
    if (!is.null(expvec)){
        if (!is.null(analyses)){
            wantedanals <- analyses[(analyses %in% names(expvec))]
            SEobj_list <- expvec[wantedanals]
        } else {
            SEobj_list <- expvec
        }
    } else {
        SEobj_list[1] <- SEobj
    }

    for (analnum in 1:length(SEobj_list)){
        currobj <- SEobj_list[[analnum]]
        if (as.character(class(currobj)[1]) != "SummarizedExperiment"){
            flog.warn("This function can only take a SummarizedExperiment object as input.")

            return(NULL)
        }
        #Obtain the feature table
        analysisname <- metadata(currobj)$analysis
        if (analysisname == "LKT"){
            flog.warn("Searching feature accessions and descriptions by keyterm will only work for non-taxonomic SummarizedExperiment objects.")
        } else {
            ftt <- rowData(currobj)
            #include accession and description
            ftt$Feature <- paste(ftt$Accession, ftt$Description, sep="-")
            #Get info by general expression matching
            for (kt in keyterms){
                tmp_ftt <- ftt[grep(kt, ftt$Feature, ignore.case = TRUE), c("Accession", "Description")]
                if (nrow(tmp_ftt) > 0){
                    tmp_ftt$Analysis <- analysisname
                    tmp_ftt <- tmp_ftt[ , c("Analysis", "Accession", "Description")]
                    Features_of_Interest_df <- rbind(Features_of_Interest_df, as.data.frame(tmp_ftt))
                    rownames(Features_of_Interest_df) <- 1:nrow(Features_of_Interest_df)
                } else {
                    flog.warn(paste("No", kt, "found within", analysisname))
                }
            }
        }
    }

    return(Features_of_Interest_df)

}
