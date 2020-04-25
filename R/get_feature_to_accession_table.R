#' get_feature_to_accession_table
#'
#' JAMSalpha function
#' @export

get_feature_to_accession_table <- function(opt = NULL, iproanalysis = NULL){

        flog.info(paste("Adding", iproanalysis, "signatures to featuredata."))
        #subset ipro to contain only applicable analysis
        if (!(iproanalysis %in% c("Interpro", "GO", "MetaCyc"))){
            iprointerest <- subset(opt$interproscanoutput, Analysis == iproanalysis)
            accessioncol <- "Accession"
            descriptioncol <- "Description"
        } else if (iproanalysis == "MetaCyc"){
            iprointerest <- subset(opt$interproscanoutput, Pathways != "none")
            accessioncol <- "Accession"
            descriptioncol <- "Description"
        } else {
            iprointerest <- subset(opt$interproscanoutput, IproAcc != "none")
            accessioncol <- "IproAcc"
            descriptioncol <- "IproDesc"
        }

        #Declare useful funcitons
        fish_pathway <- function (FeatInterest = NULL, iprointerest = NULL, pathwayspace = NULL) {
            nr_annots <- unique(unlist(strsplit(iprointerest[which(iprointerest[ , "Feature"] == FeatInterest), "Pathways"], "|", fixed = TRUE)))
            nr_annots <- nr_annots[grep(pathwayspace, nr_annots)]
            if (length(nr_annots) < 1){
                nr_annots <- "none"
            }
            nr_annots <- gsub("MetaCyc: ", "", nr_annots)
            nr_annots <- paste0(sort(unique(nr_annots)), collapse = "|")
            return(nr_annots)
        }

        if (iproanalysis == "GO"){
            #get rid of information without GO terms
            iprointerest <- subset(iprointerest, GOterms != "none")
            #If looking for GO terms, split up accessions in GOterms column, and get sorted, non-reundant list.
            featsIwant <- sapply(unique(iprointerest[ , "Feature"]), function (x) { paste0(sort(unique(unlist(strsplit(iprointerest[which(iprointerest[ , "Feature"] == x), "GOterms"], "|", fixed = TRUE)))), collapse = "|")} )
        } else if (iproanalysis == "MetaCyc"){
             #get rid of information without GO terms
            iprointerest <- subset(iprointerest, Pathways != "")
            #If looking for pathway terms, split up annotations in Pathways column, and get sorted, non-reundant list.
            featsIwant <- sapply(unique(iprointerest[ , "Feature"]), function (x) { fish_pathway(FeatInterest = x, iprointerest = iprointerest, pathwayspace = iproanalysis) } )
        } else {
            featsIwant <- sapply(unique(iprointerest[, "Feature"]), function (x) { paste0(sort(unique(iprointerest[which(iprointerest[,"Feature"] == x), accessioncol])), collapse = "|")} )
        }

    feat2acc <- data.frame(Feature = names(featsIwant), Accession = unname(featsIwant), stringsAsFactors = FALSE)
    #Remove any non-informative information
    feat2acc <- feat2acc[which(!(feat2acc$Accession %in% c("none", ""))), ]
    colnames(feat2acc)[2] <- iproanalysis

    return(feat2acc)
}
