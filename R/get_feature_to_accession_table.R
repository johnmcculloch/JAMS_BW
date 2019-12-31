#' get_feature_to_accession_table
#'
#' JAMSalpha function
#' @export

get_feature_to_accession_table <- function(opt = NULL, iproanalysis = NULL){

        flog.info(paste("Adding", iproanalysis, "signatures to featuredata."))
        #subset ipro to contain only applicable analysis
        if (!(iproanalysis %in% c("Interpro", "GO"))){
            iprointerest <- subset(opt$interproscanoutput, Analysis == iproanalysis)
            accessioncol <- "Accession"
            descriptioncol <- "Description"
        } else {
            iprointerest <- subset(opt$interproscanoutput, IproAcc != "none")
            accessioncol <- "IproAcc"
            descriptioncol <- "IproDesc"
        }

        if(iproanalysis != "GO"){
            featsIwant<-sapply(unique(iprointerest[, "Feature"]), function (x) { paste0(sort(unique(iprointerest[which(iprointerest[,"Feature"] == x), accessioncol])), collapse="|")} )
        } else {
            #get rid of information without GO terms
            iprointerest<-subset(iprointerest, GOterms != "none")
            #If looking for GO terms, split up accessions in GOterms column, and get sorted, non-reundant list.
            featsIwant<-sapply(unique(iprointerest[,"Feature"]), function (x) { paste0(sort(unique(unlist(strsplit(iprointerest[which(iprointerest[,"Feature"] == x), "GOterms"], "|", fixed=TRUE)))), collapse="|")} )
        }

    feat2acc<-data.frame(Feature=names(featsIwant), Accession=unname(featsIwant), stringsAsFactors = FALSE)
    colnames(feat2acc)[2]<-iproanalysis

    return(feat2acc)
}
