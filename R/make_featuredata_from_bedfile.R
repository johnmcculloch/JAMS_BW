#' make_featuredata_from_bedfile
#'
#' JAMSalpha function
#' @export

make_featuredata_from_bedfile <- function(opt = NULL, bedfile = NULL){

    if(is.null(bedfile)){
        bedfile <- opt$bedfile
    }

    featuredata <- fread(file = bedfile, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

    colnames(featuredata) <- c("Contig", "Start", "End", "Feature", "MapQual", "Strand", "Annotby", "FeatType", "Spin", "Annot")

    featuredata <- subset(beddata, FeatType %in% c("CDS", "tRNA",  "rRNA",  "tmRNA"))


    split_annot <- function(annot){
        annotations <- unlist(strsplit(annot, split = ";"))

        ftypes <- sapply(1:length(annotations), function (x) { unlist(strsplit(annotations[x], split = "="))[1] })
        fvalues <- sapply(1:length(annotations), function (x) { unlist(strsplit(annotations[x], split = "="))[2] })
        annotdf <- data.frame(AnnotFeatType = ftypes, AnnotFeatValues = fvalues, stringsAsFactors = FALSE)


        return (annotdf)
    }

    annot_list <- lapply(featuredata$Annot, function (x) { split_annot(x) })
    names(annot_list) <- featuredata$Feature

    retrieve_annot_value <- function (Feature = NULL, FeatType = NULL){
        annotvalue <- annot_list[[Feature]][which(annot_list[[Feature]]$AnnotFeatType == FeatType), "AnnotFeatValues"]

        if (length(annotvalue) != 1){
            annotvalue <- "none"
        }

        return(annotvalue)
    }

    featuredata$Product <- sapply(1:nrow(featuredata), function (x) { retrieve_annot_value(Feature = featuredata$Feature[x], FeatType = "product")}  )

    featuredata$GeneNameNR <- sapply(1:nrow(featuredata), function (x) { retrieve_annot_value(Feature = featuredata$Feature[x], FeatType = "gene")}  )

    featuredata$GeneName <- sapply(1:length(featuredata$GeneNameNR), function (x) { unlist(strsplit(featuredata$GeneNameNR[x], split = "_"))[1] } )

    featuredata$ECNumber <- sapply(1:nrow(featuredata), function (x) { retrieve_annot_value(Feature = featuredata$Feature[x], FeatType = "eC_number")}  )

    featuredata$ECNumber[which(featuredata$ECNumber != "none")] <- paste("EC", featuredata$ECNumber[which(featuredata$ECNumber != "none")], sep = "_")

    featuredata$LengthDNA <- (as.numeric(featuredata$End) - as.numeric(featuredata$Start))

    featuredata <- featuredata[ , c("Feature", "LengthDNA", "FeatType", "Contig", "Product", "GeneName", "GeneNameNR", "ECNumber")]

    #Trim whitespace
    featuredata <- as.data.frame(featuredata)
    featuredata <- trim_whitespace_from_df(df = featuredata)

    return(featuredata)
}
