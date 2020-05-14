#' make_featuredata_from_bedfile
#'
#' JAMSalpha function
#' @export

make_featuredata_from_bedfile <- function(opt = NULL, bedfile = NULL){

    if (is.null(bedfile)){
        bedfile <- opt$bedfile
    }

    featuredata <- fread(file = bedfile, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, sep = "\t")

    colnames(featuredata) <- c("Contig", "Start", "End", "Feature", "MapQual", "Strand", "Annotby", "FeatType", "Spin", "Annot")

    featuredata <- subset(featuredata, FeatType %in% c("CDS", "tRNA",  "rRNA",  "tmRNA"))

    annots <- tidyr::separate_rows(featuredata[ , c("Feature", "Annot")], all_of("Annot"), sep = fixed("\\;")) %>% tidyr::separate("Annot", c("Annot_type", "Annot_value"), sep = "=", remove = TRUE, convert = FALSE)
    annots <- subset(annots, Annot_type %in% c("product", "Name", "eC_number"))
    annots <- tidyr::spread(annots, "Annot_type", "Annot_value", fill = "none", convert = FALSE)

    colnames(annots) <- sapply(colnames(annots), function (x) { switch(x, "Feature" = "Feature", "product" = "Product", "Name" = "GeneNameNR", "eC_number" = "ECNumber") })

    featuredata <- left_join(featuredata, annots, by = "Feature")

    featuredata$GeneName <- sapply(1:length(featuredata$GeneNameNR), function (x) { unlist(strsplit(featuredata$GeneNameNR[x], split = "_"))[1] } )

    featuredata$ECNumber[which(featuredata$ECNumber != "none")] <- paste("EC", featuredata$ECNumber[which(featuredata$ECNumber != "none")], sep = "_")

    featuredata$LengthDNA <- (as.numeric(featuredata$End) - as.numeric(featuredata$Start))

    featuredata <- featuredata[ , c("Feature", "LengthDNA", "FeatType", "Contig", "Product", "GeneName", "GeneNameNR", "ECNumber")]

    #Trim whitespace
    featuredata <- as.data.frame(featuredata)
    featuredata <- trim_whitespace_from_df(df = featuredata)

    return(featuredata)
}
