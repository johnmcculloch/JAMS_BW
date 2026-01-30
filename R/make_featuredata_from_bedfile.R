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

    #There may be no gene names, if that is the case, the GeneNameNR will be missing.
    if ("GeneName" %in% colnames(featuredata)){
        featuredata$GeneName <- sapply(1:length(featuredata$GeneNameNR), function (x) { unlist(strsplit(featuredata$GeneNameNR[x], split = "_"))[1] } )
    } else {
        featuredata$GeneName <- ""
        featuredata$GeneNameNR <- ""
    }

    #There may be no enzymes, if that is the case, then ECNumber will be missing.
    if ("ECNumber" %in% colnames(featuredata)){
        featuredata$ECNumber[which(featuredata$ECNumber != "none")] <- paste("EC", featuredata$ECNumber[which(featuredata$ECNumber != "none")], sep = "_")
    } else {
        featuredata$ECNumber <- "none"
    }

    featuredata$LengthDNA <- (as.numeric(featuredata$End) - as.numeric(featuredata$Start))

    featuredata <- featuredata[ , c("Feature", "LengthDNA", "FeatType", "Contig", "Product", "GeneName", "GeneNameNR", "ECNumber")]

    #Trim whitespace
    featuredata <- as.data.frame(featuredata)
    featuredata <- trim_whitespace_from_df(df = featuredata)

    return(featuredata)
}
