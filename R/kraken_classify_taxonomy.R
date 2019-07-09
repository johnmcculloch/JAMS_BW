#' kraken_classify_taxonomy
#'
#' JAMSalpha function
#' @export

kraken_classify_taxonomy <- function(opt = NULL, fastafile = NULL, confidence = 0.00001){

    #Default is kraken2
    #Count kmers
    flog.info(paste("Counting k-mers of query sequences with kraken2 and using confidence score", as.character(confidence)))
    krakenargs <- c("--db", opt$workingkrakendb, "--threads", opt$threads, "--confidence", confidence, fastafile)
    kraken2taxid <- system2("kraken2", args = krakenargs, stdout = TRUE, stderr = FALSE)
    kraken2taxid <- strsplit(kraken2taxid, split = "[\t]", fixed = FALSE)
    k2out <- plyr::ldply(kraken2taxid, rbind)
    colnames(k2out) <- c("ClassFlag", "Sequence", "Taxid", "Length", "kmers")
    k2out[] <- lapply(k2out, as.character)
    JAMStaxtablefile <- file.path(opt$workingkrakendb, "JAMS_taxtable.tsv")
    if (file.exists(JAMStaxtablefile)){
        JAMStaxtable <- read.table(file = JAMStaxtablefile, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = NULL, header = TRUE)
    } else {
        #Fall back on generic taxonomy table and warn user
        flog.info("JAMS taxonomy table not found. Falling back on generic JAMS taxtable.")
        data(JAMStaxtable)
    }
    JAMStaxtable[] <- lapply(JAMStaxtable, as.character)
    krakendf <- k2out[, c("Sequence", "Taxid")]
    krakendf <- left_join(krakendf, JAMStaxtable)
    #Fill in taxids which are NOT in the database with unclassifieds
    unclassdf <- JAMStaxtable[which(JAMStaxtable$Taxid == 0), 2:ncol(JAMStaxtable)]
    krakendf[which(is.na(krakendf$Domain == TRUE)), colnames(unclassdf)] <- unname(unclassdf)

    return(krakendf)
}
