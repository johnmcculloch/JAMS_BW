#' plot_wordcloud_sample(featuredose_df = NULL, analysis = NULL, removeunclassifieds = TRUE)
#'
#' JAMSalpha function
#' @export

#Plot wordcloud of functional data
plot_wordcloud_sample <- function(featuredose_df = NULL, analysis = NULL, removeunclassifieds = TRUE){

    if (!(analysis %in% featuredose_df$Analysis)){
        flog.info(paste("There is no analysis named", analysis, "available in the featuredose object for this sample."))
    } else {
        analrelabund <- subset(featuredose_df, Analysis == analysis)
        analrelabund <- analrelabund[, c("Accession", "Description", "NumBases")]
        totbases <- sum(analrelabund$NumBases)
        analrelabund$PPM <- round((analrelabund$NumBases / totbases) * 1000000, 0)
        dffunc <- analrelabund[ , c("Description", "PPM")]

        if ((length(unique(dffunc$Description)) == 1) && (unique(dffunc$Description) == "none")){
            dffunc$Description <- rownames(dffunc)
        }

        #Eliminate unclassified signatures
        if (removeunclassifieds == TRUE){
            dffunc <- subset(dffunc, !(Description %in% c(paste(analysis, "none", sep = "_"), "none")))
        }

        if (nrow(dffunc) > 0){
            colnames(dffunc) <- c("Description", "Freq")
            dffunc <- subset(dffunc, Freq > 0)
            dffunc <- dffunc[order(dffunc$Freq, decreasing = TRUE), ]
            dffunc$CumRep <- cumsum(dffunc$Freq)/sum(dffunc$Freq)
            smallest <- min(dffunc$Freq)
            dffunc$Freq <- round((dffunc$Freq/smallest), 0)
            if (nrow(dffunc) > 200){
                dffunc <- dffunc[1:200, ]
            }
            ctit <- "black"
            representation <- round(max(dffunc$CumRep) * 100, 2)
            #Wordplot will be plot later
            set.seed(3663)
            maintit <- paste("Most common", analysis, "signatures", sep = " ")
            subtit <- paste("Top 200, or", representation, "% of signatures in sample.")
            par(oma = c(2,7,2,7) + 0.1, xpd = TRUE)
            par(mar = c(4, 7, 4, 7) + 0.1)
            wordcloud(dffunc$Description, dffunc$Freq, normalize_plurals = FALSE, scale = c(3, 0.005), random.order = FALSE)
            title(main = maintit, sub = subtit, col.main = ctit, col.sub = ctit)
        } else {
            flog.info(paste("There are no signatures for", analysis, "available in the featuredose object for this sample."))
        }
    }
}
