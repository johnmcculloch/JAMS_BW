#' calculate_contig_depth_by_length
#'
#' JAMSalpha function
#' @export

calculate_contig_depth_by_length <- function(contigsdata = opt$contigsdata, featuredata = opt$featuredata, max_length_thresh = 2000){

    contigsdata <- contigsdata[order(contigsdata$Length), ]
    totlength <- sum(contigsdata$Length)
    totbases <- sum(contigsdata$NumBases)
    totlength_in_MetaBAT_bins <- sum(contigsdata[which(contigsdata$MetaBATbin != "none"), "Length"])

    contigsdata$Length_CumSum <- cumsum(contigsdata$Length)
    contigsdata$Length_CumPct <- round((contigsdata$Length_CumSum / totlength) , 2)
    contigsdata$NumBases_CumSum <- cumsum(contigsdata$NumBases)
    contigsdata$NumBases_CumPct <- round((contigsdata$NumBases_CumSum / totbases) , 2)

    #Calculate depth rarefaction from min lengths of 500 to max_length_thresh
    Contig_minlength_stats <- data.frame(MinLength = min(contigsdata$Length):max_length_thresh, Sum_NumBases = sapply(min(contigsdata$Length):max_length_thresh, function (x) { sum(contigsdata[which(contigsdata$Length >= x) , "NumBases"]) }), Sum_Length = sapply(min(contigsdata$Length):max_length_thresh, function (x) { sum(contigsdata[which(contigsdata$Length >= x) , "Length"]) }))

    Contig_minlength_stats$Pct_Depth <- round(((Contig_minlength_stats$Sum_NumBases / totbases) * 100), 2)
    Contig_minlength_stats$Pct_Length <- round(((Contig_minlength_stats$Sum_Length / totlength) * 100), 2)

    Contig_minlength_stats$Num_Species <- sapply(min(contigsdata$Length):max_length_thresh, function (x) { length(unique(contigsdata[which(contigsdata$Length >= x) , "Species"])) })
    Contig_minlength_stats$Num_Species_Pct <- round(((Contig_minlength_stats$Num_Species / Contig_minlength_stats$Num_Species[1]) * 100), 2)

    Contig_minlength_stats$Num_LKTs <- sapply(min(contigsdata$Length):max_length_thresh, function (x) { length(unique(contigsdata[which(contigsdata$Length >= x) , "LKT"])) })
    Contig_minlength_stats$Num_LKTs_Pct <- round(((Contig_minlength_stats$Num_LKTs / Contig_minlength_stats$Num_LKTs[1]) * 100), 2)

    calc_pct_domain <- function (contigsdata = NULL, ctf = NULL){
        tmp_contigsdata <- contigsdata[which(contigsdata$Length >= ctf) , c("Domain", "NumBases")]
        tmp_agg <- aggregate(NumBases ~ Domain, data = tmp_contigsdata, FUN = sum)
        tmp_agg$Pct_Relabund <- round(((tmp_agg$NumBases / sum(tmp_agg$NumBases)) * 100), 2)
        tmp_agg$NumBases <- NULL
        rownames(tmp_agg) <- tmp_agg$Domain
        tmp_agg$Domain <- NULL
        tmp_agg <- as.data.frame(t(tmp_agg))
        tmp_agg$MinLength <- ctf

        return(tmp_agg)
    }

    tmp_domain_list <- lapply(Contig_minlength_stats$MinLength, function(x) { calc_pct_domain(contigsdata = contigsdata, ctf = x) } )
    Contig_minlength_stats <- left_join(Contig_minlength_stats, bind_rows(tmp_domain_list), by = "MinLength")

    if ("MetaBATbin" %in% colnames(contigsdata)){
        Contig_minlength_stats$Num_Contigs_in_MetaBAT_bins <- sapply(min(contigsdata$Length):max_length_thresh, function (x) { length(contigsdata[which(contigsdata$Length >= x) , "MetaBATbin"]) })

        calculate_length_pct_in_bins <- function(contigsdata = NULL, ctf = NULL){
            tmp_contigsdata <- contigsdata[which(contigsdata$Length >= x) , c("MetaBATbin", "Length")]
            length_in_bins <- sum(tmp_contigsdata[which(tmp_contigsdata$MetaBATbin != "none"), "Length"])
            length_pct_in_bins <- round(((length_in_bins / totlength_in_MetaBAT_bins) * 100), 2)

            return(length_pct_in_bins)
        }

        Contig_minlength_stats$Contigs_Length_Pct_in_MetaBAT_bins <- sapply(min(contigsdata$Length):max_length_thresh, function (x) { calculate_length_pct_in_bins(contigsdata = contigsdata, ctf = x) })
    }

    #Calculate how many product features are lost by removing short contigs.
    Contig_minlength_stats$Num_genes <- NA
    Contig_minlength_stats$Pct_hypothetical_protein <- NA
    Contig_minlength_stats$Num_genes_axed <- NA
    Contig_minlength_stats$Pct_hypothetical_protein_axed <- NA

    for (ctf in Contig_minlength_stats$MinLength){
        wanted_contigs <- contigsdata[which(contigsdata$Length >= ctf) , "Contig"]
        tmp_products <- featuredata[which(featuredata$Contig %in% wanted_contigs) , "Product"]
        axed_products <- featuredata[which(!(featuredata$Contig %in% wanted_contigs)) , "Product"]
        Contig_minlength_stats[which(Contig_minlength_stats$MinLength == ctf), "Num_genes"] <- length(tmp_products)
        Contig_minlength_stats[which(Contig_minlength_stats$MinLength == ctf), "Pct_hypothetical_protein"] <- round(((length(which(tmp_products == "hypothetical protein")) / length(tmp_products)) * 100), 2)
        Contig_minlength_stats[which(Contig_minlength_stats$MinLength == ctf), "Num_genes_axed"] <- length(axed_products)
        Contig_minlength_stats[which(Contig_minlength_stats$MinLength == ctf), "Pct_hypothetical_protein_axed"] <- round(((length(which(axed_products == "hypothetical protein")) / length(axed_products)) * 100), 2)
    }
    Contig_minlength_stats$Pct_genes_axed <- round(((Contig_minlength_stats$Num_genes_axed / length(featuredata$Product)) * 100), 2)

    return(Contig_minlength_stats)

}
