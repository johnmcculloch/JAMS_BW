#' denoise_contigs
#'
#' JAMSalpha function
#' @export

denoise_contigs <- function(opt = NULL, small_large_length_threshold = 2000) {

    ## Censor small contigs Taxonomy approach
    #Isolate dark matter contigs
    Dark_matter_contigs_df <- opt$contigsdata[which(opt$contigsdata$Taxid == "0"), , drop = FALSE]
    Surviving_contigs_df <- subset((opt$contigsdata[which(opt$contigsdata$Taxid != "0"), , drop = FALSE]), Length >= small_large_length_threshold)
    Unvetted_contigs_df <- subset((opt$contigsdata[which(opt$contigsdata$Taxid != "0"), , drop = FALSE]), Length < small_large_length_threshold)

    ## The genomic "unit" should in theory resolve at the species level, so use that to be conservative. 
    #1) If small contig Species is in large contigs, leave them alone.
    Surviving_contigs_Species <- unique(Surviving_contigs_df$Species)
    Unvetted_contigs_Species <- unique(Unvetted_contigs_df$Species)
 
    #Find out if there are shared species between large and small contigs
    large_small_shared_species <- intersect(Surviving_contigs_Species, Unvetted_contigs_Species)
 
    #Exhume these and deem them safe, if there are any shared.
    if (length(large_small_shared_species) > 0){
        Surviving_contigs_df <- rbind(Surviving_contigs_df, Unvetted_contigs_df[which(Unvetted_contigs_df$Species %in% large_small_shared_species), ])
        Unvetted_contigs_df <- Unvetted_contigs_df[which(!Unvetted_contigs_df$Species %in% large_small_shared_species), , drop = FALSE]
    }

    #2) LKT is only in small contigs, then check genome completeness. If estimated completeness > 10% and relabund > 50 PPM leave alone.
    if (nrow(Unvetted_contigs_df) > 0){
        Unvetted_contigs_est_completeness <- Unvetted_contigs_df %>% group_by(Species) %>% summarise(SumLength = sum(Length), SumBases = sum(NumBases))
        Unvetted_contigs_est_completeness <- as.data.frame(Unvetted_contigs_est_completeness)
        #add expected genome size info
        sizedict <- Unvetted_contigs_df[ , c("Species", "Taxid", "Median_taxid_genome_size", "SD_taxid_genome_size")]
        sizedict <- sizedict[!duplicated(sizedict$Species), , drop = FALSE]
        Unvetted_contigs_est_completeness <- left_join(Unvetted_contigs_est_completeness, sizedict, by = "Species")
        Unvetted_contigs_est_completeness$Est_completeness <- Unvetted_contigs_est_completeness$SumLength / Unvetted_contigs_est_completeness$Median_taxid_genome_size
        Unvetted_contigs_est_completeness$PPM <- round(((Unvetted_contigs_est_completeness$SumBases / sum(opt$contigsdata$NumBases)) * 1E6), 0)
        #Stay of execution for contigs whose estimated completeness > 10% and relabund > 50 PPM
        Species_to_exhume <- Unvetted_contigs_est_completeness[intersect(which(Unvetted_contigs_est_completeness$PPM >= 50), which(Unvetted_contigs_est_completeness$Est_completeness >= 0.1)), "Species"]
        if (length(Species_to_exhume > 0)){
            #Copy to Surviving contigs
            Surviving_contigs_df <- rbind(Surviving_contigs_df, Unvetted_contigs_df[which(Unvetted_contigs_df$Species %in% Species_to_exhume), ])
            #Eliminate from Unvetted_contigs_df
            Unvetted_contigs_df <- Unvetted_contigs_df[which(!(Unvetted_contigs_df$Species %in% Species_to_exhume)), ]
        }
    }

    if (nrow(Unvetted_contigs_df) > 1){
        #Bank original data to opt
        opt$unvetable_contigs_before_censoring <- Unvetted_contigs_df

        #Censor (i.e. default to s__Unclassified) the taxonomies of all unvettable contigs, as they are likely spurious.
        data(JAMStaxtable)
        Unvetted_contigs_df[ , colnames(Unvetted_contigs_df)[colnames(Unvetted_contigs_df) %in% colnames(JAMStaxtable)]] <- JAMStaxtable["0", colnames(Unvetted_contigs_df)[colnames(Unvetted_contigs_df) %in% colnames(JAMStaxtable)]]
    }

    #Reconstruct opt$contigsdata
    opt$contigsdata <- rbind(Surviving_contigs_df, Dark_matter_contigs_df, Unvetted_contigs_df)

    return(opt)

}