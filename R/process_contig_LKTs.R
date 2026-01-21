#' process_MAGs(opt = NULL)
#'
#' JAMSalpha function
#' @export

process_contig_LKTs <- function(opt = NULL){

    if (opt$analysis == "isolate"){

        opt <- evaluate_LKTs(opt = opt, contigsdata_name = "contigsdata", output_list_name = "taxonomic_completeness_and_counts", taxlvls = "LKT")

    } else if (opt$analysis == "metagenome"){
        #Declare final report objective
        taxonomic_space <- "Contig_LKT"
        Bin_ID_cols <- c("Completeness", "Contamination", "Completeness_Model_Used")
        Bin_QC_cols <- c("Quality", "NumBases", "PPM")
        Bin_ass_stats_cols <- c("Genome_Size", "Total_Contigs", "Contig_N50", "Max_Contig_Length", "Total_Coding_Sequences", "GC_Content", "Coding_Density")
        Taxon_info_cols <- c("Taxid", "NCBI_taxonomic_rank", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT", "Gram")
        Taxon_Reference_QC_cols <- c("Num_assemblies_in_taxid", "Num_isolate", "Proportion_of_MAG_in_taxid", "RefScore")
        Final_report_cols <- c(taxonomic_space, Bin_ID_cols, Bin_QC_cols, Bin_ass_stats_cols, Taxon_info_cols, Taxon_Reference_QC_cols)

        opt <- evaluate_LKTs(opt = opt, contigsdata_name = "contigsdata", output_list_name = "taxonomic_completeness_and_counts", taxlvls = "LKT")
        opt$abundances$taxonomic$Contig_LKT <- opt$taxonomic_completeness_and_counts$LKT
        opt$abundances$taxonomic$Contig_LKT$Contig_LKT <- opt$abundances$taxonomic$Contig_LKT$LKT

        #Calculate relative abundance in PPM
        opt$abundances$taxonomic$Contig_LKT$PPM <- round((opt$abundances$taxonomic$Contig_LKT$NumBases / sum(opt$contigsdata$NumBases)) * 1E6, 0)

        #clean up
        opt$taxonomic_completeness_and_counts <- NULL

        #Rearrange columns
        opt$abundances$taxonomic$Contig_LKT <- opt$abundances$taxonomic$Contig_LKT[ , Final_report_cols]
    }

    return(opt)

}