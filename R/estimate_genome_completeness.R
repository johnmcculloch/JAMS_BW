#' estimate_genome_completeness(taxonomic_completeness = NULL, taxon = NULL)
#'
#' JAMSalpha function. Given a JAMSalpha style taxonomic_completeness data frame, will estimate the genome completeness using the "Percent contig sum over expected genome size" method, for entities which cannot be evaluated using checkM (non-bacterial contigs or length too short.)
#' @export


estimate_genome_completeness <- function(taxonomic_completeness = NULL, taxon = NULL){
    curr_ctg_length_sum <- taxonomic_completeness[taxon, "Genome_Size"]
    curr_exp_gen_size <- max(as.numeric(taxonomic_completeness[taxon, "Median_taxid_genome_size"])) #should be the same but just err on the safe side
    curr_pct_of_egs <- round(((curr_ctg_length_sum / curr_exp_gen_size) * 100), 2)

    #If curr_exp_gen_size is NA, like for dark matter, then set completeness to 0
    if (is.na(curr_exp_gen_size)){
        curr_pct_of_egs <- 0
    }

    return(curr_pct_of_egs)
}
