#' evaluate_LKTs(opt = opt, contigsdata = "contigsdata", output_list_name = "taxonomic_completeness_and_counts")
#'
#' JAMSalpha function
#' @export

evaluate_LKTs <- function(opt = opt, contigsdata_name = "contigsdata", output_list_name = "taxonomic_completeness_and_counts", taxlvls = c("LKT", "Species")){

    #Obtain valid taxonomic table for samples
    #ensure compatibility with JAMS v < 2.0.0
    JAMStaxtablefiles <- list.files(path = opt$workingkrakendb, pattern = "\\.rd")
    if (length(JAMStaxtablefiles) > 0){
        if ("JAMStaxtable.rds" %in% JAMStaxtablefiles){
            JAMStaxtablefile <- file.path(opt$workingkrakendb, "JAMStaxtable.rds")
            JAMStaxtable <- readRDS(JAMStaxtablefile)
        } else {
            JAMStaxtablefile <- file.path(opt$workingkrakendb, "JAMStaxtable.rda")
            if (file.exists(JAMStaxtablefile)){
                load(JAMStaxtablefile)
            } else {
                #Fall back on generic taxonomy table and warn user
                flog.info("JAMS taxonomy table not found. Falling back on generic JAMS taxtable.")
                data(JAMStaxtable)
            }
        }
    }

    if (opt$analysis == "metagenome"){
        curr_contigsdata <- opt[[contigsdata_name]]
        taxonomic_completeness_and_counts <- list()
        domains_checkm_will_work_on <- c("d__2157_Archaea", "d__Archaea", "d__2_Bacteria", "d__Bacteria")
        bacterial_contigsdata <- subset(curr_contigsdata, Domain %in% domains_checkm_will_work_on)
        nonbacterial_contigsdata <- subset(curr_contigsdata, !(Domain %in% domains_checkm_will_work_on))
        for (taxlvl in taxlvls){
            flog.info(paste("Evaluating completeness of taxa found at the", taxlvl, "level."))
            taxonomic_completeness <- NULL
            if (nrow(bacterial_contigsdata) > 0) {
                wantedtaxa <- unique(bacterial_contigsdata[ , taxlvl])[!(unique(bacterial_contigsdata[ , taxlvl]) %in% c("none", NA))]
                #Eliminate "s__Unclassified" at Species taxlvl, because at the Species level, completeness of "s__Unclassified" makes no sense because it may belong to several different genera, etc.
                wantedtaxa <- wantedtaxa[wantedtaxa != "s__Unclassified"]
                flog.info(paste("There are", length(wantedtaxa), "bacterial or archaeal taxa at the", taxlvl, "taxonomic level."))
                #Write contigs to temporary bins for completeness evaluation
                curr_bin_output_folder <- file.path(opt$sampledir, "temp_bins")
                #ensure there is no standing bin output folder
                unlink(curr_bin_output_folder, recursive = TRUE)
                dir.create(curr_bin_output_folder, showWarnings = TRUE, recursive = TRUE)
                for (wantedtaxon in wantedtaxa){
                    curr_fn <- paste(wantedtaxon, "fasta", sep = ".")
                    curr_contigs_names <- bacterial_contigsdata[which(bacterial_contigsdata[ , taxlvl] == wantedtaxon), "Contig"]
                    write_contigs_to_system(opt = opt, contig_names = curr_contigs_names, filename = file.path(curr_bin_output_folder, curr_fn))
                }

                #Create a folder for checkM output and run checkm2
                curr_checkM_output_folder <- file.path(curr_bin_output_folder, "CheckM_out")
                #ensure there is no standing checkM output folder
                unlink(curr_checkM_output_folder, recursive = TRUE)
                dir.create(curr_checkM_output_folder, showWarnings = FALSE, recursive = TRUE)
                binfp <- file.path(curr_bin_output_folder, "*.fasta")
                chechmArgs <- c("predict", "--threads", opt$threads, "--input", binfp, "--output-directory", curr_checkM_output_folder)
                flog.info(paste("Evaluating quality of bacterial and archaeal taxa at the", taxlvl, "level with CheckM2"))
                system2('checkm2', args = chechmArgs, stdout = FALSE, stderr = FALSE)

                checkm_out <- fread(file = file.path(curr_checkM_output_folder, "quality_report.tsv"), data.table = FALSE)
                colnames(checkm_out)[which(colnames(checkm_out) == "Name")] <- taxlvl
                checkm_out$Additional_Notes <- NULL
                rownames(checkm_out) <- checkm_out[ , taxlvl]

                #Some checkm outputs may have failed from the contigs being too small or not having ORFs. If there are any, fall back on estimated genome completeness.
                if (nrow(checkm_out) < length(wantedtaxa)){
                    #Find out missing wantedtaxa
                    missingtaxa <- wantedtaxa[!(wantedtaxa %in% checkm_out[ , taxlvl])]
                    suppl_info <- as.data.frame(bacterial_contigsdata[which(bacterial_contigsdata[ , taxlvl] %in%  missingtaxa), ] %>% group_by_at(taxlvl) %>% summarise(Genome_Size = sum(Length), Total_Contigs = length(Length), Max_Contig_Length = max(Length)))
                    rownames(suppl_info) <- suppl_info[ , taxlvl]
                    suppl_df <- as.data.frame(matrix(data = NA, nrow = length(missingtaxa), ncol = ncol(checkm_out)))
                    colnames(suppl_df) <- colnames(checkm_out)
                    suppl_df[ , taxlvl] <- missingtaxa
                    rownames(suppl_df) <- suppl_df[ , taxlvl]
                    for (colm in c("Genome_Size", "Total_Contigs", "Max_Contig_Length")){
                        suppl_df[ , colm] <- suppl_info[rownames(suppl_df), colm]
                    }
                    suppl_df$Completeness_Model_Used <- "Percent contig sum over expected genome size"
                    checkm_out <- rbind(checkm_out, suppl_df)
                }

                taxonomic_completeness <- rbind(taxonomic_completeness, checkm_out)
                #Clean up
                unlink(curr_bin_output_folder, recursive = TRUE)
            } #end conditional for there being bacterial contigs

            if (nrow(nonbacterial_contigsdata) > 0){
                wantedtaxa <- unique(nonbacterial_contigsdata[ , taxlvl])[!unique(nonbacterial_contigsdata[ , taxlvl]) %in% c(NA, "none")]
                #Eliminate "s__Unclassified" at Species taxlvl, because at the Species level, completeness of "s__Unclassified" makes no sense because it may belong to several different genera, etc.
                wantedtaxa <- wantedtaxa[wantedtaxa != "s__Unclassified"]
                flog.info(paste("There are", length(wantedtaxa), "non-bacterial or non-archaeal taxa at the", taxlvl, "taxonomic level."))
                nonbacterial_taxonomic_completeness <- as.data.frame(matrix(data = NA, nrow = length(wantedtaxa), ncol = 13))
                rownames(nonbacterial_taxonomic_completeness) <- wantedtaxa
                colnames(nonbacterial_taxonomic_completeness) <- c(taxlvl, "Completeness", "Contamination", "Completeness_Model_Used", "Translation_Table_Used", "Coding_Density", "Contig_N50", "Average_Gene_Length", "Genome_Size", "GC_Content", "Total_Coding_Sequences", "Total_Contigs", "Max_Contig_Length")
                nonbacterial_taxonomic_completeness[ , taxlvl] <- wantedtaxa
                nonbacterial_taxonomic_completeness[ , "Completeness_Model_Used"] <- "Percent contig sum over expected genome size"
                suppl_info <- as.data.frame(nonbacterial_contigsdata %>% group_by_at(taxlvl) %>% summarise(Genome_Size = sum(Length), Total_Contigs = length(Length), Max_Contig_Length = max(Length)))
                rownames(suppl_info) <- suppl_info[ , taxlvl]

                for (colm in c("Genome_Size", "Total_Contigs", "Max_Contig_Length")){
                    nonbacterial_taxonomic_completeness[ , colm] <- suppl_info[rownames(nonbacterial_taxonomic_completeness), colm]
                }
                taxonomic_completeness <- rbind(taxonomic_completeness, nonbacterial_taxonomic_completeness)
            } #end conditional for there being non-bacterial contigs

            #Add other taxonomic data
            if (taxlvl == "LKT"){
                relevant_field_number <- 5
            } else {
                relevant_field_number <- 3
            }

            taxonomic_completeness$Taxid <- sapply(taxonomic_completeness[ , taxlvl], function (x) { unlist(strsplit(x, split = "_"))[relevant_field_number] })

            Taxinfo <- JAMStaxtable[ , c("Taxid", colnames(JAMStaxtable)[!(colnames(JAMStaxtable) %in% colnames(taxonomic_completeness))])]

            taxonomic_completeness <- left_join(taxonomic_completeness, Taxinfo, by = "Taxid")

            #Add Reference score
            taxonomic_completeness$RefScore <- (as.numeric(taxonomic_completeness$Num_isolate) * 5) + as.numeric(taxonomic_completeness$Num_MAGs)
            rownames(taxonomic_completeness) <- taxonomic_completeness[ , taxlvl]
            #Compute estimated genome completeness for taxa without checkm available
            if ("Percent contig sum over expected genome size" %in% taxonomic_completeness$Completeness_Model_Used){
                estimate_genome_completeness <- function(taxonomic_completeness = NULL, taxon = NULL){
                    curr_ctg_length_sum <- taxonomic_completeness[taxon, "Genome_Size"]
                    curr_exp_gen_size <- max(as.numeric(taxonomic_completeness[taxon, "Median_taxid_genome_size"])) #should be the same but just err on the safe side
                    curr_pct_of_egs <- round(((curr_ctg_length_sum / curr_exp_gen_size) * 100), 2)

                    return(curr_pct_of_egs)
                }

                taxa_to_compute <- rownames(taxonomic_completeness)[which(taxonomic_completeness$Completeness_Model_Used == "Percent contig sum over expected genome size")]
                taxonomic_completeness[taxa_to_compute, "Completeness"] <- sapply(taxa_to_compute, function(x) {estimate_genome_completeness(taxonomic_completeness = taxonomic_completeness, taxon = x)} )
                #Change format to checkM format of completeness capped to 100
                taxonomic_completeness <- reformat_completeness_to_CheckM_style(completeness_df = taxonomic_completeness)

            }

            #Add sequencing depth information
            aggcounts <- curr_contigsdata %>% group_by_at(taxlvl) %>% summarise(NumBases = sum(NumBases))
            taxonomic_completeness <- left_join(taxonomic_completeness, as.data.frame(aggcounts), by = taxlvl)
            rownames(taxonomic_completeness) <- taxonomic_completeness[ , taxlvl]
            taxonomic_completeness <- taxonomic_completeness[order(taxonomic_completeness$Completeness,taxonomic_completeness$Contamination, decreasing = TRUE) , ]

            taxonomic_completeness <- rate_bin_quality(completeness_df = taxonomic_completeness)

            #Stick into list
            taxonomic_completeness_and_counts[[taxlvl]] <- taxonomic_completeness

        } #end loop for taxlvl

    } else if (opt$analysis == "isolate") {

        flog.info("Evaluating completeness of isolate.")
        #Debug later.
    }

    #Bank to opt
    opt[[output_list_name]] <- taxonomic_completeness_and_counts

    return(opt)
}
