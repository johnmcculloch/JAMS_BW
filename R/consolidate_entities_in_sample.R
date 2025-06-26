#' consolidate_entities_in_sample(opt=opt)
#'
#' JAMSalpha function
#' @export


consolidate_entities_in_sample <- function(opt = opt){

    flog.info("Consolidating entities in sample")

    setwd(opt$sampledir)

    opt$contigsdata$PPM <- round(((opt$contigsdata$NumBases / sum(opt$contigsdata$NumBases)) * 1E6), 0)

    Quality_levels_pecking_order <- c("HQ", "MHQ", "MQ", "LQ", "Contaminated")

    #Ensure contigsdata rownames
    rownames(opt$contigsdata) <- opt$contigsdata$Contig
    #Copy congtigsdata to subset as contigs are vetted and accounted for.
    opt$contigsdata_unused <- opt$contigsdata
    #Mark opt$contigsdata with information of how contig was consolidated.
    opt$contigsdata$Consolidation_from <- NA
    opt$contigsdata$ConsolidatedGenomeBin <- NA

    consolidated_entities <- NULL

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

    #Declare final report objective
    taxonomic_space <- "ConsolidatedGenomeBin"
    Bin_ID_cols <- c("Completeness", "Contamination", "Completeness_Model_Used")
    Bin_QC_cols <- c("Quality", "NumBases", "PPM")
    Bin_ass_stats_cols <- c("Genome_Size", "Total_Contigs", "Contig_N50", "Max_Contig_Length", "Total_Coding_Sequences", "GC_Content", "Coding_Density")
    Taxon_info_cols <- c("Taxid", "NCBI_taxonomic_rank", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT", "Gram")
    Taxon_Reference_QC_cols <- c("Num_assemblies_in_taxid", "Num_isolate", "Proportion_of_MAG_in_taxid", "RefScore")
    Final_report_cols <- c(taxonomic_space, Bin_ID_cols, Bin_QC_cols, Bin_ass_stats_cols, Taxon_info_cols, Taxon_Reference_QC_cols)

    #Assess MAGs
    acceptable_MAG_qual_levels <- Quality_levels_pecking_order[1:2]

    if ("MB2bin" %in% names(opt$abundances$taxonomic)){
        flog.info(paste("Finding good quality MetaBAT2 bins at the", paste(acceptable_MAG_qual_levels, collapse = " and "), "level based on modified MIMAG criteria, see https://doi.org/10.1038/s41592-023-01940-w"))

        #Keep bins at HQ and MHQ levels
        if (any(c("HQ", "MHQ") %in% opt$abundances$taxonomic$MB2bin$Quality)){
            MetaBAT_bins_to_keep <- opt$abundances$taxonomic$MB2bin[which(opt$abundances$taxonomic$MB2bin$Quality %in% acceptable_MAG_qual_levels), "MB2bin"]
            surviving_MAGs_df <- opt$abundances$taxonomic$MB2bin[MetaBAT_bins_to_keep, ]
            surviving_MAGs_df$ConsolidatedGenomeBin <- surviving_MAGs_df$MB2bin

            #Mark these contigs as consolidated
            opt$contigsdata[which(opt$contigsdata$MB2bin %in% MetaBAT_bins_to_keep), "Consolidation_from"] <- "MetaBAT2"
            for (cb in rownames(surviving_MAGs_df)){
                opt$contigsdata[which(opt$contigsdata$MB2bin == cb), "ConsolidatedGenomeBin"] <- surviving_MAGs_df[cb, "ConsolidatedGenomeBin"]
            }

            #Consolidate
            surviving_MAGs_df <- surviving_MAGs_df[ , Final_report_cols]
            rownames(surviving_MAGs_df) <- surviving_MAGs_df$ConsolidatedGenomeBin
            consolidated_entities <- rbind(consolidated_entities, surviving_MAGs_df)

            #Remove these contigs from opt$contigsdata_unused
            opt$contigsdata_unused <- opt$contigsdata_unused[opt$contigsdata[which(is.na(opt$contigsdata$ConsolidatedGenomeBin)), "Contig"], , drop = FALSE]

        } else {
            flog.info(paste("No HQ or MHQ level MetaBAT2 bins found."))
        }
        #Report percentage contig depth was consolidated.
        curr_relabund_consolidated <- round((sum(opt$contigsdata[!is.na(opt$contigsdata$ConsolidatedGenomeBin), "PPM"]) / sum(opt$contigsdata$PPM)) * 100, 2)
        flog.info(paste0("At this point, ", curr_relabund_consolidated, "% of contig sequencing depth has been consolidated."))
    }

    for (curr_quality_level in Quality_levels_pecking_order){

        if (nrow(opt$contigsdata_unused) > 0){
            #Continue with consolidation via classical approach.
            curr_taxonomic_completeness_df_name <- "contigs_from_unfit_bins"
            opt <- evaluate_LKTs(opt = opt, taxlvls = "LKT", contigsdata_name = "contigsdata_unused", output_list_name = curr_taxonomic_completeness_df_name)

            flog.info("Estimating genome completeness at higher taxonomic levels for reassigning LKTs. Completeness will subsequently be verified with CheckM2.")

            relevanttaxlvls <- c("IS1", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
            agg_completeness_list <- lapply(relevanttaxlvls, function (x) { glom_completeness(completeness_df = subset(opt[[curr_taxonomic_completeness_df_name]][["LKT"]], LKT != "LKT__Unclassified"), taxlvl = x, High_contamination_threshold = 20) } )
            names(agg_completeness_list) <- relevanttaxlvls

            #Now apply an algorithm of finding best bins at quality level, from LKT to Phylum
            #Make a positional data frame of taxonomic relations of taxa, with their respective completenesses added at each taxonomic level
            colstotag <- c("Completeness", "Contamination", "Quality", "Estimated_num_entities")
            for (taxlvl in names(agg_completeness_list)){
                colnames(agg_completeness_list[[taxlvl]])[colnames(agg_completeness_list[[taxlvl]]) %in% colstotag] <- paste(taxlvl, colnames(agg_completeness_list[[taxlvl]])[colnames(agg_completeness_list[[taxlvl]]) %in% colstotag], sep = "_")
            }
            #Add completeness info at all levels to the master LKT table.
            LKT_est_completeness_df <- opt[[curr_taxonomic_completeness_df_name]][["LKT"]]
            for (taxlvl in relevanttaxlvls){
                LKT_est_completeness_df <- left_join(LKT_est_completeness_df, agg_completeness_list[[taxlvl]], by = taxlvl)
            }

            #Fix LKT data column names
            colnames(LKT_est_completeness_df)[which(colnames(LKT_est_completeness_df) %in% colstotag)] <- paste("LKT", colnames(LKT_est_completeness_df)[which(colnames(LKT_est_completeness_df) %in% colstotag)], sep = "_")

            #Reorganize columns for clarity
            LKT_est_completeness_df <- LKT_est_completeness_df[ , c("LKT", "Completeness_Model_Used", "Taxid", "NCBI_taxonomic_rank", "LKT_Completeness", "LKT_Contamination", "LKT_Quality", "IS1", "IS1_Completeness", "IS1_Contamination", "IS1_Quality", "IS1_Estimated_num_entities", "Species", "Species_Completeness", "Species_Contamination", "Species_Quality", "Species_Estimated_num_entities", "Genus", "Genus_Completeness", "Genus_Contamination", "Genus_Quality", "Genus_Estimated_num_entities", "Family", "Family_Completeness", "Family_Contamination", "Family_Quality", "Family_Estimated_num_entities", "Order", "Order_Completeness", "Order_Contamination", "Order_Quality", "Order_Estimated_num_entities", "Class", "Class_Completeness", "Class_Contamination", "Class_Quality", "Class_Estimated_num_entities", "Phylum", "Phylum_Completeness", "Phylum_Contamination", "Phylum_Quality", "Phylum_Estimated_num_entities", "Kingdom", "Kingdom_Completeness", "Kingdom_Contamination", "Kingdom_Quality", "Kingdom_Estimated_num_entities", "NumBases")]

            LKT_est_completeness_df$PPM <- round((LKT_est_completeness_df$NumBases / sum(opt$contigsdata$NumBases)) * 1E6, 0)
            LKT_est_completeness_df$Taxonomic_level_to_consolidate <- NA
            LKT_est_completeness_df$Best_consolidation_quality_at_taxonomic_level <- NA
            curr_rows_at_quality <- NULL
            curr_unconsolidated_rows <- NULL
            curr_rows_to_consolidate <- NULL

            #Loop round qualities to find the best at the most favourable taxonomic level, except if curr_quality_level is "Contaminated", then consolidate at LKT.
            if (curr_quality_level != "Contaminated"){
                for (quallvl in Quality_levels_pecking_order[Quality_levels_pecking_order != "Contaminated"]){
                    for (taxlvl in c(relevanttaxlvls, "LKT")) {
                        curr_rows_at_quality <- which(LKT_est_completeness_df[ , paste(taxlvl, "Quality", sep = "_")] == quallvl)
                        curr_unconsolidated_rows <- which(is.na(LKT_est_completeness_df[ , "Taxonomic_level_to_consolidate"]))
                        curr_rows_to_consolidate <- intersect(curr_rows_at_quality, curr_unconsolidated_rows)
                        if (length(curr_rows_to_consolidate) > 0){
                            LKT_est_completeness_df[curr_rows_to_consolidate, "Taxonomic_level_to_consolidate"] <- taxlvl
                            LKT_est_completeness_df[curr_rows_to_consolidate, "Best_consolidation_quality_at_taxonomic_level"] <- quallvl
                        }
                        #Reset variables because I am neurotic.
                        curr_rows_at_quality <- NULL
                        curr_unconsolidated_rows <- NULL
                        curr_rows_to_consolidate <- NULL
                    }
                }
            } else {
                LKT_est_completeness_df[ , "Taxonomic_level_to_consolidate"] <- "LKT"
                LKT_est_completeness_df[ , "Best_consolidation_quality_at_taxonomic_level"] <- "Contaminated"
            }

            #Consider that some binning may have higher qualities in subsequent completeness evaluation iterations, so accept any quality which is equal to the current quality level AND ABOVE.
            curr_acceptable_qual_levels <- Quality_levels_pecking_order[1:which(Quality_levels_pecking_order == curr_quality_level)]
            LKT_est_completeness_df <- LKT_est_completeness_df[which(LKT_est_completeness_df$Best_consolidation_quality_at_taxonomic_level %in% curr_acceptable_qual_levels), ]

            if (nrow(LKT_est_completeness_df) > 0){
                #Annotate contigs which survive
                for (curr_agg_tax_lvl in unique(LKT_est_completeness_df$Taxonomic_level_to_consolidate)){
                    tmp_df <- subset(LKT_est_completeness_df, Taxonomic_level_to_consolidate == curr_agg_tax_lvl)[ c("LKT", curr_agg_tax_lvl, paste(curr_agg_tax_lvl, c("Completeness", "Contamination", "Quality"), sep = "_"), "PPM")]
                    #Pick out the contigs from each unique taxon at that taxlvl
                    for (txn in unique(tmp_df[ , curr_agg_tax_lvl])){
                        curr_LKTs_to_mark <- tmp_df[which(tmp_df[ , curr_agg_tax_lvl] == txn), "LKT"]
                        curr_contigs_to_mark <- opt$contigsdata_unused[which(opt$contigsdata_unused$LKT %in% curr_LKTs_to_mark), "Contig"]
                        #Mark these contigs and check if the quality is really adequate with actual Check_M
                        opt$contigsdata_unused[curr_contigs_to_mark, "Consolidation_from"] <- paste("ContigCompleteness", curr_quality_level, sep = "_")
                        opt$contigsdata_unused[curr_contigs_to_mark, "ConsolidatedGenomeBin"] <- paste("cLKT", gsub("^LKT__", "", txn), sep = "__")
                        curr_LKTs_to_mark <- NULL
                        curr_contigs_to_mark <- NULL
                    }
                }

                #Not possible to use the evaluate_LKTs function so recycle the code
                bins_to_eval <- unique(opt$contigsdata_unused$ConsolidatedGenomeBin)[!is.na(unique(opt$contigsdata_unused$ConsolidatedGenomeBin))]
                bins_to_eval_df <- data.frame(ConsolidatedGenomeBin = bins_to_eval)
                bins_to_eval_df$Taxid <- sapply(bins_to_eval_df$ConsolidatedGenomeBin, function (x) { unlist(strsplit(x, split = "_"))[5] } )
                bins_to_eval_df <- left_join(bins_to_eval_df, JAMStaxtable, by = "Taxid")

                #Eliminate bins without taxids (Missing, Unclassified and NA)
                bins_to_eval_df <- bins_to_eval_df[sapply(bins_to_eval_df$Taxid, function (x) { can_be_made_numeric(x) }), ]

                curr_contigsdata <- opt$contigsdata_unused[which(opt$contigsdata_unused$ConsolidatedGenomeBin %in% bins_to_eval_df$ConsolidatedGenomeBin), c("Contig", "ConsolidatedGenomeBin", "Length", "NumBases")]
                curr_contigsdata <- left_join(curr_contigsdata, bins_to_eval_df, by = "ConsolidatedGenomeBin")
                taxonomic_completeness_and_counts <- list()
                domains_checkm_will_work_on <- c("d__2157_Archaea", "d__Archaea", "d__2_Bacteria", "d__Bacteria")
                bacterial_contigsdata <- subset(curr_contigsdata, Domain %in% domains_checkm_will_work_on)
                nonbacterial_contigsdata <- subset(curr_contigsdata, !(Domain %in% domains_checkm_will_work_on))

                taxonomic_completeness <- NULL
                if (nrow(bacterial_contigsdata) > 0) {
                    wantedtaxa <- unique(bacterial_contigsdata[ , "ConsolidatedGenomeBin"])[!(unique(bacterial_contigsdata[ , "ConsolidatedGenomeBin"]) %in% c("none", NA))]
                    #Eliminate "s__Unclassified" at Species taxlvl, because at the Species level, completeness of "s__Unclassified" makes no sense because it may belong to several different genera, etc.
                    wantedtaxa <- wantedtaxa[!(wantedtaxa %in% c("cLKT__d__Unclassified", "cLKT__k__Unclassified", "cLKT__p__Unclassified", "cLKT__c__Unclassified", "cLKT__o__Unclassified", "cLKT__f__Unclassified", "cLKT__g__Unclassified", "cLKT__s__Unclassified", "cLKT__is1__Unclassified", "cLKT__Unclassified"))]

                    flog.info(paste("There are", length(wantedtaxa), "bacterial or archaeal taxa to verify if quality is", curr_quality_level, "or above."))
                    #Write contigs to temporary bins for completeness evaluation
                    curr_bin_output_folder <- file.path(opt$sampledir, "temp_bins")
                    unlink(curr_bin_output_folder, recursive = TRUE)
                    dir.create(curr_bin_output_folder, showWarnings = TRUE, recursive = TRUE)
                    for (wantedtaxon in wantedtaxa){
                        curr_fn <- paste(wantedtaxon, "fasta", sep = ".")
                        curr_contigs_names <- bacterial_contigsdata[which(bacterial_contigsdata[ , "ConsolidatedGenomeBin"] == wantedtaxon), "Contig"]
                        #test if length > 5000, else, don't bother to write to system, because CheckM will fail.
                        if (sum(opt$contigsdata[curr_contigs_names, "Length"]) >= 5000){
                            write_contigs_to_system(opt = opt, contig_names = curr_contigs_names, filename = file.path(curr_bin_output_folder, curr_fn))
                        }
                    }
                    #Create a folder for checkM output and run checkm2
                    curr_checkM_output_folder <- file.path(curr_bin_output_folder, "CheckM_out")
                    #ensure there is no standing checkM output folder
                    unlink(curr_checkM_output_folder, recursive = TRUE)
                    dir.create(curr_checkM_output_folder, showWarnings = FALSE, recursive = TRUE)
                    binfp <- file.path(curr_bin_output_folder, "*.fasta")
                    chechmArgs <- c("predict", "--threads", opt$threads, "--input", binfp, "--output-directory", curr_checkM_output_folder)
                    flog.info(paste("Evaluating quality of bacterial and archaeal taxa with CheckM2"))
                    system2('checkm2', args = chechmArgs, stdout = FALSE, stderr = FALSE)

                    checkm_out <- fread(file = file.path(curr_checkM_output_folder, "quality_report.tsv"), data.table = FALSE)
                    colnames(checkm_out)[which(colnames(checkm_out) == "Name")] <- "ConsolidatedGenomeBin"
                    checkm_out$Additional_Notes <- NULL
                    rownames(checkm_out) <- checkm_out[ , "ConsolidatedGenomeBin"]

                    #Some checkm outputs may have failed from the contigs being too small or not having ORFs. If there are any, fall back on estimated genome completeness.
                    if (nrow(checkm_out) < length(wantedtaxa)){
                        #Find out missing wantedtaxa
                        missingtaxa <- wantedtaxa[!(wantedtaxa %in% checkm_out[ , "ConsolidatedGenomeBin"])]
                        suppl_info <- as.data.frame(bacterial_contigsdata[which(bacterial_contigsdata[ , "ConsolidatedGenomeBin"] %in%  missingtaxa), ] %>% group_by_at("ConsolidatedGenomeBin") %>% summarise(Genome_Size = sum(Length), Total_Contigs = length(Length), Max_Contig_Length = max(Length)))
                        rownames(suppl_info) <- suppl_info[ , "ConsolidatedGenomeBin"]
                        suppl_df <- as.data.frame(matrix(data = NA, nrow = length(missingtaxa), ncol = ncol(checkm_out)))
                        colnames(suppl_df) <- colnames(checkm_out)
                        suppl_df[ , "ConsolidatedGenomeBin"] <- missingtaxa
                        rownames(suppl_df) <- suppl_df[ , "ConsolidatedGenomeBin"]
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
                    wantedtaxa <- unique(nonbacterial_contigsdata[ , "ConsolidatedGenomeBin"])[!unique(nonbacterial_contigsdata[ , "ConsolidatedGenomeBin"]) %in% c(NA, "none")]
                    #Eliminate "s__Unclassified" at Species taxlvl, because at the Species level, completeness of "s__Unclassified" makes no sense because it may belong to several different genera, etc.
                    wantedtaxa <- wantedtaxa[!(wantedtaxa %in% c("cLKT__d__Unclassified", "cLKT__k__Unclassified", "cLKT__p__Unclassified", "cLKT__c__Unclassified", "cLKT__o__Unclassified", "cLKT__f__Unclassified", "cLKT__g__Unclassified", "cLKT__s__Unclassified", "cLKT__is1__Unclassified"))]
                    flog.info(paste("There are", length(wantedtaxa), "non-bacterial or non-archaeal taxa to verify if quality is", curr_quality_level, "or above."))
                    nonbacterial_taxonomic_completeness <- as.data.frame(matrix(data = NA, nrow = length(wantedtaxa), ncol = 13))
                    rownames(nonbacterial_taxonomic_completeness) <- wantedtaxa
                    colnames(nonbacterial_taxonomic_completeness) <- c("ConsolidatedGenomeBin", "Completeness", "Contamination", "Completeness_Model_Used", "Translation_Table_Used", "Coding_Density", "Contig_N50", "Average_Gene_Length", "Genome_Size", "GC_Content", "Total_Coding_Sequences", "Total_Contigs", "Max_Contig_Length")
                    nonbacterial_taxonomic_completeness[ , "ConsolidatedGenomeBin"] <- wantedtaxa
                    nonbacterial_taxonomic_completeness[ , "Completeness_Model_Used"] <- "Percent contig sum over expected genome size"
                    suppl_info <- as.data.frame(nonbacterial_contigsdata %>% group_by(ConsolidatedGenomeBin) %>% summarise(Genome_Size = sum(Length), Total_Contigs = length(Length), Max_Contig_Length = max(Length)))
                    rownames(suppl_info) <- suppl_info[ , "ConsolidatedGenomeBin"]
                    for (colm in c("Genome_Size", "Total_Contigs", "Max_Contig_Length")){
                        nonbacterial_taxonomic_completeness[ , colm] <- suppl_info[rownames(nonbacterial_taxonomic_completeness), colm]
                    }

                    #Add temporary reference genome sizes for completeness evaluation
                    estimate_genome_completeness <- function(taxonomic_completeness = NULL, taxon = NULL){
                        curr_ctg_length_sum <- taxonomic_completeness[taxon, "Genome_Size"]
                        curr_exp_gen_size <- max(as.numeric(taxonomic_completeness[taxon, "Median_taxid_genome_size"])) #should be the same but just err on the safe side
                        curr_pct_of_egs <- round(((curr_ctg_length_sum / curr_exp_gen_size) * 100), 2)

                        return(curr_pct_of_egs)
                    }

                    nonbacterial_taxonomic_completeness <- left_join(nonbacterial_taxonomic_completeness, bins_to_eval_df[ , c("ConsolidatedGenomeBin", "Median_taxid_genome_size")], by = "ConsolidatedGenomeBin") 
                    rownames(nonbacterial_taxonomic_completeness) <- nonbacterial_taxonomic_completeness$ConsolidatedGenomeBin

                    nonbacterial_taxonomic_completeness[ , "Completeness"] <- sapply(nonbacterial_taxonomic_completeness[ , "ConsolidatedGenomeBin"], function(x) {estimate_genome_completeness(taxonomic_completeness = nonbacterial_taxonomic_completeness, taxon = x)} )
                    #Change format to checkM format of completeness capped to 100
                    nonbacterial_taxonomic_completeness <- reformat_completeness_to_CheckM_style(completeness_df = nonbacterial_taxonomic_completeness)
                    nonbacterial_taxonomic_completeness$Median_taxid_genome_size <- NULL

                    taxonomic_completeness <- rbind(taxonomic_completeness, nonbacterial_taxonomic_completeness)
                } #end conditional for there being non-bacterial contigs

                #There may have been no entities evaluated (Taxid Missing, Unclassified, etc.)
                if (!(is.null(taxonomic_completeness) || (nrow(taxonomic_completeness) == 0))){
                    taxonomic_completeness <- rate_bin_quality(completeness_df = taxonomic_completeness)

                    #Consolidate only the bins which are at the current desired level.

                    #Bank these to consolidated_entities
                    curr_consolidated_entities <- taxonomic_completeness[which(taxonomic_completeness$Quality %in% curr_acceptable_qual_levels), ]
                    curr_consolidated_entities <- left_join(curr_consolidated_entities, bins_to_eval_df, by = "ConsolidatedGenomeBin")

                    #Censor bins below acceptable quality in opt$contigsdata_unused and consolidate OK bins to opt$contigsdata
                    bins2censor <- taxonomic_completeness[which(!(taxonomic_completeness$Quality %in% curr_acceptable_qual_levels)), "ConsolidatedGenomeBin"]
                    #Reset back to NA
                    opt$contigsdata_unused[which(opt$contigsdata_unused$ConsolidatedGenomeBin %in% bins2censor), "ConsolidatedGenomeBin"] <- NA
                    opt$contigsdata_unused[which(opt$contigsdata_unused$ConsolidatedGenomeBin %in% bins2censor), "Consolidation_from"] <- NA

                    #Add PPM and RefScore
                    curr_contigs_to_consolidate <- opt$contigsdata_unused[which(opt$contigsdata_unused$ConsolidatedGenomeBin %in% curr_consolidated_entities$ConsolidatedGenomeBin), "Contig"]

                    curr_agg_df <- as.data.frame(opt$contigsdata_unused[curr_contigs_to_consolidate, c("Contig", "ConsolidatedGenomeBin", "NumBases")] %>% group_by(ConsolidatedGenomeBin) %>% summarise(NumBases = sum(NumBases)))
                    curr_agg_df$PPM <- round((curr_agg_df$NumBases / sum(opt$contigsdata$NumBases)) * 1E6, 0)

                    curr_consolidated_entities <- left_join(curr_consolidated_entities, curr_agg_df, by = "ConsolidatedGenomeBin")
                    curr_agg_df <- NULL #Neurotic me
                    curr_consolidated_entities$RefScore <- (as.numeric(curr_consolidated_entities$Num_isolate) * 5) + as.numeric(curr_consolidated_entities$Num_MAGs)
                    rownames(curr_consolidated_entities) <- curr_consolidated_entities[ , "ConsolidatedGenomeBin"]

                    curr_consolidated_entities <- curr_consolidated_entities[ , Final_report_cols]

                    #Bank to consolidated_entities
                    consolidated_entities <- rbind(consolidated_entities, curr_consolidated_entities)

                    #Annotate opt$contigsdata
                    curr_contigs_to_consolidate <- NULL
                    for (csb in unique(curr_consolidated_entities$ConsolidatedGenomeBin)){
                        curr_contigs_to_consolidate <- opt$contigsdata_unused[which(opt$contigsdata_unused$ConsolidatedGenomeBin == csb), "Contig"]
                        opt$contigsdata[curr_contigs_to_consolidate, "ConsolidatedGenomeBin"] <- csb
                        opt$contigsdata[curr_contigs_to_consolidate, "Consolidation_from"] <- "Kraken2CheckM2"
                    }
                    curr_contigs_to_consolidate <- NULL

                    #Remove consolidated contigs from opt$contigsdata_unused
                    opt$contigsdata_unused <- opt$contigsdata_unused[which(is.na(opt$contigsdata_unused$ConsolidatedGenomeBin)), ]
                } else {
                    ##Reset back to NA
                    opt$contigsdata_unused$ConsolidatedGenomeBin <- NA
                    opt$contigsdata_unused$Consolidation_from <- NA
                }

            } else {

                flog.info(paste("No higher level bins found at the", curr_quality_level, "quality level."))

            }#end conditional that there were any aggregated LKTs at that particular quality cutoff

        }#end conditional for there being any unused contigs. They may have been used up in higher quality levels.

        #Report percentage contig depth was consolidated.
        curr_relabund_consolidated <- round((sum(opt$contigsdata[!is.na(opt$contigsdata$ConsolidatedGenomeBin), "PPM"]) / sum(opt$contigsdata$PPM)) * 100, 2)
        flog.info(paste0("At this point, ", curr_relabund_consolidated, "% of contig sequencing depth has been consolidated."))

    } #End loop for Quality levels
 
    #Lastly, consolidate dark matter
    opt$leftover <- opt$contigsdata[which(is.na(opt$contigsdata$ConsolidatedGenomeBin)), , drop = FALSE]
    if (nrow(opt$leftover) > 0){
        flog.info("Consolidating genomic obscure or dark matter.")
        opt <- evaluate_LKTs(opt = opt, taxlvls = "LKT", contigsdata_name = "leftover", output_list_name = "leftover_est_completeness")
        dark_matter_taxonomic_completeness <- opt$leftover_est_completeness$LKT
        dark_matter_taxonomic_completeness$ConsolidatedGenomeBin <- paste0("c", dark_matter_taxonomic_completeness$LKT)
        dark_matter_taxonomic_completeness$PPM <- round((dark_matter_taxonomic_completeness$NumBases / sum(opt$contigsdata$NumBases)) * 1E6, 0)
        dark_matter_taxonomic_completeness <- dark_matter_taxonomic_completeness[ , colnames(consolidated_entities)]
        rownames(dark_matter_taxonomic_completeness) <- dark_matter_taxonomic_completeness$ConsolidatedGenomeBin
        #Annotate opt$contigsdata
        curr_contigs_to_consolidate <- NULL
        for (csbLKT in unique(dark_matter_taxonomic_completeness$LKT)){
            curr_contigs_to_consolidate <- opt$leftover[which(opt$leftover$LKT == csbLKT), "Contig"]
            csb <- dark_matter_taxonomic_completeness[which(dark_matter_taxonomic_completeness$LKT == csbLKT), "ConsolidatedGenomeBin"]
            opt$contigsdata[curr_contigs_to_consolidate, "ConsolidatedGenomeBin"] <- csb
            opt$contigsdata[curr_contigs_to_consolidate, "Consolidation_from"] <- "Kraken2CheckM2"
        }
        curr_contigs_to_consolidate <- NULL
        #Bank to opt$abundances$taxonomic
        opt$abundances$taxonomic$ConsolidatedGenomeBin <- rbind(consolidated_entities, dark_matter_taxonomic_completeness)

        #Report percentage contig depth was consolidated.
        curr_relabund_consolidated <- round((sum(opt$contigsdata[!is.na(opt$contigsdata$ConsolidatedGenomeBin), "PPM"]) / sum(opt$contigsdata$PPM)) * 100, 2)
        flog.info(paste0("At this point, ", curr_relabund_consolidated, "% of contig sequencing depth has been consolidated."))
    } 

    opt$abundances$taxonomic$ConsolidatedGenomeBin <- opt$abundances$taxonomic$ConsolidatedGenomeBin[ , Final_report_cols]

    #Add consolidation information to featuredata
    opt$featuredata <- left_join(opt$featuredata, opt$contigsdata[ , c("Contig", "ConsolidatedGenomeBin")], by = "Contig")
    opt$contigsdata$PPM <- NULL

    return(opt)
}
