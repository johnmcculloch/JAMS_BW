#' evaluate_LKTs(opt = opt, contigsdata = "contigsdata", output_list_name = "taxonomic_completeness_and_counts")
#'
#' JAMSalpha function
#' @export

evaluate_LKTs <- function(opt = opt, contigsdata_name = "contigsdata", output_list_name = "taxonomic_completeness_and_counts", taxlvls = c("LKT", "Species")){

    setwd(opt$sampledir)

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

                #Check whether there is at least a single taxon whose sum length > 5000 bp. Else, do not use CheckM because there will be no output.
                SumLength_tally <- bacterial_contigsdata %>% filter(.data[[taxlvl]] %in% wantedtaxa) %>% group_by(.data[[taxlvl]]) %>% summarise(SumLength = sum(Length), .groups = "drop")
                if (length(which(SumLength_tally$SumLength >= 5000)) > 0){
                    #Write contigs to temporary bins for completeness evaluation
                    curr_bin_output_folder <- file.path(opt$sampledir, "temp_bins")
                    #ensure there is no standing bin output folder
                    unlink(curr_bin_output_folder, recursive = TRUE)
                    dir.create(curr_bin_output_folder, showWarnings = TRUE, recursive = TRUE)
                    for (wantedtaxon in wantedtaxa){
                        curr_fn <- paste(wantedtaxon, "fasta", sep = ".")
                        curr_contigs_names <- bacterial_contigsdata[which(bacterial_contigsdata[ , taxlvl] == wantedtaxon), "Contig"]
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
                    appropriatenumcores <- max(2, (opt$threads - 2))
                    checkmArgs <- c("predict", "--database_path", opt$CheckMdb, "--threads", appropriatenumcores, "--input", binfp, "--output-directory", curr_checkM_output_folder)
                    flog.info(paste("Evaluating quality of bacterial and archaeal taxa at the", taxlvl, "level with CheckM2"))
                    system2('checkm2', args = checkmArgs, stdout = FALSE, stderr = FALSE)

                    #Check an output exists before trying to read. Some garbage input may be too low even for a single assessement on a quality_report
                    if (file.exists(file.path(curr_checkM_output_folder, "quality_report.tsv"))){
                        checkm_out <- fread(file = file.path(curr_checkM_output_folder, "quality_report.tsv"), data.table = FALSE)
                        colnames(checkm_out)[which(colnames(checkm_out) == "Name")] <- taxlvl
                        checkm_out$Additional_Notes <- NULL
                        rownames(checkm_out) <- checkm_out[ , taxlvl]
                    } else {
                        checkm_out <- NULL
                    }
                }

                #Some checkm outputs may have failed from the contigs being too small or not having ORFs. If there are any, fall back on estimated genome completeness.
                if ((is.null(checkm_out) || (nrow(checkm_out) < length(wantedtaxa)))){
                    #Find out missing wantedtaxa
                    missingtaxa <- wantedtaxa[!(wantedtaxa %in% checkm_out[ , taxlvl])]
                    suppl_info <- as.data.frame(bacterial_contigsdata[which(bacterial_contigsdata[ , taxlvl] %in%  missingtaxa), ] %>% group_by_at(taxlvl) %>% summarise(Genome_Size = sum(Length), Total_Contigs = length(Length), Max_Contig_Length = max(Length)))
                    rownames(suppl_info) <- suppl_info[ , taxlvl]
                    if (!is.null(checkm_out)){
                        suppl_df <- as.data.frame(matrix(data = NA, nrow = length(missingtaxa), ncol = ncol(checkm_out)))
                        colnames(suppl_df) <- colnames(checkm_out)
                    } else {
                        suppl_df <- as.data.frame(matrix(data = NA, nrow = length(wantedtaxa), ncol = 13))
                        colnames(suppl_df) <- c(taxlvl, "Completeness", "Contamination", "Completeness_Model_Used", "Translation_Table_Used", "Coding_Density", "Contig_N50", "Average_Gene_Length", "Genome_Size", "GC_Content", "Total_Coding_Sequences", "Total_Contigs", "Max_Contig_Length")
                    }
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

        #Bank to opt
        opt[[output_list_name]] <- taxonomic_completeness_and_counts

    } else if (opt$analysis == "isolate") {

        flog.info("Evaluating completeness of isolate.")
        #Debug later.
        curr_contigsdata <- opt[[contigsdata_name]]
        taxonomic_completeness_and_counts <- list()
        domains_checkm_will_work_on <- c("d__2157_Archaea", "d__Archaea", "d__2_Bacteria", "d__Bacteria")
        bacterial_contigsdata <- subset(curr_contigsdata, Domain %in% domains_checkm_will_work_on)
        nonbacterial_contigsdata <- subset(curr_contigsdata, !(Domain %in% domains_checkm_will_work_on))

        #Warn if there are any non-bacterial contigs
        if (nrow(nonbacterial_contigsdata > 0)){
            flog.warn(paste0("WARNING: ", nrow(nonbacterial_contigsdata), "/", nrow(curr_contigsdata), " contigs were found to be non-bacterial or non archaeal. Check your bacterial/archaeal isolate for contamination."))
        }

        curr_bin_output_folder <- file.path(opt$sampledir, "temp_bins")
        #ensure there is no standing bin output folder
        unlink(curr_bin_output_folder, recursive = TRUE)
        dir.create(curr_bin_output_folder, showWarnings = TRUE, recursive = TRUE)
        curr_fn <- paste(opt$prefix, "fasta", sep = ".")
        curr_contigs_names <- bacterial_contigsdata[ , "Contig"]
        #test if length > 5000, else, don't bother to write to system, because CheckM will fail.
        if (sum(opt$contigsdata[curr_contigs_names, "Length"]) >= 5000){
            write_contigs_to_system(opt = opt, contig_names = curr_contigs_names, filename = file.path(curr_bin_output_folder, curr_fn))
        }

        #Create a folder for checkM output and run checkm2
        curr_checkM_output_folder <- file.path(curr_bin_output_folder, "CheckM_out")
        #ensure there is no standing checkM output folder
        unlink(curr_checkM_output_folder, recursive = TRUE)
        dir.create(curr_checkM_output_folder, showWarnings = FALSE, recursive = TRUE)
        binfp <- file.path(curr_bin_output_folder, "*.fasta")
        appropriatenumcores <- max(2, (opt$threads - 2))
        checkmArgs <- c("predict", "--database_path", opt$CheckMdb, "--threads", appropriatenumcores, "--input", binfp, "--output-directory", curr_checkM_output_folder)
        flog.info(paste("Evaluating quality of bacterial or archaeal contigs with CheckM2"))
        system2('checkm2', args = checkmArgs, stdout = FALSE, stderr = FALSE)

        checkm_out <- fread(file = file.path(curr_checkM_output_folder, "quality_report.tsv"), data.table = FALSE)
        colnames(checkm_out)[which(colnames(checkm_out) == "Name")] <- "Sequence"
        checkm_out$Additional_Notes <- NULL
        rownames(checkm_out) <- checkm_out[ , "Sequence"]

        #Classify isolate taxonomically as a single entity
        #concatenate multifasta file to reclassify by k-mer spectrum
        setwd(curr_bin_output_folder)
        system(paste(c("grep -v '^>'", file.path(curr_bin_output_folder, curr_fn), ">", file.path(curr_bin_output_folder, "tmp")), collapse = " "))
        cat(paste0(">", opt$prefix), sep = "\n", file = file.path(curr_bin_output_folder, "header"))
        system(paste("cat", file.path(curr_bin_output_folder, "header"), file.path(curr_bin_output_folder, "tmp") ,">>", file.path(curr_bin_output_folder, "Isolate_fullgenome.fasta")))

        Isolate_kraken_df <- kraken_classify_taxonomy(opt = opt, fastafile = file.path(curr_bin_output_folder, "Isolate_fullgenome.fasta"), confidence = 0)
        opt$abundances$taxonomic$ConsolidatedGenomeBin <- left_join(checkm_out, Isolate_kraken_df, by = "Sequence")

        #Add relative abundance data for downstream quality assessment
        opt$abundances$taxonomic$ConsolidatedGenomeBin$NumBases <- sum(bacterial_contigsdata$NumBases)
        opt$abundances$taxonomic$ConsolidatedGenomeBin$PPM <- round((opt$abundances$taxonomic$ConsolidatedGenomeBin$NumBases / sum(opt$contigsdata$NumBases)) * 1E6, 0)

        #Add Reference score
        opt$abundances$taxonomic$ConsolidatedGenomeBin$RefScore <- (as.numeric(opt$abundances$taxonomic$ConsolidatedGenomeBin$Num_isolate) * 5) + as.numeric(opt$abundances$taxonomic$ConsolidatedGenomeBin$Num_MAGs)

        #Add CGB information to opt$featuredata. This being an isolate there should be a single value for LKT
        opt$featuredata$ConsolidatedGenomeBin <- Isolate_kraken_df$LKT
        #Rename tag from LKT to CGB to reflect single-best LKT for entire genome.
        opt$featuredata$ConsolidatedGenomeBin <- gsub("^LKT__", "CGB__", opt$featuredata$ConsolidatedGenomeBin)

        #Add quality rating
        opt$abundances$taxonomic$ConsolidatedGenomeBin <- rate_bin_quality(completeness_df = opt$abundances$taxonomic$ConsolidatedGenomeBin)

        #Clean up
        unlink(curr_bin_output_folder, recursive = TRUE)
    }

    return(opt)
}
