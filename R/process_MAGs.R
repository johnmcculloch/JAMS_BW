#' process_MAGs(opt = NULL)
#'
#' JAMSalpha function
#' @export

process_MAGs <- function(opt = NULL){

    setwd(opt$sampledir)

    if ((opt$analysis != "metagenome") || (!("MetaBATbin" %in% colnames(opt$contigsdata)))){

        #Nothing to see here, MetaBAT results are not available.
        opt$abundances$taxonomic$MetaBATbin <- NULL

    } else {

        #Declare final report objective
        taxonomic_space <- "MB2bin"
        Bin_ID_cols <- c("Completeness", "Contamination", "Completeness_Model_Used")
        Bin_QC_cols <- c("Quality", "NumBases", "PPM")
        Bin_ass_stats_cols <- c("Genome_Size", "Total_Contigs", "Contig_N50", "Max_Contig_Length", "Total_Coding_Sequences", "GC_Content", "Coding_Density")
        Taxon_info_cols <- c("Taxid", "NCBI_taxonomic_rank", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT", "Gram")
        Taxon_Reference_QC_cols <- c("Num_assemblies_in_taxid", "Num_isolate", "Proportion_of_MAG_in_taxid", "RefScore")
        Final_report_cols <- c(taxonomic_space, Bin_ID_cols, Bin_QC_cols, Bin_ass_stats_cols, Taxon_info_cols, Taxon_Reference_QC_cols)

        #Process MAGs
        MAGnames <- unique(opt$contigsdata$MetaBATbin)[unique(opt$contigsdata$MetaBATbin) != "none"]
        flog.info(paste("There are", length(MAGnames), "metagenome assembled genome (MAG) bins to evaluate."))

        #Write contigs to temporary bins for completeness evaluation
        curr_bin_output_folder <- file.path(opt$sampledir, "temp_bins")
        dir.create(curr_bin_output_folder, showWarnings = TRUE, recursive = TRUE)
        for (wantedMAG in MAGnames){
            curr_fn <- paste(wantedMAG, "fasta", sep = ".")
            curr_contigs_names <- subset(opt$contigsdata, MetaBATbin == wantedMAG)[]$Contig
            write_contigs_to_system(opt = opt, contig_names = curr_contigs_names, filename = file.path(curr_bin_output_folder, curr_fn))
        }

        #Create a folder for checkM output and run checkm2
        curr_checkM_output_folder <- file.path(curr_bin_output_folder, "CheckM_out")
        dir.create(curr_checkM_output_folder, showWarnings = TRUE, recursive = TRUE)
        binfp <- file.path(curr_bin_output_folder, "*.fasta")
        appropriatenumcores <- max(2, (opt$threads - 2))
        chechmArgs <- c("predict", "--threads", appropriatenumcores, "--input", binfp, "--output-directory", curr_checkM_output_folder)
        flog.info("Evaluating quality of MAGs with CheckM2")
        system2('checkm2', args = chechmArgs, stdout = FALSE, stderr = FALSE)
        checkm_out <- fread(file = file.path(curr_checkM_output_folder, "quality_report.tsv"), data.table = FALSE)

        #Classify MAGs taxonomically as a single entity
        #concatenate all multifasta bins to reclassify by k-mer spectrum
        for (wantedMAG in MAGnames){
            curr_fn <- paste(wantedMAG, "fasta", sep = ".")
            system(paste(c("grep -v '^>'", file.path(curr_bin_output_folder, curr_fn), ">", file.path(curr_bin_output_folder, "tmp")), collapse = " "))
            cat(paste0(">", wantedMAG), sep = "\n", file = file.path(curr_bin_output_folder, "header"))
            system(paste("cat", file.path(curr_bin_output_folder, "header"), file.path(curr_bin_output_folder, "tmp") ,">>", file.path(curr_bin_output_folder, "All_concat_MAGs.fasta")))
        }
        MAG_kraken_df <- kraken_classify_taxonomy(opt = opt, fastafile = file.path(curr_bin_output_folder, "All_concat_MAGs.fasta"), confidence = 0)
        colnames(MAG_kraken_df)[which(colnames(MAG_kraken_df) == "Sequence")] <- "MetaBATbin"
        colnames(checkm_out)[which(colnames(checkm_out) == "Name")] <- "MetaBATbin"
        opt$abundances$taxonomic$MB2bin <- left_join(checkm_out, MAG_kraken_df, by = "MetaBATbin")

        #Add relative abundance data for downstream quality assessment
        opt$abundances$taxonomic$MB2bin$NumBases <- NA
        for (wantedMAG in MAGnames){
            opt$abundances$taxonomic$MB2bin[which(opt$abundances$taxonomic$MB2bin$MetaBATbin == wantedMAG), "NumBases"] <- sum(subset(opt$contigsdata, MetaBATbin == wantedMAG)[]$NumBases)
        }
        opt$abundances$taxonomic$MB2bin$PPM <- round((opt$abundances$taxonomic$MB2bin$NumBases / sum(opt$contigsdata$NumBases)) * 1E6, 0)

        #Add Reference score
        opt$abundances$taxonomic$MB2bin$RefScore <- (as.numeric(opt$abundances$taxonomic$MB2bin$Num_isolate) * 5) + as.numeric(opt$abundances$taxonomic$MB2bin$Num_MAGs)

        #Rename bins to include LKT in name
        opt$abundances$taxonomic$MB2bin$BinID <- unname(sapply(opt$abundances$taxonomic$MB2bin$MetaBATbin, function (x) { tail(unlist(strsplit(x, split = "_")), 1) }))
        opt$abundances$taxonomic$MB2bin$MB2bin <- NA
        opt$contigsdata$MB2bin <- NA

        #Construct genome bin name
        for (rn in 1:nrow(opt$abundances$taxonomic$MB2bin)){
            curr_MB2name <- paste("MB2", opt$abundances$taxonomic$MB2bin[rn, "BinID"], gsub("^LKT__", "", opt$abundances$taxonomic$MB2bin[rn, "LKT"]), sep = "__")
            #Final bin name to abundance table
            opt$abundances$taxonomic$MB2bin[rn, "MB2bin"] <- curr_MB2name
            #Mirror this information on opt$contigsdata, for use by harvesting functions
            opt$contigsdata[which(opt$contigsdata$MetaBATbin == opt$abundances$taxonomic$MB2bin[rn, "MetaBATbin"]), "MB2bin"] <- curr_MB2name
        }
        rownames(opt$abundances$taxonomic$MB2bin) <- opt$abundances$taxonomic$MB2bin$MB2bin

        #Add "none" to no bin contigs in opt$contigsdata
        opt$contigsdata$MB2bin[which(is.na(opt$contigsdata$MB2bin))] <- "none"

        #Rank bins
        opt$abundances$taxonomic$MB2bin <- rate_bin_quality(completeness_df = opt$abundances$taxonomic$MB2bin)
        #Reorder columns for standardised output
        opt$abundances$taxonomic$MB2bin <- opt$abundances$taxonomic$MB2bin[ , Final_report_cols[Final_report_cols %in% colnames(opt$abundances$taxonomic$MB2bin)]]

        #Clean up
        unlink(curr_bin_output_folder, recursive = TRUE)
    }

    #Add MB2 information to featuredata
    opt$featuredata <- left_join(opt$featuredata, opt$contigsdata[ , c("Contig", "MB2bin")], by = "Contig")
    rownames(opt$featuredata) <- opt$featuredata$Feature

    #Back to where we were
    setwd(opt$sampledir)

    return(opt)
}