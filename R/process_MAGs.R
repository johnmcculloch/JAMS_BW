#' process_MAGs(opt = NULL)
#'
#' JAMSalpha function
#' @export

process_MAGs <- function(opt = NULL){

    if ((opt$analysis != "metagenome") || (!("MetaBATbin" %in% colnames(opt$contigsdata)))){

        #Nothing to see here, MetaBAT results are not available.
        opt$MAGsdata <- NULL

    } else {

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
        chechmArgs <- c("predict", "--threads", opt$threads, "--input", binfp, "--output-directory", curr_checkM_output_folder)
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
        opt$MAGsdata <- left_join(checkm_out, MAG_kraken_df, by = "MetaBATbin")

        #Add relative abundance data for downstream quality assessment
        opt$MAGsdata$NumBases <- NA
        for (wantedMAG in MAGnames){
            opt$MAGsdata[which(opt$MAGsdata$MetaBATbin == wantedMAG), "NumBases"] <- sum(subset(opt$contigsdata, MetaBATbin == wantedMAG)[]$NumBases)
        }
        opt$MAGsdata$PPM <- round((opt$MAGsdata$NumBases / sum(opt$contigsdata$NumBases)) * 1E6, 0)

        #Clean up
        unlink(curr_bin_output_folder, recursive = TRUE)
    }

    return(opt)
}