#' sort_contigs_into_bins(opt = NULL)
#'
#' JAMSalpha function
#' @export

sort_contigs_into_bins <- function(opt = NULL, sortedbamfile = NULL){

    #Do MetaBAT if applicable
    if (opt$bin_contigs_with_metabat){
        setwd(opt$covdir)
        if (is.null(sortedbamfile)){
            sortedbamfile <- paste(paste(opt$prefix, "contigs", sep="_"), "sorted", "bam", sep=".")
        }

        flog.info("Binning contigs from assembly with metaBAT2")
        system(paste("jgi_summarize_bam_contig_depths --outputDepth depth.txt", sortedbamfile))

        #Bank contig depths for later if needed
        opt$bam_contig_depths <- fread("depth.txt", data.table = FALSE)
        rownames(opt$bam_contig_depths) <- opt$bam_contig_depths$contigName

        #Run metabat2
        system("metabat2 -i contigscov.fa -a depth.txt -o bins_dir/MAGbin")
        flog.info("Integrating binning information into JAMS")
        MAGbin_files <- list.files("bins_dir", pattern = "fa$")

        get_contigs_in_bin <- function(bin_file = NULL){
            bin_contigs <- NULL
            if (file.exists(file.path("bins_dir", bin_file))){
                bin_contigs <- system2("grep", args = c("\'^>\'", file.path("bins_dir", bin_file)), stdout = TRUE, stderr = FALSE)
                bin_contigs <- gsub("^>", "", bin_contigs)
            }
            if (length(bin_contigs) < 1){
                flog.warn(paste("MetaBAT2 bin file", bin_file, "not found, or empty."))
            }

            return(bin_contigs)
        }

        bin_contigs_list <- lapply(MAGbin_files, function(x) { get_contigs_in_bin(bin_file = x) } )
        names(bin_contigs_list) <- MAGbin_files

        opt$contigsdata$MetaBATbin <- 0

        for (mbl in 1:length(bin_contigs_list)){
            bin_number <- as.numeric(unlist(strsplit(names(bin_contigs_list)[mbl], split = "\\."))[2])
            curr_contigs <- bin_contigs_list[[mbl]]
            opt$contigsdata[which(opt$contigsdata$Contig %in% curr_contigs), "MetaBATbin"] <- bin_number
        }

        #Standardise number of charaters
        opt$contigsdata$MetaBATbin <- formatC(opt$contigsdata$MetaBATbin, width = nchar(max(opt$contigsdata$MetaBATbin)), flag = "0")
        opt$contigsdata$MetaBATbin[which(as.numeric(opt$contigsdata$MetaBATbin) == 0)] <- "none"
        opt$contigsdata$MetaBATbin[which(opt$contigsdata$MetaBATbin != "none")] <- paste(opt$prefix, "MAGbin", opt$contigsdata$MetaBATbin[which(opt$contigsdata$MetaBATbin != "none")], sep = "_")

    }

    return(opt)
}