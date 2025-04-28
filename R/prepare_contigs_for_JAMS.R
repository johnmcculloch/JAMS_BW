#' prepare_contigs_for_JAMS
#'
#' JAMSalpha function
#' @export

prepare_contigs_for_JAMS <- function(opt = NULL, fastafile = NULL, contig_minlength = 500){
    #Make sure youre in the right directory
    setwd(opt$workdir)

    #Ensure min contig size is 500 bp.
    flog.info(paste("Ensuring contigs have a minimum size of", contig_minlength, "bp."))
    mycontigs <- read.fasta(file = fastafile, seqtype = "DNA", forceDNAtolower = FALSE)
    mycontigs <- filter_sequence_by_length(sequence = mycontigs, minlength = contig_minlength)

    #Rename as so not to get any unwanted characters
    flog.info("Adjusting sequence headers.")
    mycontigs <- rename_sequences_consecutively(sequence = mycontigs, headerprefix = paste(opt$prefix, "ctg", sep="_"))

    #Write a temporary file to the system so that kraken can classify the contigs
    write.fasta(sequences = mycontigs, names = names(mycontigs), nbchar = 80, file.out = "tempcontigs.fa")

    #Kraken classify contigs
    krakendf <- kraken_classify_taxonomy(opt = opt, fastafile = "tempcontigs.fa", confidence = opt$krakenconfidencescore)
    file.remove("tempcontigs.fa")

    #Filter out any vertebrate DNA if existant
    opt$contigsdata <- subset(krakendf, !(Kingdom %in% c("k__Metazoa", "k__33208_Metazoa")))
    nonhostcontigs <- opt$contigsdata$Sequence
    #Report host contigs in opt, if applicable
    hostcontigs <- krakendf$Sequence[(!(krakendf$Sequence %in% opt$contigsdata$Sequence))]
    if (length(hostcontigs)>0){
        flog.info(paste("A total of", length(hostcontigs), "contigs out of", nrow(krakendf), "were eliminated for being classified as metazoa (host) DNA."))
        #Bank host contigs
        opt$hostcontigsdata <- subset(krakendf, Sequence %in% hostcontigs)
    }

    #Rename sequence column as being Contig
    colnames(opt$contigsdata)[which(colnames(opt$contigsdata) == "Sequence")] <- "Contig"

    #Make a Reference score
    #Each isolate reference is 5 points, each MAG reference is 1 point.
    opt$contigsdata$RefScore <- (as.numeric(opt$contigsdata$Num_isolate) * 5) + as.numeric(opt$contigsdata$Num_MAGs)

    #Bank actual sequences to opt
    opt$NHcontigs_sequence <- filter_sequence_by_name(input_sequences = mycontigs, sequencenames = nonhostcontigs, keep = TRUE)

    #Calculate sequence stats
    seqstats <- get_sequence_stats(sequences = opt$NHcontigs_sequence)
    colnames(seqstats)[which(colnames(seqstats) == "Sequence")] <- "Contig"
    opt$contigsdata <- left_join(opt$contigsdata, seqstats, by = "Contig")

    #Clean_up unwanted columns to keep it simple
    opt$contigsdata <- opt$contigsdata[ , c("Contig", "Taxid", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT", "NCBI_taxonomic_rank", "RefScore", "Median_taxid_genome_size", "SD_taxid_genome_size", "Length", "NumBases", "MetaBATbin")[c("Contig", "Taxid", "Domain", "Species", "LKT", "NCBI_taxonomic_rank", "RefScore", "Median_taxid_genome_size", "SD_taxid_genome_size", "Length", "NumBases", "MetaBATbin") %in% colnames(opt$contigsdata)]]

    #Ensure correct column object class
    numeric_cols <- c("RefScore", "Median_taxid_genome_size", "SD_taxid_genome_size", "Length", "NumBases")[c("RefScore", "Median_taxid_genome_size", "SD_taxid_genome_size", "Length", "NumBases") %in% colnames(opt$contigsdata)]
    for (colm in numeric_cols){
        opt$contigsdata[ , colm] <- as.numeric(opt$contigsdata[ , colm])
    }

    return(opt)
}
