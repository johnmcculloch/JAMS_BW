#' prepare_contigs_for_JAMS
#'
#' JAMSalpha function
#' @export

prepare_contigs_for_JAMS <- function(opt = NULL, fastafile = NULL){
    #Make sure youre in the right directory
    setwd(opt$workdir)

    #Ensure min contig size is 500 bp.
    flog.info("Ensuring contigs have a minimum size of 500 bp.")
    mycontigs <- read.fasta(file = fastafile, seqtype = "DNA", forceDNAtolower = FALSE)
    mycontigs <- filter_sequence_by_length(sequence = mycontigs, minlength = 500)

    #Rename as so not to get any unwanted characters
    flog.info("Adjusting sequence headers.")
    mycontigs <- rename_sequences_consecutively(sequence = mycontigs, headerprefix = paste(opt$prefix, "ctg", sep="_"))

    #Write a temporary file to the system so that kraken can classify the contigs
    write.fasta(sequences = mycontigs, names = names(mycontigs), nbchar = 80, file.out = "tempcontigs.fa")

    #Kraken classify contigs
    krakendf <- kraken_classify_taxonomy(opt = opt, fastafile = "tempcontigs.fa", confidence = opt$krakenconfidencescore)
    file.remove("tempcontigs.fa")

    #Filter out any vertebrate DNA if existant
    opt$contigsdata <- subset(krakendf, Phylum != "p__Chordata")
    nonhostcontigs <- opt$contigsdata$Sequence
    #Report host contigs in opt, if applicable
    hostcontigs <- krakendf$Sequence[(!(krakendf$Sequence %in% opt$contigsdata$Sequence))]
    if(length(hostcontigs)>0){
        flog.info(paste("A total of",length(hostcontigs), "contigs out of", nrow(krakendf), "were eliminated for being classified as vertebrate (host) DNA."))
        opt$hostcontigs <- hostcontigs
    }

    #Rename sequence column as being Contig
    colnames(opt$contigsdata)[which(colnames(opt$contigsdata) == "Sequence")] <- "Contig"

    opt$NHcontigs_sequence <- filter_sequence_by_name(input_sequences = mycontigs, sequencenames = nonhostcontigs, keep = TRUE)

    return(opt)
}
