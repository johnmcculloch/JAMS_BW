#' annotate_contigs
#'
#' JAMSalpha function
#' @export

annotate_contigs <- function(opt = NULL){
    #Set working directory to workdir
    setwd(opt$workdir)

    #Write a copy of contigs to be annotated to the system
    flog.info("Writing contigs for annotation.")
    write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = "contigstoannot.fa")

    #annotate with prokka
    flog.info("Annotating contigs with prokka.")
    prokkacmd <- file.path(opt$bindir, "fastannotate_JAMS.sh")
    prokkaargs <- c("-i", "contigstoannot.fa", "-p", opt$prefix, "-t", opt$threads)
    if (opt$analysis %in% c("metagenome")){
        prokkaargs <- c(prokkaargs, "-m")
    }
 
    #Use prokkaJAMS and not prokka if there are over 10,000 contigs.
    #if (length(names(opt$NHcontigs_sequence)) > 10000){
        #flog.info("Annotation will skip tbl2asn step in prokka, as this step is unnecessarily lengthy when there are more than 10,000 contigs.")
        #JAN 2025: skip tbl creation by prokka because it would rarely be used.
        prokkajams <- c("-j")
        prokkaargs <- c(prokkaargs, prokkajams)
    #}

    system(paste(prokkacmd, paste0(prokkaargs, collapse = " ")))
    flog.info("Contigs have been annotated.")
    file.remove("contigstoannot.fa")

    annotationfolder <- paste(opt$prefix, "PROKKA", sep="_")
    flog.info("Gathering predicted proteins.")
    opt$proteome <- read.fasta(file = file.path(annotationfolder, paste(opt$prefix, "faa", sep = ".")), seqtype = "AA")
    flog.info("Gathering predicted genes.")
    opt$genes <- read.fasta(file = file.path(annotationfolder, paste(opt$prefix, "ffn", sep=".")), seqtype = "DNA", forceDNAtolower = FALSE)
    opt$bedfile <- file.path(opt$workdir, annotationfolder, paste(opt$prefix, "bed", sep="."))

    #Load feature information into opt
    flog.info("Making featuredata dataframe from bedfile.")
    opt$featuredata <- make_featuredata_from_bedfile(bedfile = opt$bedfile)
    #Add taxonomy information
    tmpcontigsdata <- opt$contigsdata[, c("Contig", "Taxid", "LKT")]
    opt$featuredata <- left_join(opt$featuredata, tmpcontigsdata, by = "Contig")

    #Bank annotation files to project directory
    if (opt$workdir != opt$sampledir){
        flog.info("Banking annotation to sample directory.")
        copyargs <- c("-R", file.path(opt$workdir, annotationfolder), opt$sampledir)
        system(paste0("cp", copyargs, collapse = " "))
    }

    return(opt)
}
