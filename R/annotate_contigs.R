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
    #Use prokkaJAMS and not prokka if there are over 10,000 contigs.
    if (length(names(opt$NHcontigs_sequence)) > 10000){
        flog.info("Annotation will skip tbl2asn step in prokka, as this step is unnecessarily lengthy when there are more than 10,000 contigs.")
        prokkajams <- c("-j")
        prokkaargs <- c(prokkaargs, prokkajams)
    }

    #system2(prokkacmd, args = prokkaargs, stdout = TRUE, stderr = TRUE)
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
    feature2contig <- read.table(file = file.path(annotationfolder, "feature2contig.map"), sep = "\t", header = FALSE, skipNul = FALSE, fill = TRUE, colClasses = c("character", "numeric","character", "character"))
    colnames(feature2contig) <- c("Feature", "LengthDNA", "FeatType", "Contig")
    feature2product <- read.table(file = file.path(annotationfolder, "feature2product.map"), sep = "\t", header = FALSE, skipNul = FALSE, fill = TRUE, colClasses = c("character", "character"))
    colnames(feature2product) <- c("Feature", "Product")

    #Consider features and not genes
    feature2contig <- subset(feature2contig, FeatType %in% c("CDS", "tRNA",  "rRNA",  "tmRNA"))
    opt$featuredata <- left_join(feature2contig, feature2product)
    opt$featuredata$Product[is.na(opt$featuredata$Product)] <- "none"

    #Add enzymes if there are any (small contigs may not have any)
    CDS2EC <- read.table(file = file.path(annotationfolder, "CDS2EC.map"), sep = "\t", header = FALSE, skipNul = FALSE, fill = TRUE, colClasses = c("character", "character"))
    if (nrow(CDS2EC) > 0){
        colnames(CDS2EC) <- c("Feature", "ECNumber")
        opt$featuredata <- left_join(opt$featuredata, CDS2EC)
        opt$featuredata$ECNumber[is.na(opt$featuredata$ECNumber)] <- "none"
    } else {
        opt$featuredata$ECNumber <- "none"
    }

    #Add taxonomy information
    tmpcontigsdata <- opt$contigsdata[, c("Contig", "LKT")]
    opt$featuredata <- left_join(opt$featuredata, tmpcontigsdata)

    #Clean up
    file.remove(file.path(annotationfolder, c("feature2contig.map", "feature2product.map", "CDS2EC.map")))

    #Bank annotation files to project directory
    if (opt$workdir != opt$sampledir){
        flog.info("Banking annotation to sample directory.")
        copyargs <- c("-R", file.path(opt$workdir, annotationfolder), opt$sampledir)
        system(paste0("cp", copyargs, collapse = " "))
    }

    return(opt)
}
