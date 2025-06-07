#' launch_interpro
#'
#' JAMSalpha function
#' @export

launch_interpro <- function(opt = NULL){

    #Set working directory to sampledir
    setwd(opt$sampledir)

    #Discard proteins with stop codons in them, as this blocks interpro
    CDSswithstop <- names(which(sapply(names(opt$proteome), function (x) { length(grep("\\*", opt$proteome[[x]])) } ) > 0))

    proteomeforipro <- filter_sequence_by_name(input_sequences = opt$proteome, sequencenames = CDSswithstop, keep = FALSE)

    #Write a copy of contigs to be annotated to the system
    flog.info("Writing proteome for Interproscan annotation.")
    write.fasta(sequences = proteomeforipro, names = names(proteomeforipro), nbchar = 80, file.out = "proteome.faa")
    if(length(opt$proteome) > 10000){
        chksize = 1000
    } else {
        chksize = 500
    }

    #annotate with prokka
    flog.info("Launching Interproscan analysis of proteome.")
    interprocmd <- file.path(opt$bindir, "dointerproBW.sh")
    interproargs <- c("-i", "proteome.faa", "-p", opt$prefix, "-c", chksize)
    system2(interprocmd, args = interproargs, stdout = TRUE, stderr = TRUE)

    opt$iprodir <- file.path(opt$sampledir, "interproscan")

    file.remove("proteome.faa")

    return(opt)
}
