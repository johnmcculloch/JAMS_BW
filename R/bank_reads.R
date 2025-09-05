#' bank_reads
#'
#' JAMSalpha function
#' @export

bank_reads <- function(opt = opt){

   #Check if reads are in tmpdir
    if (opt$workdir != opt$sampledir){
        readsworkdir <- file.path(opt$workdir, "reads")
    } else {
        readsworkdir <- opt$readsdir
    }
    setwd(readsworkdir)

    if (opt$host != "none"){

        flog.info("Banking Non-Aligned to Host Species (microbiota) reads to outdir.")
        fastqstobank <- sort(opt$nahsreads)
        fnsuffixes <- paste((paste0("R", 1:length(fastqstobank))), "fastq", "gz", sep = ".")
        outfilenames <- paste(opt$prefix, "Microbiota", "reads", fnsuffixes, sep = "_")

    } else {

        flog.info("Banking QC trimmed reads to outdir.")
        fastqstobank <- sort(opt$trimreads)
        fnsuffixes <- paste((paste0("R", 1:length(fastqstobank))), "fastq", "gz", sep = ".")
        outfilenames <- paste(opt$prefix, "QCtrimmed", "reads", fnsuffixes, sep="_")

    }

    #Execute banking commands. Files are gz compressed, so just copy
    for (fln in 1:length(fastqstobank)){
        file.copy(from = file.path(opt$readsdir, fastqstobank[fln]), to = file.path(opt$outdir, outfilenames[fln]))
    }

    #Back to sample folder
    setwd(opt$sampledir)
}
