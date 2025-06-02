#' bank_reads
#'
#' JAMSalpha function
#' @export

bank_reads <- function(opt = opt, maketarball = FALSE){

   #Check if reads are in tmpdir
    if(opt$workdir != opt$sampledir){
        readsworkdir <- file.path(opt$workdir, "reads")
    } else {
        readsworkdir <- opt$readsdir
    }
    setwd(readsworkdir)

    if(opt$host != "none"){
        flog.info("Banking Non-Aligned to Host Species (microbiota) reads to outdir.")
        fastqstobank <- sort(opt$nahsreads)
        if (maketarball){
            outfilename <- paste(paste(opt$prefix, "Microbiota", "reads", sep="_"), "tar", "gz", sep=".")
            outfilepath <- file.path(opt$outdir, outfilename)
            readbankcmds <- paste("tar", "cf", "-", paste0(fastqstobank, collapse=" "), "|", "pigz", "-9", "-p", opt$threads, ">", outfilepath, collapse=" ")
        } else {
            fnsuffixes <- paste((paste0("R", 1:length(fastqstobank))), "fastq", "gz", sep = ".")
            outfilenames <- paste(opt$prefix, "Microbiota", "reads", fnsuffixes, sep="_")
            readbankcmds <- sapply(1:length(fastqstobank), function (x) { paste("pigz --best --stdout --processes", opt$threads, fastqstobank[x], ">", file.path(opt$outdir, outfilenames[x]), collapse = " ") })
        }

        sapply(readbankcmds, function (x) { system(x) } )


    } else {
        flog.info("Banking QC trimmed reads to outdir.")
        if (maketarball){
            fastqstobank <- sort(opt$trimreads)
            outfilename <- paste(paste(opt$prefix, "QCtrimmed", "reads", sep="_"), "tar", "gz", sep = ".")
            readbankcmds <- paste("tar", "cf", "-", paste0(fastqstobank, collapse=" "), "|", "pigz", "-9", "-p", opt$threads, ">", outfilepath, collapse=" ")
        } else {
            fnsuffixes <- paste((paste0("R", 1:length(fastqstobank))), "fastq", "gz", sep = ".")
            outfilenames <- paste(opt$prefix, "QCtrimmed", "reads", fnsuffixes, sep="_")
            readbankcmds <- sapply(1:length(fastqstobank), function (x) { paste("pigz --best --stdout --processes", opt$threads, fastqstobank[x], ">", file.path(opt$outdir, outfilenames[x]), collapse = " ") })
        }
    }

    system(tarcmd)

}
