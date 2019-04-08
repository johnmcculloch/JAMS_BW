#' bank_reads
#'
#' JAMSalpha function
#' @export

bank_reads<-function(opt=opt){

   #Check if reads are in tmpdir
    if(opt$workdir != opt$sampledir){
        readsworkdir<-file.path(opt$workdir, "reads")
    } else {
        readsworkdir<-opt$readsdir
    }
    setwd(readsworkdir)

    if(opt$host != "none"){
        flog.info("Banking Non-Aligned to Host Species (microbiota) reads to outdir.")
        fastqstobank<-opt$nahsreads
        outfilename<-paste(paste(opt$prefix, "Microbiota", "reads", sep="_"), "tar", "gz", sep=".")
    } else {
        flog.info("Banking QC trimmed reads to outdir.")
        fastqstobank<-opt$trimreads
        outfilename<-paste(paste(opt$prefix, "QCtrimmed", "reads", sep="_"), "tar", "gz", sep=".")
    }

    outfilepath <- file.path(opt$outdir, outfilename)
    tarcmd<-paste("tar", "cf", "-", paste0(fastqstobank, collapse=" "), "|", "pigz", "-9", "-p", opt$threads, ">", outfilepath, collapse=" ")
    system(tarcmd)

}
