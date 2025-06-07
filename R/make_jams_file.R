#' make_jams_file
#'
#' JAMSalpha function
#' @export

make_jams_file <- function(opt = opt, objectstoinclude = c("contigsdata", "featuredata", "abundances", "taxonomic_Quality_split_list", opt$blastanalyses, "bam_contig_depths", "TNF_contigs", "readstats", "fastqstats", "projinfo"), outdir = opt$outdir){

    setwd(opt$sampledir)

    flog.info("Making .jams file")
    JAMSobj <- objectstoinclude[objectstoinclude %in% names(opt)]

    filesuffix <- "rds"

    JAMSobjfiles <- paste(paste(opt$prefix, JAMSobj, sep = "_"), filesuffix, sep = ".")
    jamsfile <- file.path(outdir, paste(opt$prefix, "jams", sep = "."))
    jamsargs <- c("-zcvf", jamsfile, JAMSobjfiles)

    system2('tar', args = jamsargs, stdout = TRUE, stderr = TRUE)

}
