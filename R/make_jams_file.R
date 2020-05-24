#' make_jams_file
#'
#' JAMSalpha function
#' @export

make_jams_file <- function(opt = opt, outdir = opt$outdir){
    setwd(opt$sampledir)
    file.info("Making .jams file")
    JAMSobj <- c("projinfo", "contigsdata", "featuredata", "assemblystats", "LKTdose", "featuredose", "readstats", "ucobias", "taxa_16S_cons", "resfinder", "plasmidfinder", "abricate", "vfdb")[c("projinfo","contigsdata", "featuredata", "assemblystats", "LKTdose", "featuredose", "readstats", "ucobias","taxa_16S_cons", "resfinder", "plasmidfinder", "abricate", "vfdb") %in% names(opt)]
    JAMStsvs <- paste(paste(opt$prefix, JAMSobj, sep = "_"), "tsv", sep = ".")
    jamsfile <- file.path(outdir, paste(opt$prefix, "jams", sep = "."))
    jamsargs <- c("-zcvf", jamsfile, JAMStsvs)

    system2('tar', args = jamsargs, stdout = TRUE, stderr = TRUE)
}
