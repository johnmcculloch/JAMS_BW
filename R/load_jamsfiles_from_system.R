#' load_jamsfiles_from_system(path = ".", recursive = TRUE, onlysamples = NULL, loadfromscratch = TRUE, list.data = NULL, threads = 8))
#' JAMSbeta function
#'
#' Loads all JAMS files from system
#' @export

load_jamsfiles_from_system <- function(path = ".", recursive = TRUE, onlysamples = NULL, threads = 8, multithread_decomp = TRUE){
    fljams <- list.files(path = path, pattern = "*.jams$", full.names = TRUE, recursive = recursive)

    if (length(fljams) < 1){
        stop("There are no .jams files in the specified path.")
    }

    currpath <- getwd()
    #Create temp files directory if needed.
    jamstempfilespath <- file.path(currpath, "jamstempfiles")
    if (!(dir.exists(jamstempfilespath))){
        flog.info("Creating directory to hold raw data.")
        dir.create(jamstempfilespath, showWarnings = FALSE, recursive = FALSE)
    }

    #Make a dataframe with jamsfiles whereabouts
    extract_prefix <- function(fullpath){
        prefix <- tail(unlist(strsplit(fullpath, split = "/")), n = 1)
        prefix <- gsub(".jams", "", prefix)

        return(prefix)
    }
    jamsfilesdf <- data.frame(Prefix = sapply(1:length(fljams), function(x) { extract_prefix(fljams[x]) }), FullPath = fljams, stringsAsFactors = FALSE)

    #Restrict loading to only certain samples if required.
    if (!is.null(onlysamples)){
        flog.info(paste("There are", length(onlysamples), "samples to load."))
        flog.info(paste("Will only load .jams files for samples", paste0(onlysamples, collapse = ", ")))
        jamsfilesdfwant <- subset(jamsfilesdf, Prefix %in% onlysamples)
        samplesIhavejams <- jamsfilesdfwant$Prefix
        if (length(onlysamples %in% samplesIhavejams) < length (onlysamples)){
            missingjamsfiles <- paste0(onlysamples[which(!(onlysamples %in% samplesIhavejams))], collapse = ", ")
            stop(paste("Could not find .jams files for samples", missingjamsfiles))
        }
    }

    fpfljams <- jamsfilesdfwant$FullPath

    #If there are .jams files, expand them.
    if (length(fpfljams) > 0){
        flog.info("Please be patient. Depending on how much data you have this might take a while.")
        flog.info(paste("There are", length(fpfljams), "objects to expand."))

        decompress_jamsfile <- function (x){
            untar(jamsfilesdfwant$FullPath[x], list = FALSE, exdir = jamstempfilespath, compressed = TRUE, verbose = TRUE)
        }
        if (multithread_decomp == TRUE){
            appropriatenumcores <-  max((threads - 2), 1)
            explist <- mclapply(1:nrow(jamsfilesdfwant), function (x) { decompress_jamsfile(x) }, mc.cores = appropriatenumcores)
        } else {
            for (f in 1:nrow(jamsfilesdfwant)){
                flog.info(paste("Expanding ", jamsfilesdfwant$Prefix[f], ", file ", f, "/", nrow(jamsfilesdfwant), "...", sep = ""))
                untar(jamsfilesdfwant$FullPath[f], list = FALSE, exdir = jamstempfilespath, compressed = TRUE, verbose = TRUE)
            }
        }
    }

    fl <- list.files(path = jamstempfilespath, pattern = "*.tsv$")
    on <- gsub("*.tsv$", "", fl)
    flwp <- file.path(jamstempfilespath, fl)

    flog.info("Please be patient. Depending on how much data you have this might take a while.")

    lastpos <- length(names(list.data))

    flog.info(paste("There are", length(fl), "objects to load."))

    loadtsv <- function(tsvfn = NULL){
        jamsdf <- fread(data.table = FALSE, file = tsvfn, sep = "\t", header = TRUE, quote = "", fill = TRUE, integer64 = "numeric", logical01 = FALSE, stringsAsFactors = FALSE, nThread = 1)

        return(jamsdf)
    }

    appropriatenumcores <- max(1 , (min((opt$threads - 2), length(flwp))))
    flog.info(paste("Using", appropriatenumcores, "CPUs to load objects."))
    list.data <- mclapply(1:length(flwp), function (x) { loadtsv(tsvfn = flwp[x]) }, mc.cores = appropriatenumcores)

    names(list.data) <- on
    unlink(jamstempfilespath, recursive = TRUE)

    flog.info("Finished loading all objects.")

    return(list.data)
}
