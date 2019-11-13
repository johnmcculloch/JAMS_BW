#' load_jamsfiles_from_system(path = ".", recursive = TRUE, onlysamples = NULL, loadfromscratch = TRUE, list.data = NULL, threads = 8))
#' JAMSbeta function
#'
#' Loads all JAMS files from system
#' @export

load_jamsfiles_from_system <- function(path = ".", recursive = TRUE, onlysamples = NULL, loadfromscratch = TRUE, list.data = NULL, threads = 8){
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

    #If not loading from scratch, check for existing data to append to.
    if (loadfromscratch != TRUE){
        if (is.null(list.data)){
            stop("If you want to append data to an existing project, you must supply a list.data object to append to.  If you want to load everything from scratch, then set loadfromscratch = TRUE.")
        }
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
        for (f in 1:nrow(jamsfilesdfwant)){
            flog.info(paste("Expanding ", jamsfilesdfwant$Prefix[f], ", file ", f, "/", nrow(jamsfilesdfwant), "...", sep = ""))
            untar(jamsfilesdfwant$FullPath[f], list = FALSE, exdir = jamstempfilespath, compressed = TRUE, verbose = TRUE)
        }
    }

    fl <- list.files(path = jamstempfilespath, pattern = "*.tsv$")
    on <- gsub("*.tsv$", "", fl)
    flwp <- file.path(jamstempfilespath, fl)

    if (is.null(list.data)){
        list.data <- list()
    }

    flog.info("Please be patient. Depending on how much data you have this might take a while.")

    lastpos <- length(names(list.data))

    flog.info(paste("There are", length(fl), "objects to load."))
    for (i in 1:length(fl)){
        flog.info(paste("Loading ", fl[i], ", file ", i, "/", length(fl), "...", sep = ""))
        list.data[[i+lastpos]] <- read.table(file = flwp[i], sep = "\t", header = TRUE, quote = "", skipNul = FALSE, fill = TRUE, stringsAsFactors = FALSE)
        #list.data[[i+lastpos]] <- fread(file = flwp[i], sep = "\t", header = TRUE, quote = "", fill = TRUE, stringsAsFactors = FALSE, nThread = threads)
        names(list.data)[i + lastpos] <- on[i]
    }

    unlink(jamstempfilespath, recursive = TRUE)

    return(list.data)
}
