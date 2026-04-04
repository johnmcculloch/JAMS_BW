#' load_jamsfiles_from_system(path = ".", recursive = TRUE, onlysamples = NULL, onlyobjects = NULL, threads = 8, multithread_decomp = TRUE, multithread_load = TRUE, use_exactly_n_threads = NULL)
#' JAMSbeta function
#'
#' Loads all JAMS files from system
#' @export

load_jamsfiles_from_system <- function(path = ".", recursive = TRUE, onlysamples = NULL, onlyobjects = NULL, threads = 8, multithread_decomp = TRUE, multithread_load = TRUE, use_exactly_n_threads = NULL, tempdir = NULL){

    require(parallel)
    flog.info("Searching for jams files.")
    fljams <- list.files(path = path, pattern = "*.jams$", full.names = TRUE, recursive = recursive, include.dirs = TRUE)
    if (length(fljams) < 1){
        stop("There are no .jams files in the specified path.")
    }

    fljamsnames <- sapply(1:length(fljams), function(x) { tail(unlist(strsplit(fljams[x], split = "/")), n = 1) })

    #Make a dataframe with all the info
    jamsfilesdf <- data.frame(Filename = fljamsnames, FullPath = fljams, stringsAsFactors = FALSE)
    jamsfilesdf$FullPath <- sapply(1:length(jamsfilesdf$FullPath), function(x) { fixrelpath(jamsfilesdf$FullPath[x]) })
    jamsfilesdf$Prefix <- gsub(".jams$", "", jamsfilesdf$Filename)

    if (is.null(tempdir)){
        currpath <- getwd()
    } else {
        currpath <- tempdir
    }

    #Create temp files directory if needed.
    jamstempfilespath <- file.path(currpath, "jamstempfiles")
    #if (!(dir.exists(jamstempfilespath))){
        flog.info(paste("Creating directory to hold raw data at", jamstempfilespath))
        dir.create(jamstempfilespath, showWarnings = FALSE, recursive = FALSE)
    #}

    #Restrict loading to only certain samples if required.
    if (!is.null(onlysamples)){
        flog.info(paste("There are", length(onlysamples), "samples to load."))
        #Avoid being overly verbose, only report sample names if < 50 samples.
        if (length(onlysamples) <= 50){
            flog.info(paste("Will only load .jams files for samples", paste0(onlysamples, collapse = ", ")))
        }
        jamsfilesdfwant <- subset(jamsfilesdf, Prefix %in% onlysamples)
        samplesIhavejams <- jamsfilesdfwant$Prefix
        if (sum(onlysamples %in% samplesIhavejams) < length (onlysamples)){
            missingjamsfiles <- paste0(onlysamples[which(!(onlysamples %in% samplesIhavejams))], collapse = ", ")
            flog.warn(paste("Could not find .jams files for sample(s):", missingjamsfiles))
        }
    } else {
        jamsfilesdfwant <- jamsfilesdf
    }

    fpfljams <- jamsfilesdfwant$FullPath

    #If there are .jams files, expand them.
    if (length(fpfljams) > 0){
        #flog.info("Please be patient. Depending on how much data you have this might take a while.")
        flog.info(paste("There are", length(fpfljams), "objects to expand."))

        decompress_jamsfile <- function (x){
            untar(jamsfilesdfwant$FullPath[x], list = FALSE, exdir = jamstempfilespath, verbose = FALSE)
        }

        if (multithread_decomp == TRUE){
            if (is.null(use_exactly_n_threads)){
                appropriatenumcores <- max((threads - 2), 1)
            } else {
                appropriatenumcores <- as.numeric(use_exactly_n_threads)
            }
            explist <- mclapply(1:nrow(jamsfilesdfwant), function (x) { decompress_jamsfile(x) }, mc.cores = appropriatenumcores)
        } else {
            for (f in 1:nrow(jamsfilesdfwant)){
                flog.info(paste("Expanding ", jamsfilesdfwant$Prefix[f], ", file ", f, "/", nrow(jamsfilesdfwant), "...", sep = ""))
                untar(jamsfilesdfwant$FullPath[f], list = FALSE, exdir = jamstempfilespath, compressed = TRUE, verbose = TRUE)
            }
        }
    }

    fl <- list.files(path = jamstempfilespath, pattern = "\\.rds$")
    ## No longer load ucobias, TNF_contigs, TNF_features, taxa_16S_cons as these are not used.
    ##
    fl <- fl[grep("_ucobias\\.", fl, invert = TRUE)]
    #fl <- fl[grep("_TNF_contigs\\.", fl, invert = TRUE)]
    fl <- fl[grep("_TNF_features\\.", fl, invert = TRUE)]
    fl <- fl[grep("_taxa_16S_cons\\.", fl, invert = TRUE)]

    if (!is.null(onlyobjects)){
        flog.warn(paste("Only loading JAMS objects matching", paste0(onlyobjects, collapse = ", ")))
        fl <- unlist(lapply(paste0("_", onlyobjects, "\\.rds"), function (x) { fl[grep(x, fl)] }))
    }

    on <- gsub("*.rds$", "", fl)
    if (length(which(duplicated(on) == TRUE)) > 0){
        stop("There are duplicated object names. Were different versions of the same jams file loaded? Check files and try again. Aborting now.")
        q()
    }
    flwp <- file.path(jamstempfilespath, fl)
    names(flwp) <- on

    flog.info("Please be patient. Depending on how much data you have this might take a while.")

    flog.info(paste("There are", length(fl), "objects to load."))

    if (threads < 8){
        flog.warn("It is ill-advised to load JAMS objects with multiple threads with less than 8 CPUs. Defaulting to single threaded loading. Please be patient while all objects are loaded into memory.")
        multithread_load <- FALSE
    }

    if (multithread_load) {

        if (is.null(use_exactly_n_threads)){
            appropriatenumcores <- max(1 , (min((threads - 2), length(flwp))))
        } else {
            appropriatenumcores <- as.numeric(use_exactly_n_threads)
        }

        flog.info(paste("Using", appropriatenumcores, "CPUs to load objects."))
        #Revert to this commented line in case the makeCluster approach below doesn't work out.
        #list.data <- mclapply(flwp, function (x) { loadobj(objfn = x, never_used_columns = c("AntiFam", "CDD", "Coils", "FunFam", "Gene3D", "PANTHER", "Phobius", "ProSitePatterns", "SFLD", "SignalP_EUK", "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE", "SMART", "TIGRFAM", "MetaCyc")) }, mc.cores = appropriatenumcores, mc.preschedule = FALSE, mc.cleanup = TRUE)
 
        # Make cluster
        cl <- parallel::makeCluster(appropriatenumcores)
        # Explicitly export the necessary function and objects to each worker.
        parallel::clusterExport(cl, varlist = c("loadobj", "fread", "readRDS", "flog.warn", "setDTthreads"))
        list.data <- parallel::parLapply(cl, flwp, function(x) {
            data.table::setDTthreads(1)
            loadobj(objfn = x)
        })
        # Stop the cluster
        suppressWarnings(parallel::stopCluster(cl))

    } else {
        list.data <- lapply(flwp, function (x) { loadobj(objfn = x) })
    }

    #names(list.data) <- on #no longer needed because on was attributed as names of vector flwp.

    #Flush out any object that was not loaded properly
    empties <- sapply(1:length(on), function(x){ is.null(list.data[[on[x]]])})
    if (sum(empties) > 0){
        flog.warn(paste("Unable to load the following objects:", paste0(on[empties], collapse = ", ")))
    }

    validon <- on[!empties]
    list.data <- list.data[validon]

    #Check if objects were loaded into list correctly. Not necessary, but one is better safe than sorry.
    prefixespresent <- sapply(validon, function (x) { tail(rev(unlist(strsplit(x, split = "_"))), n = 1)} )
    prefixespresent <- unique(unname(prefixespresent))
    #minobjlist <- as.vector(sapply(c("projinfo", "contigsdata", "featuredata", "LKTdose", "featuredose"), function(x) { paste(prefixespresent, x, sep = "_") }))
    #Because of MetaPhlAnn
    minobjlist <- as.vector(sapply(c("projinfo", "abundances"), function(x) { paste(prefixespresent, x, sep = "_") }))

    if (!all(minobjlist %in% names(list.data))){
        missingobjects <- minobjlist[!minobjlist %in% names(list.data)]
        incompleteprefixes <- unique(unname(sapply(missingobjects, function (x) { tail(rev(unlist(strsplit(x, split = "_"))), n = 1)} )))
        flog.warn(paste("Sample(s)", paste0(incompleteprefixes, collapse = ", "), "will be omitted because the following necessary objects are missing:", paste0(missingobjects, collapse = ", "), ". Please check the integrity of your jams files and try again."))
        objectstoremove <- as.vector(sapply(c("projinfo", "contigsdata", "featuredata", "LKTdose", "featuredose"), function(x) { paste(incompleteprefixes, x, sep = "_") }))
        validon <- validon[!(validon %in% objectstoremove)]
        minobjlist <- minobjlist[minobjlist %in% validon]
        list.data <- list.data[validon]
    }

    file.remove(flwp)
    unlink(jamstempfilespath, recursive = TRUE)
    gc()

    flog.info("Finished loading all objects.")

    return(list.data)
}
