#' load_jamsfiles_from_system(path = ".", recursive = TRUE, onlysamples = NULL, threads = 8, multithread_decomp = TRUE, multithread_load = TRUE, use_exactly_n_threads = NULL)
#' JAMSbeta function
#'
#' Loads all JAMS files from system
#' @export

load_jamsfiles_from_system <- function(path = ".", recursive = TRUE, onlysamples = NULL, threads = 8, multithread_decomp = TRUE, multithread_load = TRUE, use_exactly_n_threads = NULL){

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

    currpath <- getwd()
    #Create temp files directory if needed.
    jamstempfilespath <- file.path(currpath, "jamstempfiles")
    if (!(dir.exists(jamstempfilespath))){
        flog.info("Creating directory to hold raw data.")
        dir.create(jamstempfilespath, showWarnings = FALSE, recursive = FALSE)
    }

    #Restrict loading to only certain samples if required.
    if (!is.null(onlysamples)){
        flog.info(paste("There are", length(onlysamples), "samples to load."))
        flog.info(paste("Will only load .jams files for samples", paste0(onlysamples, collapse = ", ")))
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
        flog.info("Please be patient. Depending on how much data you have this might take a while.")
        flog.info(paste("There are", length(fpfljams), "objects to expand."))

        decompress_jamsfile <- function (x){
            untar(jamsfilesdfwant$FullPath[x], list = FALSE, exdir = jamstempfilespath, compressed = TRUE, verbose = TRUE)
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

    fl <- list.files(path = jamstempfilespath, pattern = "\\.tsv$|\\.rds$")
    ## No longer load ucobias, TNF_contigs, TNF_features, taxa_16S_cons as these are not used.
    ##
    fl <- fl[grep("_ucobias\\.", fl, invert = TRUE)]
    fl <- fl[grep("_TNF_contigs\\.", fl, invert = TRUE)]
    fl <- fl[grep("_TNF_features\\.", fl, invert = TRUE)]
    fl <- fl[grep("_taxa_16S_cons\\.", fl, invert = TRUE)]

    on <- gsub("*.tsv$", "", fl)
    on <- gsub("*.rds$", "", on)
    if (length(which(duplicated(on) == TRUE)) > 0){
        stop("There are duplicated object names. Were different versions of the same jams file loaded? Check files and try again. Aborting now.")
        q()
    }
    flwp <- file.path(jamstempfilespath, fl)

    flog.info("Please be patient. Depending on how much data you have this might take a while.")

    flog.info(paste("There are", length(fl), "objects to load."))

    loadobj <- function(objfn = NULL, never_used_columns = c("AntiFam", "CDD", "Coils", "FunFam", "Gene3D", "PANTHER", "Phobius", "ProSitePatterns", "SFLD", "SignalP_EUK", "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE", "SMART", "SUPERFAMILY", "TIGRFAM", "MetaCyc")){
        #decide whether it is a TSV or RDS file and load appropriately.
        if ((length(grep("\\.rds$", objfn)) > 0)){
            tryCatch((jamsdf <- readRDS(objfn)), error = function() { flog.warn(paste("Unable to read file", objfn)) } )
            #file is rds
        } else {
            #assume tsv
            jamsdf <- fread(data.table = FALSE, file = objfn, sep = "\t", header = TRUE, quote = "", fill = FALSE, integer64 = "numeric", logical01 = FALSE, stringsAsFactors = FALSE, nThread = 1)
        }

        #Prune unwanted columns
        if (any(never_used_columns %in% colnames(jamsdf))){
            wantedcols <- colnames(jamsdf)[!(colnames(jamsdf) %in% never_used_columns)]
            jamsdf <- jamsdf[ , wantedcols]
        }

        if ("Analysis" %in% colnames(jamsdf)){
            jamsdf <- jamsdf[!(jamsdf$Analysis %in% never_used_columns), ]
        }

        return(jamsdf)
    }

    if (threads < 16){
        flog.warn("It is ill-advised to load JAMS objects with multiple threads with less than 16 CPUs. Defaulting to single threaded loading. Please be patient while all objects are loaded into memory.")
        multithread_load <- FALSE
    }

    if (multithread_load) {

        if (is.null(use_exactly_n_threads)){
            appropriatenumcores <- max(1 , (min((threads - 2), length(flwp))))
        } else {
            appropriatenumcores <- as.numeric(use_exactly_n_threads)
        }

        flog.info(paste("Using", appropriatenumcores, "CPUs to load objects."))
        list.data <- mclapply(flwp, function (x) { loadobj(objfn = x) }, mc.cores = appropriatenumcores, mc.preschedule = FALSE, mc.cleanup = TRUE)

    } else {
        list.data <- lapply(flwp, function (x) { loadobj(objfn = x) })
    }

    names(list.data) <- on
    #Flush out any object that was not loaded properly
    empties <- sapply(1:length(on), function(x){ is.null(list.data[[on[x]]])})
    if(sum(empties) > 0){
        flog.warn(paste("Unable to load the following objects:", paste0(on[empties], collapse = ", ")))
    }

    validon <- on[!empties]
    list.data <- list.data[validon]

    #Check if objects were loaded into list correctly. Not necessary, but one is better safe than sorry.
    prefixespresent <- sapply(validon, function (x) { tail(rev(unlist(strsplit(x, split = "_"))), n = 1)} )
    prefixespresent <- unique(unname(prefixespresent))
    #minobjlist <- as.vector(sapply(c("projinfo", "contigsdata", "featuredata", "LKTdose", "featuredose"), function(x) { paste(prefixespresent, x, sep = "_") }))
    #Because of MetaPhlAnn
    minobjlist <- as.vector(sapply(c("projinfo", "LKTdose"), function(x) { paste(prefixespresent, x, sep = "_") }))

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
