#' ExpObjVetting(ExpObj = NULL)
#'
#' Performs vetting of a SummarizedExperiment object for use in several functions
#' @export

ExpObjVetting <- function(ExpObj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, glomby = NULL, variables_to_fix = NULL, class_to_ignore = NULL){

        #Get appropriate object to work with
        if (as.character(class(ExpObj)) != "SummarizedExperiment"){
            stop("This function can only take a SummarizedExperiment object as input. For using a metagenomeSeq object (deprecated), please use plot_relabund_heatmap_mgseq.")
        }
        obj <- ExpObj

        if (!(is.null(glomby))){
            obj <- agglomerate_features(ExpObj = obj, glomby = glomby)
        }

        #Exclude samples and features if specified
        if (!(is.null(samplesToKeep))){
            samplesToKeep <- samplesToKeep[samplesToKeep %in% colnames(obj)]
            obj <- obj[, samplesToKeep]
        }

        if (!(is.null(featuresToKeep))){
            featuresToKeep <- featuresToKeep[featuresToKeep %in% rownames(obj)]
            obj <- obj[featuresToKeep, ]
        }

        obj <- suppressWarnings(filter_sample_by_class_to_ignore(SEobj = obj, variables = variables_to_fix, class_to_ignore = class_to_ignore))

    return(obj)
}

#' declare_filtering_presets(analysis = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, maxl2fc = NULL, minl2fc = NULL)
#'
#' Performs vetting of a SummarizedExperiment object for use in several functions
#' @export
declare_filtering_presets <- function(analysis = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, maxl2fc = NULL, minl2fc = NULL){

    if ((analysis != "LKT") && (!(is.null(GenomeCompletenessCutoff)))){
        warning("Genome completeness only makes sense for taxa. Please choose a taxonomic (non functional) analysis.")
        GenomeCompletenessCutoff <- NULL
    }

    presetlist <- list()

    if (!is.null(applyfilters)){
        if (applyfilters == "stringent"){
            if (analysis == "LKT"){
                presetlist$featcutoff <- c(2000, 15)
                presetlist$GenomeCompletenessCutoff <- c(30, 10)
                presetlist$PctFromCtgscutoff <- c(50, 10)
                presetlist$minl2fc <- 2
            } else {
                presetlist$featcutoff <- c(50, 15)
                presetlist$minl2fc <- 2.5
            }
        } else if (applyfilters == "moderate"){
            if (analysis == "LKT"){
                presetlist$featcutoff <- c(250, 15)
                presetlist$GenomeCompletenessCutoff <- c(10, 5)
                presetlist$PctFromCtgscutoff <- c(25, 10)
                presetlist$minl2fc <- 1
            } else {
                presetlist$featcutoff <- c(10, 5)
                presetlist$minl2fc <- 1
            }
        } else if (applyfilters == "light"){
            if (analysis == "LKT"){
                presetlist$featcutoff <- c(50, 5)
                presetlist$GenomeCompletenessCutoff <- c(5, 5)
                presetlist$PctFromCtgscutoff <- c(25, 10)
                presetlist$minl2fc <- 1
            } else {
                presetlist$featcutoff <- c(5, 5)
                presetlist$minl2fc <- 1
            }
        }
    }

    #Replace with any values explicitly set by the user
    argstoset <- c("featcutoff", "GenomeCompletenessCutoff", "PctFromCtgscutoff", "maxl2fc", "minl2fc",)[!unlist(lapply(list(featcutoff, GenomeCompletenessCutoff, PctFromCtgscutoff, maxl2fc, minl2fc), is.null))]

    if (length(argstoset) > 0){
        for (ats in argstoset){
            presetlist[[ats]] <- get(ats)
        }
    }

    #Generate a filtration message
    presetlist$filtermsg <- NULL
    #Discard features which do not match certain criteria
    if (!(is.null(presetlist$featcutoff))){
        presetlist$thresholdPPM <- presetlist$featcutoff[1]
        presetlist$sampcutoffpct <- min(presetlist$featcutoff[2], 100)
        presetlist$filtermsg <- paste("Feature must be >", presetlist$thresholdPPM, "PPM in at least ", presetlist$sampcutoffpct, "% of samples", sep = "")
    } else {
        presetlist$filtermsg <- NULL
        presetlist$featcutoff <- c(0, 0)
    }

    if (!(is.null(presetlist$PctFromCtgscutoff))){
        presetlist$thresholdPctFromCtgs <- presetlist$PctFromCtgscutoff[1]
        presetlist$sampcutoffpctPctFromCtgs <- min(presetlist$PctFromCtgscutoff[2], 100)
        presetlist$filtermsg <- paste(presetlist$filtermsg, (paste("Taxonomy information must come from >", presetlist$thresholdPctFromCtgs, "% contigs in at least ", presetlist$sampcutoffpctPctFromCtgs, "% of samples", sep = "")), sep = "\n")
    }

    if (!(is.null(presetlist$GenomeCompletenessCutoff))){
        presetlist$thresholdGenomeCompleteness <- presetlist$GenomeCompletenessCutoff[1]
        presetlist$sampcutoffpctGenomeCompleteness <- min(presetlist$GenomeCompletenessCutoff[2], 100)
        presetlist$filtermsg <- paste(presetlist$filtermsg, (paste("Taxon genome completeness must be >", presetlist$thresholdGenomeCompleteness, "% in at least ", presetlist$sampcutoffpctGenomeCompleteness, "% of samples", sep = "")), sep = "\n")
    }

    return(presetlist)
}


#' log2PPMtoPct(log2PPM)
#'
#' Returns a rounded percentage given a log2 transformed PPM value
#' @export

log2PPMtoPct <- function(log2PPM= NULL, signifdigits = 1){
    PPM <- ((2 ^ log2PPM) - 1)
    Pct <- round((PPM / 10000), signifdigits)

    return(Pct)
}

#' Pct2log2PPM(Pct)
#'
#' Returns a log2 transformed PPM value given a percentage
#' @export

Pct2log2PPM <- function(Pct){
    PPM <- Pct * 10000
    log2PPM <- log2(PPM + 1)

    return(log2PPM)
}


#' convert_matrix_log2(mat = NULL, transformation = NULL)
#'
#' Wrapper for transforming a matrix to and from log2 formats
#' @export

convert_matrix_log2 <- function(mat = NULL, transformation = NULL){

    if (!(transformation %in% c("to_log2", "from_log2"))){
        stop("Please choose either \"to_log2\", or \"from_log2\" as transformation to apply to matrix.")
    }

    if (class(mat) != "matrix"){
        stop("Object to transform must be a matrix.")
    }

    if (transformation == "to_log2"){
        tfun <- function (x) {
            fx <- log2(x + 1)
            return(fx)
        }
    } else if (transformation == "from_log2"){
        tfun <- function (x) {
            fx <- ((2 ^ x) - 1)
            return(fx)
        }
    }

    transmat <- base::apply(mat, MARGIN = 2, FUN = tfun)

    return(transmat)

}


#' detectHardwareResources()
#' Generic function for getting number of CPUs and total memory on Biowulf or otherwise
#' Returns a named vector with available CPUs and Memory in bytes
#' @export

detectHardwareResources <- function(){
    #First off, detect if on Slurm type cluster
    #Get slurm job ID
    currslurmjobid <- as.character(Sys.getenv("SLURM_JOB_ID"))

    if(nchar(currslurmjobid) < 3){
       #Define appropriate functions for non-slurm system
       detectBatchCPUs <- function() {
            ncores <- detectCores()
            if (is.na(ncores)) {
                stop("Could not determine how many CPUs you have. Aborting.")
            }
            return(ncores)
        }

        detectAvailRAM <- function(){
            totmembytes<-as.numeric(get_ram())

            return(totmembytes)
        }

    } else {
        #Define appropriate functions for slurm system
        detectBatchCPUs <- function() {
            slurmjobid <- as.character(Sys.getenv("SLURM_JOB_ID"))
            ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

            if (is.na(ncores)) {
                #Try plan B
                sacctraw <- system2("sacct", args = c("-j", slurmjobid, "-X", "-P"), stdout = TRUE)
                jobinforaw <- sacctraw[2]
                jobinfoheaders <- sacctraw[1]
                jobinfo <- unlist(strsplit(jobinforaw, split="\\|"))
                names(jobinfo) <- unlist(strsplit(jobinfoheaders, split="\\|"))
                ncores <- as.integer(jobinfo["AllocCPUS"])
                print(jobinfo)
                if (is.na(ncores)) {
                    stop("Could not determine how many CPUs you have. Aborting.")
                }
            }

            return(ncores)
        }


        detectAvailRAM <- function(){
            mempercpu <- as.integer(Sys.getenv("SLURM_MEM_PER_CPU"))
            mempernode <- as.integer(Sys.getenv("SLURM_MEM_PER_NODE"))
            cpuspertask <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

            if(!(is.na(mempernode))){
                totmem <- mempernode
            } else {
                totmem <- mempercpu * cpuspertask
            }

            totmembytes <- totmem * 1000000

            return(totmembytes)
        }
    }
    hardwareRes <- NULL
    hardwareRes[1] <- detectBatchCPUs()
    hardwareRes[2] <- detectAvailRAM()
    names(hardwareRes) <- c("threads", "memory")

    return(hardwareRes)
}


#' whoopsieplot(msg = NULL)
#' Shuts down the device and gives a message on error, useful for when making reports.
#' @export
whoopsieplot <- function(msg = "trying to do this."){
    flog.info(paste("Whoops, something went wrong while", msg))
    dev.off()
}


#' filetype(path)
#' Wrapper for returning class of system file.
#' @export
filetype <- function(path){
    f <- file(path)
    ext <- summary(f)$class
    close.connection(f)
    ext
}

#' countfastq(fastqfile)
#' Wrapper for counting number of reads and number of bases of a fastq file in the system.
#' Returns a vector with read counts and base counts.
#' Caveat: Input fastq file must have exactly four lines per sequence. If your fastq files do not fit this criterion, you are totally bonkers.
#' @export

countfastq <- function(fastqfile){

    countargs <- c(fastqfile, "|", "paste", "-", "-", "-", "-", "|", "cut", "-f", "2", "|", "wc", "-lc")
    fastqstats <- system2('cat', args = countargs, stdout = TRUE, stderr = FALSE)
    fastqstats <- unlist(strsplit(fastqstats, split = " "))
    fastqstats <- fastqstats[which(fastqstats != "")]
    fastqstats <- as.numeric(fastqstats)

    return(fastqstats)
}


#' countfastq_files(fastqfiles = NULL, threads = NULL)
#' Wrapper for applying the countfastq() function to a vector of filenames using multiple threads.
#' Returns a dataframe with read counts and base counts for each input fastq file.
#' Caveat: Input fastq file must have exactly four lines per sequence. If your fastq files do not fit this criterion, you are totally bonkers.
#' @export

countfastq_files <- function(fastqfiles = NULL, threads = NULL){
    #fastqstatslist <- lapply(1:length(fastqfiles), function (x) { countfastq(fastqfiles[x]) })
    fastqstatslist <- mclapply(1:length(fastqfiles), function (x) { countfastq(fastqfiles[x]) }, mc.cores = threads)
    readcounts <- sapply(1:length(fastqstatslist), function (x) { fastqstatslist[[x]][1] })
    basecounts <- sapply(1:length(fastqstatslist), function (x) { fastqstatslist[[x]][2] })
    fastqstatsdf <- data.frame(Reads = fastqfiles, Count = readcounts, Bases = basecounts, stringsAsFactors = FALSE)
    fastqstatsdf$Readlength <- round((fastqstatsdf$Bases / fastqstatsdf$Count), 0)
    fastqstatsdf$Bases <- format(fastqstatsdf$Bases, scientific = FALSE)
    fastqstatsdf$Count <- format(fastqstatsdf$Count, scientific = FALSE)
    fastqstatsdf$Readlength <- format(fastqstatsdf$Readlength, scientific = FALSE)
    fastqstatsdf$Reads <- as.character(fastqstatsdf$Reads)
    fastqstatsdf$Count <- as.numeric(fastqstatsdf$Count)
    fastqstatsdf$Bases <- as.numeric(fastqstatsdf$Bases)
    fastqstatsdf$Readlength <- as.numeric(fastqstatsdf$Readlength)

    return(fastqstatsdf)
}


#' trim_whitespace_from_df(df)
#' Trims leading and trailing whitespace from a dataframe
#' @export

trim_whitespace_from_df <- function(df = NULL){
    #check if input is reasonable
    if (class(df) != "data.frame"){
        stop("Input must be a data.frame")
    }

    for (colm in 1:ncol(df)){
        df[ , colm] <- trimws(df[ , colm])
    }

    return(df)
}

#' exporttabletsv(dataobj = NULL, project = NULL, basename = NULL, row.names = TRUE, col.names = TRUE, path = NULL)
#' Generic function to export a dataframe to tsv format
#' @export

exporttabletsv <- function(dataobj = NULL, project = NULL, basename = NULL, row.names = TRUE, col.names = TRUE, path = NULL){
    flname <- paste(paste("JAMS", project, basename, sep = "_"), "tsv", sep = ".")
    if (!is.null(path)){
        flname <- paste(path, flname, sep = "/")
    }
    write.table(dataobj, file = flname, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = row.names, col.names = col.names)
}

#' quiet(function(x))
#' Suppresses output messages
#' By Hadley Wickham
#' http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
#'
#' @export

quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}


#' can_be_made_numeric(x, cats_to_ignore = NULL)
#' Returns logical for if a vector can be coerced into numeric.
#'
#' @export

can_be_made_numeric <- function(x, cats_to_ignore = NULL){

    if (!is.null(cats_to_ignore)){
        x <- x[which(!(x %in% cats_to_ignore))]
    }
    numtest <- length(which(is.na(suppressWarnings(as.numeric(x))) == TRUE))
    if (numtest == 0){
        testresult <- TRUE
    } else {
        testresult <- FALSE
    }

    return(testresult)
}


#' filter_correlations(corrmat = NULL, mincorrelcoeff = NULL)
#' Given a pairwise correlation matrix, eliminate features which do not present an absolute correlation coefficient smaller than mincorrelcoeff with all other features other than itself.
#'
#' @export
filter_correlations <- function(corrmat = NULL, mincorrelcoeff = NULL){

     if(nrow(corrmat) != ncol(corrmat)){
         stop("Correlation matrix must have equal numbers of rows and columns.")
     }

     featsIwant <- NULL

     for (rw in 1:nrow(corrmat)){
         featint <- rownames(corrmat)[rw]
         #print(paste("Checking:", featint))
         correlations <- corrmat[which(rownames(corrmat) != featint), featint]

         if(max(abs(correlations)) >= mincorrelcoeff){
             feat <- featint
         } else {
             feat <- NULL
         }

         featsIwant <- append(featsIwant, feat)

     }

     corrmat <- corrmat[featsIwant, featsIwant]

     return(corrmat)
 }

 #' IO_jams_workspace_image(opt = NULL, workspaceimage = NULL, threads = 8, operation = c("save", "load"))
 #' Safe way of either loading or saving an R workspace image. If argument workspaceimage is null, workspace image file will be searched for in opt (opt$projimage). If that is also NULL, saving or loading is aborted. If the fastSave package () is installed, multi-threaded loading or saving will be used. If opt is passed, number of CPUs will be set to opt$threads, trumping the threads argument.
 #'
 #' @export

IO_jams_workspace_image <- function(opt = NULL, workspaceimage = NULL, threads = 8, operation = NULL, verbose = FALSE){
    if (is.null(workspaceimage)){
        workspaceimage <- opt$projimage
    }

    if (is.null(workspaceimage)){
        stop("Workspace image not found or specified.")
    }

    if (!is.null(opt)){
        threads <- opt$threads
        RAMbytesavail <- opt$totmembytes
    } else {
        RAMbytesavail <- NULL
    }

    flog.info(paste("Workspace image is", workspaceimage))

    if (verbose){
        flog.info("Random Access Memory (RAM) status")
        print(RAMbytes_status(RAMbytesavail = RAMbytesavail))
    }

    if (operation == "save"){

        if ("fastSave" %in% rownames(installed.packages())){
            flog.info(paste("Saving project workspace image using fastSave package with", threads, "CPUs"))
            save.image.pigz(file = workspaceimage, n.cores = threads)
        } else {
            flog.info("Saving project workspace image. Please be patient...")
            save.image(workspaceimage)
        }

    } else if (operation == "load"){

        if ("fastSave" %in% rownames(installed.packages())){
            flog.info(paste("Loading project workspace image using fastSave package with", opt$threads, "CPUs"))
            load.pigz(file = workspaceimage, verbose = TRUE)
        } else {
            flog.info("Loading project workspace image. Please be patient...")
            save.image(workspaceimage)
        }

    } else {
        flog.info("You must choose between \"load\" or \"save\" as an operation.")
    }

    if (verbose){
        flog.info("Random Access Memory (RAM) status")
        print(RAMbytes_status(RAMbytesavail = RAMbytesavail))
    }
}

#' spew_heatmap_report(c("analysis", "report"))
#' Wrapper for launching a report for a SINGLE analysis
#'
#' @export

spew_heatmap_report <- function(hmcomb){
    analysis <- hmcomb[1]
    report <- hmcomb[2]
    setwd(file.path(opt$outdir, "Reports", analysis))

    if (report == "comparative"){

        tryCatch((make_heatmap_report(report = "comparative", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", showonlypbelow = 0.05, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating comparative heatmaps for", analysis, "analysis")))

    } else if (report == "exploratory") {

        tryCatch((make_heatmap_report(report = "exploratory", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 2, numthreads = 1, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating exploratory heatmaps for", analysis, "analysis")))

    } else if (report == "correlation") {

        tryCatch((make_heatmap_report(report = "correlation", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", minabscorrcoeff = 0.55, ntopvar = 250)), error = function(e) whoopsieplot(paste("generating correlation heatmaps for", analysis, "analysis")))

    } else if (report == "PA") {

        tryCatch((make_heatmap_report(report = "PA", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", showonlypbelow = 0.05, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating presence/absence heatmaps for", analysis, "analysis")))

    } else if ((report %in% c("tUMAP", "PCA", "tSNE"))){

        tryCatch((make_ordination_report(algorithm = report, project = project, expvec = expvec, usefulexp = analysis, applyfilters = "none", appendtofilename = paste(analysis, report, sep = "_"), variable_list = variable_list, doreads = opt$doreads, cdict = cdict, threads = 1, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating", report, "ordination plots for", analysis, "analysis")))

    }  else if (report == "alpha") {

        tryCatch(make_alpha_report(project = project, expvec = expvec, usefulexp = basicexp, variable_list = variable_list, cdict = cdict), error = function(e) whoopsieplot(paste("generating alpha diversity plots for", analysis, "analysis")))

    }

    setwd(opt$outdir)
}


#' RAMbytes_status(RAMbytesavail = NULL)
#' Reports how much RAM memory is being or has maximally been used.
#'
#' @export
RAMbytes_status <- function(RAMbytesavail = NULL){

    if (is.null(RAMbytesavail)){
        RAMbytesavail <- detectHardwareResources()["memory"]
    }

    memstats <- gc(full = TRUE)
    usedMbcolm <- (which(colnames(memstats) == "used") + 1)
    maxusedMbcolm <- (which(colnames(memstats) == "max used") + 1)
    usedRAMbytes <- sum(memstats[ , usedMbcolm]) * (1024 ^ 2)
    maxusedRAMbytes <- sum(memstats[ , maxusedMbcolm]) * (1024 ^ 2)
    RAMbytes <- unname(c(usedRAMbytes, maxusedRAMbytes, RAMbytesavail))
    RAMdf <- data.frame(Bytes = RAMbytes, Gbytes = round((RAMbytes / 1e9), 1), stringsAsFactors = FALSE)
    rownames(RAMdf) <- c("Used", "MaxUsed", "Available")
    RAMdf$ProportionAvail <- round((RAMdf$Bytes / RAMbytesavail), 2)

    return(RAMdf)
}


#' update_ExpObj_metadata(ExpObj = NULL, phenotable = NULL)
#'
#' Updates phenotable of a SummarizedExperiment object
#' @export

update_ExpObj_metadata <- function(ExpObj = NULL, phenotable = NULL){

        #Get appropriate object to work with
        if (as.character(class(ExpObj)) != "SummarizedExperiment"){
            stop("This function can only take a SummarizedExperiment object as input.")
        }

        #Let's see what we have here
        pheno_original <- colData(ExpObj)
        Samples_have <- rownames(pheno_original)
        if (!all(rownames(phenotable) %in% Samples_have)){
            stop("New phenotable contains samples not present in SummarizedExperiment. Phenotable sample names must match those in SummarizedExperiment. Check your new metadata and try again.")
        }

        #The SummarizedExperiment package does not permit coercing different metadata, so object has to be re-assembled from scratch
        assays_ExpObj <-(assays(ExpObj))
        ftt <- rowData(ExpObj)
        metadata_ExpObj <- metadata(ExpObj)

        ##Create SummarizedExperiment
        newExpObj <- SummarizedExperiment(assays = assays_ExpObj, rowData = ftt, colData = phenotable)
        metadata(newExpObj) <- metadata_ExpObj

        return(newExpObj)
}


#' fixrelpath(JAMSpath = NULL)
#'
#' Fixes path relativity and returns absolute path
#' @export

fixrelpath <- function(JAMSpath = NULL){
    require(R.utils)
    if (!(isAbsolutePath(JAMSpath))){
        fixedpath <- getAbsolutePath(JAMSpath)
    } else {
        fixedpath <- JAMSpath
    }

    return(fixedpath)
}

#' name_samples(list.data = NULL)
#'
#' Given a list.data object, returns a vector of sample names present in the object.
#' @export

name_samples <- function(list.data = NULL){
    loadedsamples <- gsub("_projinfo", "", (names(list.data)[grep("_projinfo", names(list.data))]))

    return(loadedsamples)
}
