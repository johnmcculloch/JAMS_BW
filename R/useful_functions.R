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
    analysis = hmcomb[1]
    report = hmcomb[2]
    setwd(file.path(opt$outdir, "Reports", analysis))

    if (report == "comparative"){

        tryCatch((make_heatmap_report(report = "comparative", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, list.data = list.data, mgSeqnorm = FALSE, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", showonlypbelow = 0.05, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating comparative heatmaps for", analysis, "analysis")))

    } else if (report == "exploratory") {

        tryCatch((make_heatmap_report(report = "exploratory", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, list.data = list.data, mgSeqnorm = FALSE, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 2, numthreads = 1, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating exploratory heatmaps for", analysis, "analysis")))

    } else if (report == "correlation") {

        tryCatch((make_heatmap_report(report = "correlation", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, list.data = list.data, mgSeqnorm = FALSE, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", minabscorrcoeff = 0.55, ntopvar = 250)), error = function(e) whoopsieplot(paste("generating correlation heatmaps for", analysis, "analysis")))

    } else if (report == "PA") {

        tryCatch((make_heatmap_report(report = "PA", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, list.data = list.data, mgSeqnorm = FALSE, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", showonlypbelow = 0.05, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating presence/absence heatmaps for", analysis, "analysis")))

    } else if ((report %in% c("tUMAP", "PCA", "tSNE"))){

        tryCatch((make_ordination_report(algorithm = report, project = project, expvec = expvec, usefulexp = analysis, appendtofilename = paste(analysis, report, sep = "_"), mgSeqnorm = FALSE, variable_list = variable_list, list.data = list.data, doreads = opt$doreads, cdict = cdict, threads = 1, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating", report, "ordination plots for", analysis, "analysis")))

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
