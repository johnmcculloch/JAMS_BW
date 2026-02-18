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

#' is.redundant(vec)
#'
#' Returns logical for a vector being redundant, i.e. has at least one element repeated in it.
#' @export
is.redundant <- function(vec){
    propunique <- length(unique(vec)) / length(vec)
    if (propunique < 1){
        redundant <- TRUE
    } else {
        redundant <- FALSE
    }

    return(redundant)
}


#' is.useful(vec)
#'
#' Returns logical for a vector having more than one class in it
#' @export

is.useful <- function(vec){
    numunique <- length(unique(vec))
    if (numunique < 2){
        useful <- FALSE
    } else {
        useful <- TRUE
    }

    return(useful)
}


#' convert_matrix_log2(mat = NULL, transformation = NULL)
#'
#' Wrapper for transforming a matrix to and from log2 formats
#' @export

convert_matrix_log2 <- function(mat = NULL, transformation = NULL){

    if (!(transformation %in% c("to_log2", "from_log2"))){
        stop("Please choose either \"to_log2\", or \"from_log2\" as transformation to apply to matrix.")
    }

    if (class(mat)[1] != "matrix"){
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


#' convert_matrix_PA(mat = NULL, threshPA = 0)
#'
#' Wrapper for transforming a matrix to and from log2 formats
#' @export

convert_matrix_PA <- function(mat = NULL, threshPA = 0){

    if (class(mat)[1] != "matrix"){
        stop("Object to transform must be a matrix.")
    }

    PAvec <- NULL
    for (rownum in 1:nrow(mat)){
        PAvec <- mat[rownum, ]
        if (threshPA != 0){
            PAvec[which(PAvec < threshPA)] <- 0
            PAvec[which(PAvec >= threshPA)] <- 1
        } else {
            PAvec[which(PAvec > 0)] <- 1
        }
        mat[rownum, ] <- PAvec
    }

    return(mat)

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


#' trim_whitespace_from_df(df)
#' Trims leading and trailing whitespace from a dataframe
#' @export

trim_whitespace_from_df <- function(df = NULL){
    #check if input is reasonable
    if (class(df)[1] != "data.frame"){
        stop("Input must be a data.frame")
    }

    for (colm in 1:ncol(df)){
        df[ , colm] <- trimws(df[ , colm])
    }

    return(df)
}

#' replace_NAs_with_character(df)
#' Replaces NAs with a specific character in a data frame
#' @export

replace_NAs_with_character <- function(df = NULL, replacement = "N_A"){
    #check if input is reasonable
    if (class(df)[1] != "data.frame"){
        stop("Input must be a data.frame")
    }

    for (colm in 1:ncol(df)){
        df[is.na(df[ , colm]) , colm] <- replacement
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

spew_heatmap_report <- function(hmcomb = NULL, outdir = NULL, expvec = NULL, applyfilters = NULL, variable_list = NULL, scaled = NULL, cdict = NULL, makespreadsheets = TRUE, makeheatmaps = TRUE, project = NULL){
    analysis <- hmcomb[1]
    report <- hmcomb[2]
    if (!file.exists(file.path(outdir, "Reports", analysis))){
        dir.create(file.path(outdir, "Reports", analysis), showWarnings = FALSE, recursive = TRUE)
    }
    setwd(file.path(outdir, "Reports", analysis))

    if (report == "comparative"){

        tryCatch((make_heatmap_report(report = "comparative", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, scaled = scaled, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", showonlypbelow = 0.05, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating comparative heatmaps for", analysis, "analysis")))

    } else if (report == "exploratory") {

        tryCatch((make_heatmap_report(report = "exploratory", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, scaled = scaled, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 2, numthreads = 1, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating exploratory heatmaps for", analysis, "analysis")))

    } else if (report == "correlation") {

        tryCatch((make_heatmap_report(report = "correlation", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, scaled = scaled, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", minabscorrcoeff = 0.55, ntopvar = 250)), error = function(e) whoopsieplot(paste("generating correlation heatmaps for", analysis, "analysis")))

    } else if (report == "PA") {

        tryCatch((make_heatmap_report(report = "PA", project = project, expvec = expvec, usefulexp = analysis, appendtofilename = analysis, applyfilters = applyfilters, variable_list = variable_list, scaled = scaled, cdict = cdict, makespreadsheets = TRUE, makeheatmaps = TRUE, maxnumheatmaps = 3, numthreads = 1, adjustpval = "auto", showonlypbelow = 0.05, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating presence/absence heatmaps for", analysis, "analysis")))

    } else if ((report %in% c("tUMAP", "PCA", "tSNE"))){

        tryCatch((make_ordination_report(algorithm = report, project = project, expvec = expvec, usefulexp = analysis, applyfilters = "none", appendtofilename = paste(analysis, report, sep = "_"), variable_list = variable_list, doreads = opt$doreads, cdict = cdict, threads = 1, class_to_ignore = opt$class_to_ignore)), error = function(e) whoopsieplot(paste("generating", report, "ordination plots for", analysis, "analysis")))

    }  else if (report == "alpha") {

        tryCatch(make_alpha_report(project = project, expvec = expvec, usefulexp = analysis, variable_list = variable_list, measures = c("Observed", "InvSimpson", "GeneCounts"), cdict = cdict, makespreadsheets = TRUE, stratify_by_kingdoms = TRUE, applyfilters = NULL, appendtofilename = paste(analysis, report, sep = "_"), GenomeCompletenessCutoff = c(5, 5), PPM_normalize_to_bases_sequenced = FALSE, max_pairwise_cats = 4, ignoreunclassified = TRUE, class_to_ignore = opt$class_to_ignore), error = function(e) whoopsieplot(paste("generating alpha diversity plots for", analysis, "analysis")))

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
        if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
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
    #Chomp a "/" from the end of paths
    fixedpath <- gsub("/$", "", fixedpath)

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


#' sample_by_category
#'
#' Sample by category
#' @export

sample_by_category <- function(rows, frac = 0.25) {
    #counts <- table(rows$category) %>% .[.>= n]
    counts <- table(rows$category)
    result <- data.frame()
    for (name in names(counts)) {
        result <- rbind(result, sample_frac(rows[rows$category==name,], frac))
    }

    return(result)
}


#' add_shape_to_plot_safely(p = NULL, shapevec = NULL, shapeby = NULL, cdict = NULL, use_letters_as_shapes = FALSE)
#'
#' #Given a vector of classes for adding shapes to a ggplot, attributes shapes safely in the presence or absence of a cdict containing shape info
#' @export

add_shape_to_plot_safely <- function (p = NULL, shapevec = NULL, shapeby = NULL, cdict = NULL, use_letters_as_shapes = FALSE){

    numshapes <- length(unique(shapevec))
    if (!use_letters_as_shapes & (numshapes <= 27)){
        shape_pecking_order <- c(19, 17, 15, 8, 12, 13, 18, 10, 3, 4, 11, 0, 1, 2, 5, 6, 7, 6, 35, 36, 38, 64, 163, 165, 167, 169, 198)
    } else {
        #Use letters
        shape_pecking_order <- c(65:90, 97:122, 35:38, 134:140)
    }

    if (!is.null(cdict)){
        st <- cdict[[shapeby]]
        if ("Shape" %in% colnames(st)){
            groupshapes <- setNames(as.numeric(st$Shape), as.character(st$Name))
            p <- p + scale_shape_manual(values = groupshapes)
        } else {
            p <- p + scale_shape_manual(values = shape_pecking_order[1:numshapes])
        }
    } else {
        p <- p + scale_shape_manual(values = shape_pecking_order[1:numshapes])
    }

    return(p)
}


#' find_clusters_in_matrix(input_matrix = NULL, kmeans_split = NULL)
#'
#' #Splits a matrix with counts in to k clusters by k-means clustering
#' @export

find_clusters_in_matrix <- function(input_matrix = NULL, kmeans_split = NULL){

    kmclusters <- stats::kmeans(input_matrix, centers = kmeans_split)
    clusterdf <- data.frame(Feature = names(fitted(kmclusters, method = "classes")), Cluster_Num = fitted(kmclusters, method = "classes"), stringsAsFactors = FALSE)
    clusterdf <- clusterdf[order(clusterdf$Cluster_Num), ]
    clusterdf$Cluster_Name <- paste("Cluster", as.character(sprintf("%02d", clusterdf$Cluster_Num)), sep = "_")

    return(clusterdf)

}


#' Author: Phillip Burger
#' Date: Sept 3, 2014
#' Purpose: Get the ordinal rank of a number.
#' http://www.phillipburger.net/wordpress/ordinal-number-suffix-function-in-r
#' @export

getOrdinalNumber1 <- function(num) {
    result <- ""
    if (!(num %% 100 %in% c(11, 12, 13))) {
        result <- switch(as.character(num %% 10),
            "1" = {paste0(num, "st")},
            "2" = {paste0(num, "nd")},
            "3" = {paste0(num, "rd")},
            paste0(num, "th"))
    } else {
        result <- paste0(num, "th")
    }
    result
}


#' install_biocon_deps()
#'
#' Installs Bioconductor JAMS dependencies
#' @export

install_biocon_deps <- function(){
bioconductordeps <- c("dada2", "ComplexHeatmap", "HybridMTest", "genefilter", "SummarizedExperiment")
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(bioconductordeps)
}


#' make_phenolabels_from_phenotable(phenotable = NULL)
#'
#' Does exactly what the name says.
#' @export

make_phenolabels_from_phenotable <- function (phenotable = NULL){
    Var_label <- colnames(phenotable)
    infer_column_type <- function(colm){
        if (colm %in% c("Sample", "sample")){
            colmtype <- "Sample"
        } else {
            if (can_be_made_numeric( (phenotable[, colm] ), cats_to_ignore = class_to_ignore)){
                colmtype <- "continuous"
            } else {
                colmtype <- "discrete"
            }
        }
        return(colmtype)
    }
    Var_type <- sapply(Var_label, function (x) { infer_column_type(x) } )
    phenolabels <- data.frame(Var_label = unname(Var_label), Var_type = unname(Var_type), stringsAsFactors = FALSE)

    return(phenolabels)

}

#' generate_filename(title = "ABC_look_at_me", add_date = TRUE, suffix = "pdf")
#'
#' Wrapper for generating file names.
#' @export

generate_filename <- function(title = "ABC_look_at_me", add_date = TRUE, suffix = "pdf"){

    if (add_date){
        dte <- Sys.Date()
    } else {
        dte <- NULL
    }

    fn <- paste(paste(title, dte, sep = "_"), suffix, sep = ".")

    return(fn)
}

#' analyze_fastq_filename(fn = NULL, maintain_Illumina_format = FALSE)
#'
#' Given a fastq filename, will return a list with facts about its anatomy.
#' @export

analyze_fastq_filename <- function(fn = NULL, maintain_Illumina_format = FALSE){

    filefacts <- list()
    fnsplit <- unlist(strsplit(fn, split = "\\."))
    #find out if gzipped and what the filename body is
    if (tail(fnsplit, n = 1) == "gz"){
        filefacts$gunzipped <- TRUE
        fnbody <- paste0(rev(rev(fnsplit)[3:length(fnsplit)]), collapse = ".")
        filefacts$suffix <- paste0(tail(fnsplit, n = 2), collapse = ".")
    } else {
        filefacts$gunzipped <- FALSE
        fnbody <- paste0(rev(rev(fnsplit)[2:length(fnsplit)]), collapse = ".")
        filefacts$suffix <- tail(fnsplit, n = 1)
    }

    #Now delve into fnbody and see what it it about
    #Does it have lane info? Must have L00X in the last 4 elements.
    contains_lane_info <- any(c("L001", "L002", "L003", "L004", "L005", "L006", "L007", "L008") %in% tail(unlist(strsplit(fnbody, split = "_")), n = 4))
    #Does it have paired read info? Must have RX in the last 4 elements.
    contains_paired_read_info <- any(c("R1", "R2") %in% tail(unlist(strsplit(fnbody, split = "_")), n = 4))
    #Is it Illumina output filename style?
    is_Illumina_output <- (tail(unlist(strsplit(fnbody, split = "_")), n = 1) == "001")
    #Get number of body segments
    numfnbodysegments <- length(unlist(strsplit(fnbody, split = "_")))

    #Determine the orgininal sample prefix
    #Case 1, not Illumina style
    if (!is_Illumina_output){
        if (contains_paired_read_info){
            filefacts$OriPrefix <- paste0(unlist(strsplit(fnbody, split = "_"))[1:(numfnbodysegments - 1)], collapse = "_")
            filefacts$Read <- tail(unlist(strsplit(fnbody, split = "_")), n = 1)
            filefacts$Appendage <-  paste(filefacts$Read, filefacts$suffix, sep = ".")
        } else {
            filefacts$OriPrefix <- fnbody
            filefacts$Read <- ""
            filefacts$Appendage <- filefacts$suffix
        }

    } else {
        #OK, is Illumina style prefix
        #Das Read ist immer an zweiter position
        filefacts$Read <- tail(unlist(strsplit(fnbody, split = "_")), n = 2)[1]

        if (contains_lane_info){
            #Case 2, is MiSeq, HiSeq or NextSeq file
            filefacts$OriPrefix <- paste0(unlist(strsplit(fnbody, split = "_"))[1:(numfnbodysegments - 4)], collapse = "_")
            filefacts$Lane <- tail(unlist(strsplit(fnbody, split = "_")), n = 3)[1]
            filefacts$SampNum <- tail(unlist(strsplit(fnbody, split = "_")), n = 4)[1]
            if (maintain_Illumina_format){
                filefacts$Appendage <- paste(paste(filefacts$SampNum, filefacts$Lane, filefacts$Read, "001", sep = "_"), filefacts$Container, sep = ".")
            } else {
                filefacts$Appendage <- paste(filefacts$Read, filefacts$suffix, sep = ".")
            }
        } else {
            #Case 3, is NovaSeq joined file, or contains no lane info
            filefacts$OriPrefix <- paste0(unlist(strsplit(fnbody, split = "_"))[1:(numfnbodysegments - 3)], collapse = "_")
            filefacts$Lane <- ""
            filefacts$SampNum <- tail(unlist(strsplit(fnbody, split = "_")), n = 3)[1]
            if (maintain_Illumina_format){
                filefacts$Appendage <- paste(paste(filefacts$SampNum, filefacts$Read, "001", sep = "_"), filefacts$suffix, sep = ".")
            } else {
                filefacts$Appendage <- paste(filefacts$Read, filefacts$suffix, sep = ".")
            }
        }
    }

    return(filefacts)

}


#' infer_column_type(phenotable = NULL, colm = NULL, class_to_ignore = "N_A")
#'
#' From a phenotable, infers what kind of a column it is, discrete or continuous.
#' @export

infer_column_type <- function(phenotable = NULL, colm = NULL, class_to_ignore = "N_A"){
    if (colm %in% c("Sample", "sample")){
        colmtype <- "Sample"
    } else {
        if (can_be_made_numeric( (phenotable[, colm] ), cats_to_ignore = class_to_ignore)){
            colmtype <- "continuous"
        } else {
            colmtype <- "discrete"
        }
    }

    return(colmtype)

}

#' compute_x_font(nx = 1, upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 8, cex = 0.7)
#'
#' Compute the appropriate size of font to be used in JAMS heatmap column (sample) annotations.
#' @export

compute_x_font <- function(nx = 1, upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 8, cex = 0.7){

    m <- ((lower_fs - upper_fs) / (upper_n - lower_n)) * -1
    cnst <- lower_fs - (lower_n * m)
    fontsizex <- round(((m * nx) + cnst), 2)
    fontsizex <- fontsizex * cex

    return(fontsizex)
}

#' hm_fontsize_computer(mat_rownames = NULL, mat_colnames = NULL, upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 8, cex = 0.7)
#'
#' Returns a vector of length 2 with values for the appropriate sizes of fonts to be used in JAMS heatmap columns (samples) and rows (features), respectively.
#' @export

hm_fontsize_computer <- function(mat_rownames = NULL, mat_colnames = NULL, upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 8, cex = 0.7){

    cexy <- cex
    #If not LKT and featname is split, then decrease row font size
    if (!(any(length(grep("^LKT_", mat_rownames) > 0)))){
        if (max(nchar(mat_rownames)) > 60){
            cexy <- cex * 0.7
        }
    }
    fontsizex <- compute_x_font(nx = length(mat_colnames), upper_n = upper_n, upper_fs = upper_fs, lower_n = lower_n, lower_fs = lower_fs, cex = cex)
    fontsizey <- compute_x_font(nx = length(mat_rownames), upper_n = upper_n, upper_fs = upper_fs, lower_n = lower_n, lower_fs = lower_fs, cex = cexy)
    #floor rownames to 0.1
    fontsizey <- max(0.1, fontsizey)
    xy_fs <- c(fontsizex, fontsizey)

    return(xy_fs)
}

#' split_featname(featname = NULL, thresh_featname_split = 40)
#'
#' For a feature name (or any string, for that matter) will add a carriage return to the middle of the string if it exceeds the threshold number of characters. This is used in the heatmap function.
#' @export

split_featname <- function(featname = NULL, thresh_featname_split = 40){
    if (nchar(featname) >= thresh_featname_split){
        len1 <- round((nchar(featname) / 2), 0)
        len2 <- nchar(featname) - len1
        part1 <- paste0(unlist(strsplit(featname, split = ""))[1:len1], collapse = "")
        part2 <- paste0(unlist(strsplit(featname, split = ""))[(1:len2) + len1], collapse = "")
        featname <- paste0(part1, "-\n", part2)
    }

    return(featname)
}


#' post_message(msg = NULL, face = "plain", just = "centre", color = "black", size = 20)
#'
#' Wrapper function for plotting text to a full page on a pdf.
#' @export

post_message <- function(msg = NULL, face = "plain", just = "centre", color = "black", size = 20){
    text_msg <- paste0(msg, collapse = "\n")
    tgrob <- text_grob(text_msg, face = face, color = color, just = just, size = size)
    print(as_ggplot(tgrob))
}


#' print_jams_ASCII(palette = "neon")
#'
#' Prints to terminal a retro ASCII JAMS logo.
#' @export

print_jams_ASCII <- function (palette = NULL){

    palettes <- list(
        #apple2 = c(34, 11, 202, 196, 127, 20),
        apple2rev = rev(c(34, 11, 202, 196, 127, 20)),
        neon = c(27, 93, 129, 201, 207, 99),
        sunset = c(202, 208, 214, 220, 226, 228),
        ocean = c(25, 31, 38, 44, 50, 56),
        fire = c(196, 202, 208, 214, 220, 226),
        grad_dark = c(16, 17, 18, 19, 20, 21),
        candy = c(21, 57, 93, 129, 165, 201),
        reds = 52:57,
        reds2 = (36 * 0:5) + 17,
        skies = (36 * 0:5) + 37,
        green_blue = 22:27,
        greys = 244:249,
        dark2 = c(16, 17, 18, 54, 53, 52),
        mono1 = rep(5, 6),
        mono2 = rep(4, 6),
        mono3 = rep(124, 6),
        mono4 = rep(246, 6)
     )

    if (is.null(palette)){
        #choose random palette
        palette <- sample(names(palettes), 1, replace = TRUE)
        rowvec <- palettes[[palette]]
    } else {
        rowvec <- palettes[[palette]]
    }

    ASCII_lines <- c(sprintf("\033[38;5;%dm         ██╗\033[38;5;%dm     █████╗\033[38;5;%dm     ███╗   ███╗\033[38;5;%dm     ███████╗\033[0m", rowvec[1], rowvec[1], rowvec[1], rowvec[1]), sprintf("\033[38;5;%dm        ██║\033[38;5;%dm    ██╔══██╗\033[38;5;%dm    ████╗ ████║\033[38;5;%dm    ██╔════╝\033[0m", rowvec[2], rowvec[2], rowvec[2], rowvec[2]), sprintf("\033[38;5;%dm       ██║\033[38;5;%dm    ███████║\033[38;5;%dm    ██╔████╔██║\033[38;5;%dm    ███████╗\033[0m", rowvec[3], rowvec[3], rowvec[3], rowvec[3]), sprintf("\033[38;5;%dm  ██   ██║\033[38;5;%dm    ██╔══██║\033[38;5;%dm    ██║╚██╔╝██║\033[38;5;%dm    ╚════██║\033[0m", rowvec[4], rowvec[4], rowvec[4], rowvec[4]),sprintf("\033[38;5;%dm ╚█████╔╝\033[38;5;%dm    ██║  ██║\033[38;5;%dm    ██║ ╚═╝ ██║\033[38;5;%dm    ███████║\033[0m", rowvec[5], rowvec[5], rowvec[5], rowvec[5]), sprintf("\033[38;5;%dm  ╚════╝ \033[38;5;%dm   ╚═╝  ╚═╝\033[38;5;%dm    ╚═╝     ╚═╝\033[38;5;%dm    ╚══════╝\033[0m", rowvec[6], rowvec[6], rowvec[6], rowvec[6]), paste("Microbial genomics analysis software.", paste("Version", packageVersion("JAMS"))))

    cat(paste0(ASCII_lines, collapse = "\n"), "\n")
}

#' extract_NCBI_taxid_from_featname(Taxon = NULL, NCBI_taxonomic_rank = NULL)
#'
#' Securely extracts the NCBI taxid from a JAMS2-style taxonomic feature name
#' @export

extract_NCBI_taxid_from_featname <- function(Taxon = NULL){

    #MetaBAT bin, so determine whether it is one which contains sample name (JAMSbeta style) or not (JAMSalpha style).
    #If LKT or cLKT, this will also work
    split_elements <- unlist(strsplit(Taxon, split = "_"))
    #Get earliest position matching a taxonomic tag
    taxtag_pos <- which(split_elements %in% c("d", "k", "p", "c", "o", "f", "g", "s", "is1"))[1]
    NCBI_taxid <- split_elements[(taxtag_pos + 2)]

    #If not found, then there is no tag, so convert NA to taxid 0 (root)
    if (is.na(NCBI_taxid)){
        NCBI_taxid <- "0"
    }
    #If is not numeric, like "missing", or "Unclassified", return "0"
    if (!can_be_made_numeric(NCBI_taxid)){
         NCBI_taxid <- "0"
    }

    return(NCBI_taxid)
}


#' check_for_banked_samples(JAMSarchive_path = NULL)
#'
#' Given an absolute path to a JAMSarchive style directory (wherein jamsfiles and JAMStarballs are), this function will return a vector containing the names of samples for which there are both a .jams file and a JAMS alpha run tarball. This may be useful to quickly ascertiain which samples have been banked or not, in the scope of project management.
#' @export

check_for_banked_samples <- function(JAMSarchive_path = NULL){
    jfl <- list.files(path = file.path(JAMSarchive_path, "jamsfiles"), pattern = ".jams")
    jtbl <- list.files(path = file.path(JAMSarchive_path, "JAMStarballs"), pattern = ".tar.gz")
    sh <- unname(sapply(jfl, function (x) { unlist(strsplit(x, split = "\\."))[1] } ))
    tbh <- gsub("_JAMS.tar.gz", "", jtbl)
    if (length(sh) != length(tbh)){
        flog.warn("Mismatch between .jams files and JAMStarballs!")
        jamfiles_without_tarballs <- sh[!(sh %in% tbh)]
        tarballs_without_jamsfiles <- tbh[!(tbh %in% sh)]
        if (length(jamfiles_without_tarballs) > 0){
            flog.warn(paste("The following .jams files have no tarballs:", paste0(jamfiles_without_tarballs, collapse = ", ")))
        }
        if (length(tarballs_without_jamsfiles) > 0){
            flog.warn(paste("The following tarballs have no .jams files:", paste0(tarballs_without_jamsfiles, collapse = ", ")))
        }
    }
    return(intersect(sh, tbh))
}

#' retrieve_GCA_accession_URL(assemblytargetname = NULL)
#'
#' Given a single valid GCA assembly accession number, this function will return the https URL for downloading the fasta.gz file for that assembly accession.
#' @export

retrieve_GCA_accession_URL <- function(assemblytargetname){

    #First, check you have esearch
    missingdep <- FALSE
    for (dep in c("esearch", "esummary", "xtract")) {
        if (system(str_c("which ", dep), ignore.stdout = TRUE) ==1 ) { 
            flog.warn(str_c("You are missing ", dep))
            missingdep <- TRUE
        }
    }

    if (!missingdep){
        esearchcmd <- paste("esearch -db assembly -query", assemblytargetname, "| esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank >", paste(assemblytargetname, "url", sep = "."))
        system(esearchcmd)
        assemblyURL <- readLines(paste(assemblytargetname, "url", sep = "."), n = 1L, warn = FALSE)
        assemblyURL <- sub("^ftp:", "https:", assemblyURL)
        assemblyURL <- paste(paste(assemblyURL, basename(assemblyURL), sep = "/"), "genomic.fna.gz", sep = "_")

        #if GCA accession is spoof, return NULL
        if (assemblyURL == "_genomic.fna.gz"){
            assemblyURL <- NULL
        }

        return(assemblyURL)

    } else {
        #You don't have what it takes to obtain the accession URL
        return(NULL)
    }
}

#' pretty_time(seconds)
#'
#' Automatically reformats numbers of seconds to a more human-readable format.
#' @export

pretty_time <- function(seconds) {
    if (seconds < 60) {
        sprintf("%.1f s", seconds)
    } else if (seconds < 3600) {
        sprintf("%.1f min", seconds / 60)
    } else {
        sprintf("%.1f h", seconds / 3600)
    }
}


#' tally_metadata(md = NULL, column_to_group_by = NULL, columns_to_tally = NULL)
#'
#' Tallies the number of classes present in a metadata df and returns a data frame with counts.
#' @export

tally_metadata <- function(md = NULL, column_to_group_by = NULL, columns_to_tally = NULL){

    tally_column <- function(md = NULL, column_to_group_by = NULL, col2tally = NULL){ 
        tmp_tally_df <- md %>% dplyr::count(across(all_of(column_to_group_by)),across(all_of(col2tally))) %>% tidyr::pivot_wider(names_from = all_of(col2tally),values_from = n, values_fill = 0) %>% dplyr::rename_with( ~ paste0("Num_Samples_in_", col2tally, ".", .x), -all_of(column_to_group_by)) %>% as.data.frame()

        return(tmp_tally_df)
    }
 
    tallied_column_list <- lapply(columns_to_tally, function (x) { tally_column(md = md, column_to_group_by = column_to_group_by, col2tally = x) })
   
    #Stitch all together into a single data frame
    md_tally_df <- purrr::reduce(tallied_column_list, dplyr::left_join, by = column_to_group_by)

    return(md_tally_df)
}



#' plot_table_a4_autofit_paged(dataframe = NULL, title = NULL, orientation = NULL, rows_per_page = 35, start_fontsize = 10, min_fontsize = 5, margins_portrait_mm = c(top=6, right=6, bottom=6, left=6), margins_landscape_mm = c(top=6, right=8, bottom=6, left=8), gap_mm = 2, cell_padding_mm = c(v = 0.6, h = 1.2), header_padding_mm = c(v = 0.8, h = 1.2), lineheight = 0.95, wrap = FALSE, wrap_chars = NULL, wrap_min = 12, wrap_max = 60)
#' EXPERIMENTAL novel plotting function for printing data frames into a pdf.
#' Usage:
#' pdf("my_landscape.pdf", height = 210/25.4, width = 297/25.4)  # A4 landscape
#' plot_table_a4_autofit_paged(dataframe = opt$featuredata[1:80, ], title="Featuretable", rows_per_page = 35, wrap = TRUE, wrap_chars = 20, start_fontsize = 15, min_fontsize = 2,  orientation = "landscape")
#' dev.off()
#'
#' @export

plot_table_a4_autofit_paged <- function(dataframe = NULL, title = NULL, orientation = "landscape", rows_per_page = 35, start_fontsize = 10, min_fontsize = 5, margins_portrait_mm = c(top=6, right=6, bottom=6, left=6), margins_landscape_mm = c(top=6, right=8, bottom=6, left=8), gap_mm = 2, cell_padding_mm = c(v = 0.6, h = 1.2), header_padding_mm = c(v = 0.8, h = 1.2), lineheight = 0.95, wrap = FALSE, wrap_chars = NULL, wrap_min = 12, wrap_max = 60) {

    if (!is.data.frame(dataframe)) {
        stop("`dataframe` must be a data.frame")
    }
    orientation <- match.arg(orientation)

    # define margins based on orientation
    margins_mm <- if (orientation == "portrait") margins_portrait_mm else margins_landscape_mm

    if (length(margins_mm) == 1) {
        margins_mm <- c(top=margins_mm, right=margins_mm, bottom=margins_mm, left=margins_mm)
    } else if (is.null(names(margins_mm))) {
        stop("margins_*_mm must be named (top/right/bottom/left) or a single number")
    }

    # optional wrapping (applied once to full df so it stays consistent across pages)
    df <- dataframe
  if (wrap) {
    if (is.null(wrap_chars)) {
      nc <- max(1, ncol(df))
      base <- round(90 / nc)
      bump <- round((start_fontsize - 8) * 2)
      wrap_chars <- max(wrap_min, min(wrap_max, base + bump))
    }
    wrap_one <- function(x) {
      x <- as.character(x)
      vapply(x, function(s) paste(strwrap(s, width = wrap_chars), collapse = "\n"), character(1))
    }
    df[] <- lapply(df, wrap_one)
  }

  # device size in mm (device must be open)
  dev_in <- dev.size("in")
  page_w_mm <- dev_in[1] * 25.4
  page_h_mm <- dev_in[2] * 25.4

  avail_w_mm <- page_w_mm - (margins_mm["left"] + margins_mm["right"])
  avail_h_mm <- page_h_mm - (margins_mm["top"] + margins_mm["bottom"])
  if (avail_w_mm <= 0 || avail_h_mm <= 0) stop("Margins too large for device size")

  # build title grob once; show it on every page by default
  title_g <- NULL
  title_h_mm <- 0
  if (!is.null(title) && nzchar(title)) {
    title_g <- textGrob(title, gp = gpar(fontsize = max(10, start_fontsize + 2)))
    title_h_mm <- convertHeight(grobHeight(title_g), "mm", valueOnly = TRUE)
  }

  make_theme <- function(fs) {
    ttheme_default(
      base_size = fs,
      core = list(
        padding = unit(c(cell_padding_mm["v"], cell_padding_mm["h"]), "mm"),
        fg_params = list(lineheight = lineheight)
      ),
      colhead = list(
        padding = unit(c(header_padding_mm["v"], header_padding_mm["h"]), "mm"),
        fg_params = list(lineheight = lineheight)
      )
    )
  }

  measure_tbl_mm <- function(tbl) {
    w <- convertWidth(sum(tbl$widths), "mm", valueOnly = TRUE)
    h <- convertHeight(sum(tbl$heights), "mm", valueOnly = TRUE)
    c(width = as.numeric(w), height = as.numeric(h))
  }

  # split into pages (35 rows per page if needed)
  n <- nrow(df)
  if (n == 0) {
    warning("Empty data frame: nothing to plot.")
    return(invisible(NULL))
  }

  starts <- seq(1, n, by = rows_per_page)
  ends <- pmin(starts + rows_per_page - 1, n)

  # For each page, fit fontsize (can vary per page, but we’ll keep it consistent:
  # choose the *smallest* fontsize required across all pages so pages match visually.)
  fitted_fontsizes <- integer(length(starts))

  # First pass: compute needed fontsize per chunk
  for (i in seq_along(starts)) {
    chunk <- df[starts[i]:ends[i], , drop = FALSE]
    fs <- start_fontsize
    repeat {
      tbl <- tableGrob(chunk, rows = NULL, theme = make_theme(fs))
      dims <- measure_tbl_mm(tbl)

      avail_table_h_mm <- avail_h_mm - title_h_mm - ifelse(is.null(title_g), 0, gap_mm)
      if (avail_table_h_mm <= 0) avail_table_h_mm <- avail_h_mm - title_h_mm

      fits <- (dims["width"] <= avail_w_mm) && (dims["height"] <= avail_table_h_mm)

      if (fits || fs <= min_fontsize) break
      fs <- fs - 1
    }
    fitted_fontsizes[i] <- fs
  }

  # Use a single fontsize across pages for consistency (smallest that worked)
  final_fs <- min(fitted_fontsizes)

  # Second pass: draw pages
  for (i in seq_along(starts)) {
    chunk <- df[starts[i]:ends[i], , drop = FALSE]
    tbl <- tableGrob(chunk, rows = NULL, theme = make_theme(final_fs))

    grid.newpage()

    # margins in npc
    left_npc   <- margins_mm["left"]   / page_w_mm
    right_npc  <- margins_mm["right"]  / page_w_mm
    top_npc    <- margins_mm["top"]    / page_h_mm
    bottom_npc <- margins_mm["bottom"] / page_h_mm

    avail_w_npc <- 1 - left_npc - right_npc

    # title
    title_h_npc <- 0
    gap_npc <- 0
    if (!is.null(title_g)) {
      title_h_npc <- convertHeight(grobHeight(title_g), "npc", valueOnly = TRUE)
      gap_npc <- gap_mm / page_h_mm

      pushViewport(viewport(
        x = left_npc + avail_w_npc / 2,
        y = 1 - top_npc - title_h_npc / 2,
        width = avail_w_npc,
        height = title_h_npc,
        just = c("center", "center")
      ))
      grid.draw(title_g)
      popViewport()
    }

    # table area below title
    table_top <- 1 - top_npc - title_h_npc - gap_npc
    table_bottom <- bottom_npc
    table_h <- table_top - table_bottom
    table_center_y <- table_bottom + table_h / 2

    pushViewport(viewport(
      x = left_npc + avail_w_npc / 2,
      y = table_center_y,
      width = avail_w_npc,
      height = table_h,
      just = c("center", "center"),
      clip = "on"
    ))
    grid.draw(tbl)
    popViewport()
  }

  invisible(list(
    pages = length(starts),
    rows_per_page = rows_per_page,
    final_fontsize = final_fs,
    wrap_chars = if (wrap) wrap_chars else NA
  ))
}

