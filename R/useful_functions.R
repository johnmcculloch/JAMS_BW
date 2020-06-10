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

#' retrieve_features_by_taxa(FuncExpObj = NULL, wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, PPMthreshold = 0)
#'
#' Returns a long form data frame of stratification by taxa of the relative abundance or number of bases wanted of functional features in wanted samples, given allfeaturesbytaxa_matrix and allfeaturesbytaxa_index metadata present in a JAMS SummarizedExperiment functional object.
#' @export

retrieve_features_by_taxa <- function(FuncExpObj = NULL, wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, PPMthreshold = 0){

    allfeaturesbytaxa_matrix <- metadata(FuncExpObj)$allfeaturesbytaxa_matrix
    allfeaturesbytaxa_index <- metadata(FuncExpObj)$allfeaturesbytaxa_index
    curr_pt <- colData(FuncExpObj)

    #Get appropriate rows
    rowsinterestdf <- subset(allfeaturesbytaxa_index, Sample %in% wantedsamples)
    rowsinterestdf <- subset(rowsinterestdf, Accession %in% wantedfeatures)

    allfeaturesbytaxa_interest <- as.matrix(allfeaturesbytaxa_matrix[rowsinterestdf$RowNumber, ])

    if (length(rowsinterestdf$RowNumber) == 1){
        allfeaturesbytaxa_interest <- t(allfeaturesbytaxa_interest)
    }

    allfeaturesbytaxa_interest <- as.data.frame(allfeaturesbytaxa_interest[, which(colSums(allfeaturesbytaxa_interest) != 0)])
    allfeaturesbytaxa_interest$RowNumber <- as.numeric(rownames(allfeaturesbytaxa_interest))

    allfeaturesbytaxa_interest <- left_join(allfeaturesbytaxa_interest, rowsinterestdf, by = "RowNumber")
    allfeaturesbytaxa_interest$RowNumber <- NULL
    allfeaturesbytaxa_interest <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", (sort(colnames(allfeaturesbytaxa_interest)[which(!colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession"))])))]

    if (asPPM){
        taxsplit <- allfeaturesbytaxa_interest
        #Transform to PPM
        if (PPM_normalize_to_bases_sequenced == TRUE){
            totbases <- "TotalBasesSequenced"
        } else {
            totbases <- "TotalBasesSequencedinAnalysis"
        }
        numbases2sampl <- as.data.frame(t(metadata(FuncExpObj)[[totbases]]))
        numbases2sampl$Sample <- rownames(numbases2sampl)
        taxsplit <- left_join(taxsplit, numbases2sampl, by = "Sample")
        LKTcolumns <- colnames(taxsplit)[!(colnames(taxsplit) %in% c("Sample", "Accession", "NumBases"))]

        #Transform to PPM
        for(colm in LKTcolumns){
            taxsplit[ , colm] <- round(((taxsplit[ , colm] / taxsplit$NumBases) * 1000000), 0)
        }

        #Denoise
        LKTsMaxima <-sapply(LKTcolumns, function(x) {max(taxsplit[ , x])})
        LKTsToKeep <- names(which(LKTsMaxima > PPMthreshold))
        sample2metadata <- as.data.frame(curr_pt)
        nonSamplecolms <- colnames(sample2metadata)[colnames(sample2metadata) != "Sample"]
        taxsplit <- left_join(taxsplit, sample2metadata, by = "Sample")
        taxsplit <- taxsplit[ , c("Sample", "Accession", "NumBases", nonSamplecolms, LKTsToKeep)]
        taxsplit$NumBases <- NULL

        return(taxsplit)

    } else {

        return(allfeaturesbytaxa_interest)

    }
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


#' shrink_perbasecoverage(perbasecoverage=NULL, percentage=2)
#'
#' #Subset down to n% of the dataset to make it calculatable. To be used only by JAMSalpha.
#' @export

shrink_perbasecoverage <-function(perbasecoverage = NULL, percentage = 2){
    numbases <- 1:nrow(perbasecoverage)
    breaks <- numbases[seq(1, length(numbases), (100/percentage))]
    perbasecoverage_reduced <- perbasecoverage[breaks, ]

    return(perbasecoverage_reduced)
}

#' add_shape_to_plot_safely(p = NULL, shapevec = NULL, cdict = NULL)
#'
#' #Given a vector of classes for adding shapes to a ggplot, attributes shapes safely in the presence or absence of a cdict containing shape info
#' @export

add_shape_to_plot_safely <- function (p = NULL, shapevec = NULL, shapeby = NULL, cdict = NULL){

    shape_pecking_order <- c(19, 17, 15, 8, 12, 13, 18, 10, 3, 4, 11, 0, 1, 2, 5, 6, 7, 6, 35, 36, 38, 64)
    numshapes <- length(unique(shapevec))

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
