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


 #' getN50(contigs2length)
 #' JAMSalpha function to get N50 from contigs2length
 #'
 #' @export
 getN50 <- function(contigs2length){
     quantile50 <- (sum(as.numeric(contigs2length$Length))) * 0.5
     N50 <- max(subset(contigs2length, Lcum > quantile50)$Length)
     return(N50)
 }

 #' getL50(contigs2length)
 #' JAMSalpha function to get L50 from contigs2length
 #'
 #' @export
 getL50 <- function(contigs2length){
     quantile50 <- (sum(as.numeric(contigs2length$Length))) * 0.5
     L50 <- length(subset(contigs2length, Lcum > quantile50)$Length)
     return(L50)
 }

 #' getN90(contigs2length)
 #' JAMSalpha function to get N90 from contigs2length
 #'
 #' @export
 getN90 <- function(contigs2length){
     quantile90 <- (sum(as.numeric(contigs2length$Length))) * 0.9
     N90 <- max(subset(contigs2length, Lcum > quantile90)$Length)
     return(N90)
 }


 #' estimate_genome_completeness(contigs2length)
 #' JAMSalpha function to estimate genome completeness from contigs2length
 #'
 #' @export
 estimate_genome_completeness <- function(contigs2length = NULL){
     taxinteresttype <- gsub("^LKT__", "", unique(contigs2length$Taxon))
     expgensize <- as.integer(JAMSMedian_Genome_Sizes[taxinteresttype, "Genome_Size_Median"])
     if ((length(expgensize) < 1) | (is.na(expgensize)) | (taxinteresttype %in% c("Unclassified", "Missing", "Unclassified__Unclassified"))){
         expgensize <- as.integer(JAMSMedian_Genome_Sizes["k__Bacteria", "Genome_Size_Median"])
     }
     probnumgen <- round((sum(contigs2length$Length) / expgensize), 2)
     return(probnumgen)
 }

 #' find16SrRNA(contigs2length = NULL, opt = NULL)
 #' JAMSalpha function to find number of 16S genes from contigs2length
 #'
 #' @export
 find16SrRNA <- function(contigs2length = NULL, opt = NULL){
     featuredata = opt$featuredata
     #Find how many 16S rRNAs which have been classified can be found in the contigs were looking at.
     num16S <- length(which((subset(featuredata, Contig %in% contigs2length$Contig)[]$SixteenSid) != "none"))

     return(num16S)
 }


 #' getcontigs2length(opt = NULL, taxlevel = NULL, taxinterest = NULL)
 #' JAMSalpha function to compute total contig length of a taxon from contigs2length
 #'
 #' @export
 getcontigs2length <- function(opt = NULL, taxlevel = NULL, taxinterest = NULL){
     contigsdata <- opt$contigsdata
     contigs2length <- subset(contigsdata, contigsdata[,taxlevel] == taxinterest)
     contigs2length <- contigs2length[, c(taxlevel, "Contig", "Length")]
     colnames(contigs2length) <- c("Taxon", "Contig", "Length")
     contigs2length <- contigs2length[order(contigs2length$Length, decreasing = TRUE), ]
     contigs2length$Lcum <- cumsum(contigs2length$Length)
     return(contigs2length)
 }


 #' get_assembly_stats_by_taxon(opt = NULL, taxlevel = NULL, taxinterest = NULL)
 #' JAMSalpha function to compute assembly statistics by taxon
 #'
 #' @export
 get_assembly_stats_by_taxon <- function(opt = NULL, taxlevel = NULL, taxinterest = NULL){
     contigs2length <- getcontigs2length(opt=opt, taxlevel=taxlevel, taxinterest = taxinterest)
     assemblystats <- data.frame(TaxLevel = taxlevel, Taxon = taxinterest, NumContigs = nrow(contigs2length), ContigSum = sum(contigs2length$Length), LargestContigSize = max(contigs2length$Length), N50 = getN50(contigs2length), L50 = getL50(contigs2length), N90 = getN90(contigs2length), ProbNumGenomes = estimate_genome_completeness(contigs2length), Num16S = find16SrRNA(contigs2length = contigs2length, opt = opt))
     return(assemblystats)
 }

 #' get_assembly_stats_by_taxlevel(opt = NULL, taxlevel = NULL)
 #' JAMSalpha function to compute assembly statistics by taxonomic level
 #'
 #' @export
 get_assembly_stats_by_taxlevel <- function(opt = NULL, taxlevel = NULL){
     contigsdata <- opt$contigsdata
     #Get a list of non-redundant tax entities
     Taxa <- unique(contigsdata[ , taxlevel])
     assemblystats_taxlevel <- lapply(Taxa, function(t) { get_assembly_stats_by_taxon(opt = opt, taxlevel = taxlevel, taxinterest = t) })
     assemblystats_taxlevel <- plyr::ldply(assemblystats_taxlevel, rbind)
     assemblystats_taxlevel <- assemblystats_taxlevel[order(assemblystats_taxlevel$ProbNumGenomes, decreasing = TRUE), ]
     return(assemblystats_taxlevel)
 }


 #' IO_jams_workspace_image(opt = NULL, workspaceimage = NULL, threads = 8, operation = c("save", "load"))
 #' Safe way of either loading or saving an R workspace image. If argument workspaceimage is null, workspace image file will be searched for in opt (opt$projimage). If that is also NULL, saving or loading is aborted. If the fastSave package () is installed, multi-threaded loading or saving will be used. If opt is passed, number of CPUs will be set to opt$threads, trumping the threads argument.
 #'
 #' @export
IO_jams_workspace_image <- function(opt = NULL, workspaceimage = NULL, threads = 8, operation = NULL){
    if (is.null(workspaceimage)){
        workspaceimage <- opt$projimage
    }

    if (is.null(workspaceimage)){
        stop("Workspace image not found or specified.")
    }

    if (!is.null(opt)){
        threads <- opt$threads
    }

    flog.info(paste("Workspace image is", workspaceimage))

    if (operation == "save"){

        if ("fastSave" %in% rownames(installed.packages())){
            flog.info(paste("Saving project workspace image using fastSave package with", opt$threads, "CPUs"))
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
}
