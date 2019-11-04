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
