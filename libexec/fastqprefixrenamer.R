#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(library(optparse)))
suppressWarnings(suppressPackageStartupMessages(library(futile.logger)))
suppressWarnings(suppressPackageStartupMessages(library(benchmarkme)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(RCurl)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(openxlsx)))

#####################################
# Define System-specific Functions ##
#####################################
if ((.Platform$OS.type) != "unix"){
    stop("This script only works on UNIX. Install Linux and try again.")
}

#Decide which kind of system you are on.
detectHardwareResources <- function(){
    #First off, detect if on Slurm type cluster
    #Get slurm job ID
    slurmjobid <- as.character(Sys.getenv("SLURM_JOB_ID"))

    if(nchar(slurmjobid) < 3){
       #Define appropriate functions for non-slurm system
       detectBatchCPUs <- function() {
            ncores <- detectCores()
            if (is.na(ncores)) {
                stop("Could not determine how many CPUs you have. Aborting.")
            }
            return(ncores)
        }

        detectAvailRAM <- function(){
            totmembytes <- as.numeric(get_ram())

            return(totmembytes)
        }

    } else {
        #Define appropriate functions for slurm system
        detectBatchCPUs <- function() {
            ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
            if (is.na(ncores)) {
                ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
            }
            if (is.na(ncores)) {
                stop("Could not determine how many CPUs you have. Aborting.")
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

            totmembytes<-totmem * 1000000

            return(totmembytes)
        }
    }
    hardwareRes <- NULL
    hardwareRes[1] <- detectBatchCPUs()
    hardwareRes[2] <- detectAvailRAM()
    names(hardwareRes) <- c("threads", "memory")

    return(hardwareRes)
}

############################
## Define other functions ##
############################
filetype <- function(path){
    f = file(path)
    ext = summary(f)$class
    close.connection(f)
    ext
}

# get path of running script
getScriptPath <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
        return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
    } else {
        return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
}

#########################
# Get options from args #
#########################
#Define defaults
defopt <- list()
defopt$verstr <- paste0("fastqprefixrenamer v", "1.0")
defopt$readsfolder <- getwd()
defopt$destination <- getwd()
defopt$basespaceformat <- FALSE
defopt$move <- FALSE
defopt$simulate <- FALSE
defopt$threads <- as.numeric(detectHardwareResources()["threads"])

option_list <- list(
    make_option(c("-r", "--readsfolder"), default = defopt$readsfolder, action="store",
                help = str_c("path to root directory where reads are (default: ", defopt$readsfolder, ")")),

    make_option(c("-d", "--destination"), default = defopt$destination, action="store",
                help = str_c("path to root directory where reads are (default: ", defopt$destination, ")")),

    make_option(c("-t", "--substitutionmap"), default = NULL, action="store",
                help = str_c("Tab limited text file for sample name substitution in the format: oldname[TAB]newname")),

    make_option(c("-x", "--substitutionspreadsheet"), default = NULL, action="store",
                help = str_c("Tab limited xlsx spreadsheet for sample name substitution")),

    make_option(c("-b", "--basespaceformat"), default = defopt$basespaceformat,
                action="store_true", help = str_c("Original prefixes are in Basespace format, i.e. OldPrefix_S264_R1_001.fastq.gz. (default: FALSE. Prefixes are in OldPrefix_R1.fastq format.)")),

    make_option(c("-m", "--move"), default = defopt$move, action="store_true",
                help = str_c("Copy to destination rather than move. (default: FALSE. Prefixes are in OldPrefix_R1.fastq format.)")),

    make_option(c("-c", "--checkprefixexists"), default = NULL,
                action="store", help = str_c("Optional path to folder with banked .jams files. If prefix already exists, a warning is issued, will abort before start of renaming.")),

    make_option(c("-s", "--simulate"), default = defopt$simulate, action="store_true",
                help = str_c("Simulate and print commands for renaming, do not actually rename.")),

    make_option(c("-v", "--version"), action="store_true",
                help ="report version")
)

# parse the options
args <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(OptionParser(option_list = option_list), args)
opt <- merge.list(opt, defopt)

#####################
## Set environment ##
#####################

# print version & exit if -v
if (!is.null(opt$version)) {
    print(opt$verstr)
    quit()
}

#Get Script path
opt$bindir <- getScriptPath()

#Fix path relativity
fixrelpath <- function(JAMSpath = NULL){
    require(R.utils)
    if (!(isAbsolutePath(JAMSpath))){
        fixedpath <- getAbsolutePath(JAMSpath)
    } else {
        fixedpath <- JAMSpath
    }

    return(fixedpath)
}

for (pathtofix in c("readsfolder", "destination", "substitutionmap", "substitutionspreadsheet", "checkprefixexists")){
    if (!is.null(opt[[pathtofix]])){
        opt[[pathtofix]] <- fixrelpath(opt[[pathtofix]])
    }
}

# give help if needed input option not provided
if (is.null(opt$readsfolder)) {
    print("You must supply a folder with reads to rename.")
    parse_args(OptionParser(option_list = option_list), c("-h"))
    q()
}


#Define useful functions
is.redundant <- function(vec){
    propunique <- length(unique(vec))/length(vec)
    if (propunique < 1){
        redundant = TRUE
    } else {
        redundant = FALSE
    }

    return(redundant)
}

find_container <- function(x) {
    prefsuff <- unlist(strsplit(rawfastqsdf$Filenames[x], split = "\\."))
    container <- paste(prefsuff[(which(prefsuff == "fastq")):(length(prefsuff))], collapse = ".")

    return(container)
}

###################
## Main Function ##
###################

if (opt$destination != opt$readsfolder){
    flog.info("Destination directory and origin are different.")
    if (!(file.exists(opt$destination))){
        flog.info("Creating directory to hold renamed output reads.")
        dir.create(opt$destination, recursive = TRUE)
    }
}

setwd(opt$destination)
opt$projimage <- file.path(opt$destination, ".RData")
save.image(opt$projimage)

#Find fastq files at the origin
rawfastqpaths <- list.files(path = opt$readsfolder, pattern = "fastq", full.names = TRUE, recursive = TRUE)
rawfastqnames <- sapply(1:length(rawfastqpaths), function(x) { tail(unlist(strsplit(rawfastqpaths[x], split = "/")), n = 1) })
datasetids <- sapply(1:length(rawfastqpaths), function(x) { tail(unlist(strsplit(rawfastqpaths[x], split = "/")), n = 2)[1] })

#Make a dataframe with all the info
rawfastqsdf <- data.frame(Filenames = rawfastqnames, Filepaths = rawfastqpaths, DatasetIDs = datasetids, stringsAsFactors = FALSE)

#Find out read F and R
if (opt$basespaceformat == TRUE){
    rawfastqsdf$Read <- sapply(1:nrow(rawfastqsdf), function(x) { tail(unlist(strsplit(rawfastqsdf$Filenames[x], split = "_")), n = 2)[1] } )
    rawfastqsdf$OriPrefix <- sapply(1:nrow(rawfastqsdf), function(x) { paste(unlist(strsplit(rawfastqsdf$Filenames[x], split = "_"))[1:((which(unlist(strsplit(rawfastqsdf$Filenames[x], split = "_")) == tail(unlist(strsplit(rawfastqsdf$Filenames[x], split = "_")), n = 4)[1])) - 1)], collapse = "_") } )
} else {
    flog.info("Still working on other input formats. Check for future versions.")
    q()
}

save.image(opt$projimage)

#Load substitution map
#Adjust prefix if there is a substitution map
if (!(is.null(opt$substitutionmap))){
    submap <- read.table(opt$substitutionmap, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(submap) <- c("OriPrefix", "Prefix")
} else if(!(is.null(opt$substitutionspreadsheet))){
    submap <- read.xlsx(opt$substitutionspreadsheet)
    colnames(submap) <- c("OriPrefix", "Prefix")
} else {
    flog.info("If you want to rename samples, you must provide a substitution map. Aborting now.")
    q()
}

#Test for redundancy in substitution map.
redtestsubmap <- apply(submap, is.redundant, MARGIN=2)
if (any(redtestsubmap)){
    flog.warn(paste("Values in substitution table column(s)" , paste0(names(which(redtestsubmap == TRUE)), collapse = " and "), "are not unique. Check table and try again. Aborting now"))
    q()
}

rawfastqsdf <- subset(rawfastqsdf, OriPrefix %in% (submap$OriPrefix))

if (nrow(rawfastqsdf) < 1){
    flog.info("No fastq files bearing the original prefixes in the substitution table provided were found. Aborting now.")
    q()
}

#Tack on the new prefixes
rawfastqsdf <- left_join(rawfastqsdf, submap, by = "OriPrefix")
#Tack on the destinations
rawfastqsdf$Container <- sapply(1:nrow(rawfastqsdf), function (x) { find_container(x) } )
rawfastqsdf$NewFN <- paste(paste(rawfastqsdf$Prefix, rawfastqsdf$Read, sep = "_"), rawfastqsdf$Container, sep = ".")


#Check if they are duplicates
rawfastqsdfR1 <- subset(rawfastqsdf, Read == "R1")
if (is.redundant(rawfastqsdfR1$OriPrefix)){
    flog.warn("There is more than one input fastq file with the same prefix. Check input reads folder substitution table and try again. Aborting now.")
    dupes <- rawfastqsdfR1$OriPrefix[duplicated(rawfastqsdfR1$OriPrefix)]
    flog.warn(paste("The following original prefixes found in the reads folder are duplicated:", paste0(dupes, collapse = ", ")))
    dupefns <- rawfastqsdf$Filepaths[which(rawfastqsdf$OriPrefix %in% dupes)]
    flog.warn(paste("The following input filenames have the same prefix:\n", paste0(dupefns, collapse = "\n")))
    q()
}

if (is.redundant(rawfastqsdfR1$Prefix)){
    flog.warn("There is more than one output fastq file with the same prefix. Check input reads folder substitution table and try again. Aborting now.")
    dupes <- rawfastqsdfR1$Prefix[duplicated(rawfastqsdfR1$Prefix)]
    flog.warn(paste("The following new prefixes are duplicated:", paste0(dupes, collapse = ", ")))
    dupefns <- rawfastqsdf$NewFN[which(rawfastqsdf$Prefix %in% dupes)]
    flog.warn(paste("The following proposed output filenames are identical:", paste0(dupefns, collapse = ", ")))
    q()
}

#Check if new prefix has already been used.
if (!is.null(opt$checkprefixexists)){
    flog.info(paste("Checking if new prefixes proposed in substitution table coincide with prefixes found in:", as.character(opt$checkprefixexists)))
    jamsfns <- list.files(path = opt$checkprefixexists, pattern=".jams$")
    usedprefixes <- sapply(1:length(jamsfns), function(x) { unlist(strsplit(jamsfns[x], split="\\."))[1] } )
}

prefixestouse <- unique(rawfastqsdf$Prefix)
if (!is.null(opt$checkprefixexists)){
    if (length(usedprefixes) > 0){
        alreadyusedprexixes <- prefixestouse[(prefixestouse %in% usedprefixes)]
        if (length(alreadyusedprexixes) > 0){
            flog.warn(paste("WARNING: THE FOLLOWING PREFIXES HAVE BEEN FOUND IN", as.character(opt$checkprefixexists)))
            flog.warn(paste0(alreadyusedprexixes, collapse = ", "))
            flog.warn("Check your filenames and substitution table. ABORTING NOW.")
            q()
        } else {
            flog.info("None of the proposed new prefixes have been used previously")
        }
    }
}
save.image(opt$projimage)

rawfastqsdf$NewFullFN <- file.path(opt$destination, rawfastqsdf$NewFN)

#Copy or move the files
if (opt$move == TRUE){
    manipcmd <- "mv"
} else {
    manipcmd <- "cp"
}

rencommands <- sapply(1:nrow(rawfastqsdf), function(x) { paste(manipcmd, rawfastqsdf$Filepaths[x], rawfastqsdf$NewFullFN[x], collapse = " ") } )

if (opt$simulate){
    flog.info("This is a simulation. Commands for renaming when not simulating are:")
    print(rencommands)
} else {
    for (cmd in rencommands){
        flog.info(paste("Renaming with command", cmd))
        system(cmd)
    }
}

save.image(opt$projimage)
