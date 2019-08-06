#' load_jamsfiles_from_system(path = ".", onlysamples=NULL, loadfromscratch = TRUE, list.data = NULL))
#' JAMSbeta function
#'
#' Loads all JAMS files from system
#' @export

load_jamsfiles_from_system <- function(path = ".", onlysamples = NULL, loadfromscratch = TRUE, list.data = NULL){
    fljams = list.files(path = path, pattern="*.jams$")

    if (length(fljams)<1){
        stop("There are no .jams files in the specified path.")
    }

    currpath <- getwd()
    #Create temp files directory if needed.
    jamstempfilespath <- paste(currpath, "jamstempfiles", sep = "/")
    if(!(dir.exists(jamstempfilespath))){
        print("Creating directory to hold raw data.")
        dir.create(jamstempfilespath, showWarnings = FALSE, recursive = FALSE)
    }

    #If not loading from scratch, check for existing data to append to.
    if(loadfromscratch != TRUE){
        if(missing(list.data)){
            stop("If you want to append data to an existing project, you must supply a list.data object to append to.  If you want to load everything from scratch, then set loadfromscratch = TRUE.")
        }
    }

    #Restrict loading to only certain samples if required.
    if(!missing(onlysamples)){
        print(paste("Will only load .jams files for samples", paste0(onlysamples, collapse=", ")))
        samplesIhavejams<-unique(gsub(".jams$", "", fljams))
        if(length(onlysamples %in% samplesIhavejams) < length (onlysamples)){
            missingjamsfiles<-paste0(onlysamples[which(!(onlysamples %in% samplesIhavejams))], collapse=", ")
            stop(paste("Could not find .jams files for samples", missingjamsfiles))
        }
        fljams<-paste(onlysamples, "jams", sep=".")
    }

    #If there are .jams files, expand them.
    if(length(fljams)>0){
        print("Please be patient. Depending on how much data you have this might take a while.")
        print(paste("There are", length(fljams), "objects to expand."))
        fljamswp<-paste(path, fljams, sep="/")
        for(f in 1:length(fljams)){
            print(paste("Expanding ", fljams[f], ", file ", f, "/", length(fljams), "...", sep = ""))
            untar(fljamswp[f], list = FALSE, exdir = jamstempfilespath, compressed = TRUE, verbose = TRUE)
        }
    }

    fl = list.files(path = jamstempfilespath, pattern="*.tsv$")
    on = gsub("*.tsv$", "", fl)
    flwp<-paste(jamstempfilespath, fl, sep="/")

    if(missing(list.data)){
        list.data<-list()
    }

    print("Please be patient. Depending on how much data you have this might take a while.")

    lastpos<-length(names(list.data))

    print(paste("There are", length(fl), "objects to load."))
    for (i in 1:length(fl)){
        print(paste("Loading ", fl[i], ", file ", i, "/", length(fl), "...", sep = ""))
        list.data[[i+lastpos]]<-read.table(file=flwp[i], sep="\t", header=TRUE, quote="", skipNul=FALSE, fill=TRUE)
        names(list.data)[i+lastpos]<-on[i]
    }

    unlink(jamstempfilespath, recursive = TRUE)

    return(list.data)
}
