#' download_SRA_reads(SRAaccessions = NULL, outfolder = NULL, tempfolder = NULL, prefetch = TRUE, threads = 8)
#'
#' Wrapper for safely downloading SRA reads within the JAMS framework. You will need to have SRAtoolkit installed on your system.
#' @export

download_SRA_reads <- function(SRAaccessions = NULL, outfolder = NULL, tempfolder = NULL, prefetch = TRUE, threads = 8){
    flog.info("Retrieving fastq files from SRA.")
    flog.info(paste("There are", length(SRAaccessions), "SRA accessions to download. Will download and process these serially to minimize storage space use."))
    curr_wd <- getwd()
    #Check or make final destination folder
    if (is.null(outfolder)){
        flog.warn("Output folder not set in arguments. Using current folder as output folder.")
        outfolder <- getwd()
    }

    if (!file.exists(outfolder)){
        flog.info(paste("Creating final output folder to hold processed reads at", outfolder))
        dir.create(outfolder, recursive = TRUE)
    }

    #Check temporary folder
    if (!is.null(tempfolder)){
        if (file.exists(tempfolder)){
            flog.info(paste("Temporary download and processing folder is set to", tempfolder))
        } else {
            flog.info(paste("Unable to find folder", tempfolder, "so temporary folder will NOT be used."))
            tempfolder <- NULL
        }
    }

    #Process each accession serially
    for (SRAacc in SRAaccessions){
        flog.info(paste("Processing", SRAacc))

        if (!is.null(tempfolder)){
            workfolder <- tempfolder
        } else {
            workfolder <- outfolder
        }

        setwd(workfolder)

        if (prefetch == TRUE) {
            fetchcmd <- paste("prefetch", SRAacc, "-O", workfolder, "--max-size 1t", sep = " ")
            system(fetchcmd)
            #go to prefetch folder and do fasterq-dump
            setwd(file.path(workfolder, SRAacc))
            commandtorun <- paste("fasterq-dump --split-files --skip-technical --threads", threads, "--outdir", workfolder, SRAacc, sep = " ")
            system(commandtorun)
        } else {
            commandtorun <- paste("fasterq-dump --split-files --skip-technical --threads", threads, "--outdir", workfolder, SRAacc, sep = " ")
            system(commandtorun)
        }

        #Process reads, eliminate leftovers and gz compress.
        setwd(workfolder)
        for (RN in c(1, 2)){
            curr_SRA_fn <- paste(paste(SRAacc, RN, sep = "_"), "fastq", sep = ".")
            #Check if file exists
            if (file.exists(curr_SRA_fn)){
                file.rename(from = curr_SRA_fn, to = paste(paste(SRAacc, paste0("R", RN), sep = "_"), "fastq", sep = "."))
                #Compress
                system2('pigz', args = paste("-p", threads, paste(paste(SRAacc, paste0("R", RN), sep = "_"), "fastq", sep = ".")))
                #Move files to final destination if applicable
                if (workfolder != outfolder){
                    file.copy(from = paste(paste(SRAacc, paste0("R", RN), sep = "_"), "fastq", "gz", sep = "."), to = file.path(outfolder, paste(paste(SRAacc, paste0("R", RN), sep = "_"), "fastq", "gz", sep = ".")))
                    file.remove(paste(paste(SRAacc, paste0("R", RN), sep = "_"), "fastq", "gz", sep = "."))
                }
            } else {
                flog.warn(paste("Could not find", curr_SRA_fn))
            }
        }

        #Cleanup prefetch or intermediate files
        if (file.exists(paste0(SRAacc, "_pass.fastq"))){
            file.remove(paste0(SRAacc, "_pass.fastq"))
        }

        if (file.exists(file.path(workfolder, SRAacc))){
            unlink(file.path(workfolder, SRAacc), recursive = TRUE)
        }

    }

    #Return to original folder
    setwd(curr_wd)

}
