#' get_reads
#'
#' JAMSalpha function
#' @export

get_reads <- function(opt = NULL){

    setwd(opt$workdir)

    #Make final readsdir
    opt$readsdir <- file.path(opt$sampledir, paste0(opt$prefix, "_reads"))
    flog.info("Creating directory to hold reads.")
    dir.create(opt$readsdir, showWarnings = FALSE, recursive = TRUE)

    #Make readsdir in tempfile if applicable
    if (opt$workdir != opt$sampledir){
        readsworkdir <- file.path(opt$workdir, "reads")
        flog.info("Creating temporary directory to hold reads.")
        dir.create(readsworkdir, showWarnings = FALSE, recursive = TRUE)
    } else {
        readsworkdir <- opt$readsdir
    }

    setwd(readsworkdir)
    #If there is an accession number supplied, get reads from there
    if (!(is.null(opt$sraaccession))){
        #Download reads from accession
        opt$readorigin <- "sraaccession"
        flog.info(paste("Downloading from NCBI SRA reads pertaining to run", opt$sraaccession))
        #if tempdir exists use that for fasterq-dump temporary folder, as fasterq uses up to 7 x the final fastq size in intermediate files.
        if (!is.null(opt$tempdir)){
            #Prefetch to tempdir, much, much faster
            #set temporary folder to contain both the SRA accession AND the prefix. This is safer in case the same temporary folder is used for multiple runs concomitantly. Even if two runs in parallel for the same SRA accession are run, there would be no clash, as long as prefixes are different. Who would launch two exact same JAMSalpha commands in parallel? Right?
            SRAtempdir <- file.path(opt$tempdir, opt$sraaccession)
            fetchcmd <- paste("prefetch", opt$sraaccession, "-O", SRAtempdir, "--max-size 1t", sep = " ")
            system(fetchcmd)
            #go to prefetch folder and do fasterq-dump
            setwd(SRAtempdir)
            commandtorun <- paste("fasterq-dump --split-files --skip-technical --threads", opt$threads, "--outdir", opt$readsdir, opt$sraaccession, sep = " ")
            system(commandtorun)
            #Return to where we were
            setwd(readsworkdir)
            #Purge temporary SRA folder
            unlink(SRAtempdir, recursive = TRUE)
        } else {
            commandtorun <- paste("fasterq-dump --split-files --skip-technical --threads", opt$threads, opt$sraaccession, sep = " ")
            system(commandtorun)
        }

        #Eliminate unpaired read because split-e was used. This means that unpaired reads from a paired Run will be dumped into a third file.
        if (length(list.files(pattern = "fastq") > 2)){
            file.remove(paste0(opt$sraaccession, "_pass.fastq"))
        }
        #gzip reads to save intermediate file size and conserve disk space
        flog.info("Compressing SRA reads for maximizing speed and minimizing storage space.")
        for (RN in c(1, 2)){
            if (file.exists(paste(opt$prefix, paste(RN, "fastq", sep = "."), sep = "_"))){
                gzcmd <- paste("pigz -p", opt$threads, paste(opt$prefix, paste(RN, "fastq", sep = "."), sep = "_"), sep = " ")
                system(gzcmd)
                file.rename(from = paste(opt$prefix, paste(RN, "fastq.gz", sep = "."), sep = "_"), to = paste(opt$prefix, paste(paste0("R", RN), "fastq.gz", sep = "."), sep = "_"))
            }
        }

    } else {
        opt$readorigin <- "supplied"
        #Copy supplied reads into readsdir
        if (!(is.null(opt$readstarball))){
            #Frowned upon, but will still maintain this option for the time being.
            #Copy tarball
            if ((summary(file(opt$readstarball))$class) != "gzfile" ){
                flog.info("You chose a tarball as input but it does not look like a tar.gz file. Aborting now.")
                q()
            }
            targetfilename <- file.path(readsworkdir, paste(opt$prefix, "rawreads.tar.gz", sep = "_"))
            system2('cp', args = c(opt$readstarball, targetfilename), stdout = TRUE, stderr = TRUE)
            #Decompress using pigz
            commandtorun <- paste("pigz -p", opt$threads, "-dc", targetfilename, "| tar xf -", collapse = " ")
            flog.info("Decompressing tarball with reads.")
            system(commandtorun)
            file.remove(targetfilename)
            #Now, compress the individual files
            list.files(pattern = "fastq")
            for (fn in list.files(pattern = "fastq")){
                system(paste("pigz -p", opt$threads, fn))
            }

        } else {
            #Copy R1, R2 and U reads
            fastqsinopt <- names(opt)[which(names(opt) %in% c("R1", "R2", "SE"))]
            if (length(fastqsinopt) < 1){
                flog.info("You must supply either reads to assemble or contigs or a GenBank SRA accession number as input. None of these were found. Check arguments passed JAMSalpha. Aborting now.")
                q()
            }
            flog.info(paste("You supplied", paste(fastqsinopt, collapse = ", "), "reads."))
            myfastqs <- opt[fastqsinopt]

            #Check if files actually exist
            if (!(all(file.exists(as.character(myfastqs))))){
                flog.info("At least one of the input files is missing. Aborting now.")
                q()
            }

            #Copy reads
            for (f in 1:length(myfastqs)){
                rtype <- names(myfastqs)[f]
                if (filetype(myfastqs[[f]]) == "gzfile"){
                    suffix <- "fastq.gz"
                } else if (filetype(myfastqs[[f]]) == "bzfile"){
                    suffix <- "fastq.bz2"
                } else {
                    suffix <- "fastq"
                }

                #Set target filename 
                targetfilename <- file.path(readsworkdir, paste(paste(opt$prefix, rtype, sep = "_"), "fastq.gz", sep="."))
                #In case files are ungzipped, gzip them to save intermediate file space
                if (suffix == "fastq"){
                    flog.info(paste("Copying and compressing", rtype, myfastqs[[f]]))
                    commandtorun <- paste("pigz --processes", opt$threads, "-c", myfastqs[[f]], ">", targetfilename, collapse = " ")
                    system(commandtorun)
                } else if (suffix == "fastq.bz2"){
                    #What is wrong with you, you maniac?
                    flog.warn("Currently, JAMS only takes fastq or fastq.gz files as read inputs. Please re-select or convert your input read sequence files. Aborting now.")
                    q()
                } else if (suffix == "fastq.gz") {
                    #copy and rename file. Best to copy rather than read from origin, as disk access can be slow depending on original storage location
                    system2('cp', args = c(myfastqs[[f]], targetfilename), stdout = TRUE, stderr = TRUE)
                }
            }
        }
    }

    #find out number of files and rename accordingly
    rawfiles <- sort(list.files(pattern = "fastq.gz$"))
    if (length(rawfiles) == 1){
        flog.info("Found a single fastq file. Will assume this is of unpaired reads.")
        file.rename(rawfiles, paste0(opt$prefix, "_SE.fastq.gz"))
        opt$libstructure <- "singleend"
    } else if (length(rawfiles) == 2){
        flog.info("Found two fastq files. Will assume these are paired reads.")
        file.rename(rawfiles, paste0(opt$prefix, c("_R1.fastq.gz", "_R2.fastq.gz")))
        opt$libstructure <- "pairedend"
    } else if (length(rawfiles) == 0){
        flog.warn("Found NO files. Aborting now.")
        q()
    } else {
        flog.info(paste("Found", length(rawfiles), "files. Currently, JAMS supports only either paired reads OR single-end reads, but not both kinds simultaneously. Aborting now."))
        q()
    }

    opt$rawreads <- sort(list.files(pattern = "*.fastq.gz$"))
    opt$fastqstats <- countfastq_files(fastqfiles = opt$rawreads, threads = opt$threads, type_tag = "Raw")

    setwd(opt$workdir)

    return(opt)
}
