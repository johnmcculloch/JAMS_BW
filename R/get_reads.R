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
    if(opt$workdir != opt$sampledir){
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
        #If in Biowulf, module load sra-toolkit
        slurmjobid <- as.character(Sys.getenv("SLURM_JOB_ID"))
        #if (nchar(slurmjobid) > 3){
        #    flog.info("You are probably on Biowulf. Will use NIH HPC version of SRA toolkit to avoid caching issues.")
        #    commandtorun <- paste(file.path(opt$bindir, "getreadsSRAonBW.sh"), opt$sraaccession, sep =  " ")
        #} else {
            #commandtorun <- paste("fastq-dump --skip-technical --readids --read-filter pass --dumpbase --split-e --clip", opt$sraaccession, sep = " ")
            commandtorun <- paste("fasterq-dump --split-files --skip-technical --threads", opt$threads, opt$sraaccession, sep = " ")
        #}
        system(commandtorun)

        #Eliminate unpaired read because split-e was used. This means that unpaired reads from a paired Run will be dumped into a third file.
        if (length(list.files(pattern = "fastq") > 2)){
            file.remove(paste0(opt$sraaccession, "_pass.fastq"))
        }
    } else {
        opt$readorigin <- "supplied"
        #Copy supplied reads into readsdir
        if (!(is.null(opt$readstarball))){
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

        } else {
            #Copy R1, R2 and U reads
            fastqsinopt <- names(opt)[which(names(opt) %in% c("R1", "R2", "SE"))]
            if(length(fastqsinopt) < 1){
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

                targetfilename <- file.path(readsworkdir, paste(paste(opt$prefix, rtype, sep="_"), suffix, sep="."))
                flog.info(paste("Copying", rtype, myfastqs[[f]]))
                system2('cp', args = c(myfastqs[[f]], targetfilename), stdout = TRUE, stderr = TRUE)
                #In case files are gzipped, ungzip them
                if (suffix == "fastq.gz"){
                    flog.info(paste("Decompressing gzipped file", rtype))
                    commandtorun <- paste("pigz -d", targetfilename, collapse = " ")
                    system(commandtorun)
                } else if (suffix == "fastq.bz2"){
                    flog.info(paste("Decompressing bzipped file", rtype))
                    commandtorun <- paste("bzip2 -d", targetfilename, collapse = " ")
                    system(commandtorun)
                }
            }
        }
    }

    #find out number of files and rename accordingly
    rawfiles <- sort(list.files())
    if(length(rawfiles) == 1){
        flog.info("Found a single fastq file. Will assume this is of unpaired reads.")
        file.rename(rawfiles, paste0(opt$prefix, "_SE.fastq"))
        opt$libstructure <- "singleend"
    } else if (length(rawfiles) == 2){
        flog.info("Found two fastq files. Will assume these are paired reads.")
        file.rename(rawfiles, paste0(opt$prefix, c("_R1.fastq", "_R2.fastq")))
        opt$libstructure <- "pairedend"
    } else if (length(rawfiles) == 0){
        flog.info("Found NO files. Aborting now.")
        q()
    } else {
        flog.info(paste("Found", length(rawfiles), "files. Currently, JAMS supports only either paired reads OR single-end reads, but not both kinds simultaneously. Aborting now."))
        q()
    }

    rawfastqs <- sort(list.files(pattern="*.fastq$"))
    opt$rawreads <- rawfastqs

    return(opt)
}
