#' check_resources(opt=opt)
#'
#' JAMSalpha function
#' @export

check_resources<-function(opt=opt){
    opt$abort<-FALSE

    #Check for non-R software dependencies
    deps <- c("pigz", "wget", "sratoolkit", "samtools", "trimmomatic", "bowtie2", "megahit", "spades", "kraken2", "convert2bed", "bedtools")
    missingdep <- FALSE
    for (dep in deps) {
        cmd <- dep
        if (cmd == "sratoolkit") cmd = "fastq-dump"
        if (cmd == "spades") cmd = "spades.py"
        if (system(str_c("which ", cmd), ignore.stdout = TRUE)==1) { # not found
            flog.info(str_c("You are missing ", dep))
            missingdep = TRUE
        }
    }
    if (missingdep == TRUE) {
        flog.info("Please install missing dependencies before using JAMSalpha. Aborting.")
        opt$abort<-TRUE
    }

    #Check for kraken2 database
    flog.info("Checking for presence of kraken2 database")
    k2dbfiles <- list.files(opt$dbdir, recursive=TRUE, full.names=TRUE, include.dirs=FALSE, pattern="*.k2d")
    if (length(k2dbfiles) != 3){
        flog.info(paste0("A kraken2 database was not found in", opt$dbdir,". Check your database structure. Aborting."))
        opt$abort <- TRUE
    }
    opt$krakendb <- paste0(unlist(strsplit(k2dbfiles[1], split="/"))[1:length(unlist(strsplit(k2dbfiles[1], split="/")))-1], collapse="/")
    if (!(dir.exists(opt$krakendb))){
        flog.info("JAMS kraken db directory not found. Aborting.")
        opt$abort <- TRUE
    } else {
        flog.info(paste("JAMS kraken2 database found at", opt$krakendb))
        krakendbsize <- file.size(file.path(opt$krakendb, "hash.k2d"))
        if (is.na(krakendbsize)){
            flog.info("Kraken2 database supplied seems to be corrupt or not valid. Check your database structure. Aborting.")
            opt$abort <- TRUE
        } else if (krakendbsize > opt$totmembytes) {
            flog.info("It seems there is not enough RAM memory to use the kraken2 database supplied. Try again on a system with more RAM. Aborting.")
            opt$abort <- TRUE
        } else {
            if (file.exists(file.path(opt$krakendb, "JAMSKdb.ver"))){
                opt$JAMS_Kdb_Version <- as.character(read.table(file.path(opt$krakendb, "JAMSKdb.ver"), header=FALSE, quote=NULL)[1,1])
            } else {
                opt$JAMS_Kdb_Version <- "UNDETERMINED"
            }
            flog.info(paste("The JAMS kraken2 database supplied is version", opt$JAMS_Kdb_Version, "and has", round((krakendbsize/1000000000), 1), "Gb."))
        }
    }

    #Check for blast databases
    blastfiles <- list.files(opt$dbdir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="\\.[p|n]hr")
    if (length(blastfiles) > 0){
        opt$blastanalyses <- sapply(1:length(blastfiles), function (x) { unlist(strsplit(blastfiles[x], split="/"))[(length(unlist(strsplit(blastfiles[x], split="/"))) - 1)] })
        opt$blastpath <- paste0(unlist(strsplit(blastfiles[1], split="/"))[1:(length(unlist(strsplit(blastfiles[1], split="/")))-2)], collapse = "/")
        flog.info(paste("Found blast databases for the following blast analyses:", paste0(opt$blastanalyses, collapse = ", ")))
        flog.info(paste("Blast databases are in:", opt$blastpath))
    } else {
        opt$blastanalyses <- NULL
        flog.info("No BLAST databases found in the JAMSdb folder supplied.")
    }

    #Set host species to none if isolates or if using contigs.
    if (opt$analysis=="isolate"){
        opt$host <- "none"
    }

    if (!(is.null(opt$contigs))){
        opt$host <- "none"
    }

    #Check for host bowtie indices
    if (opt$host %in% c("mouse", "human")){
        flog.info("Checking for the presence of pre-compiled host genome bowtie2 indices.")
        #Set host prefix
        if (opt$host == "mouse"){
           opt$hostspecies <- "Mmusculus"
        } else if (opt$host == "human"){
           opt$hostspecies <- "Hsapiens"
        }
        indexfiles <- list.files(opt$dbdir, recursive=TRUE, full.names=TRUE, include.dirs=FALSE, pattern="*.bt2$")
        opt$relindexfiles <- indexfiles[grep(opt$hostspecies, indexfiles)]
        if (length(opt$relindexfiles) == 6){
            opt$indexpath <- paste0(unlist(strsplit(opt$relindexfiles[1], split="/"))[1:length(unlist(strsplit(opt$relindexfiles[1], split="/")))-1], collapse="/")
            flog.info(paste0("Pre-built bowtie2 indices for the ", opt$host, " genome were found at ", opt$indexpath))
        } else {
            flog.info(paste0("Pre-built bowtie2 indices for the ", opt$host, " genome were not found. Check your database structure. Aborting."))
            opt$abort <- TRUE
        }
    } else if (opt$host == "none") {
        flog.info("No need to eliminate sequencing reads coming from a host species.")
    } else {
        flog.info("Host species supplied is not human, mouse or none.")
        if (file.exists(opt$host) == TRUE){
            flog.info("Host species bowtie index has been supplied as a path.")
            opt$relindexfiles <- list.files(opt$host, recursive=TRUE, full.names=TRUE, include.dirs=FALSE, pattern="*.bt2$")
            #Need 6 index files
            if (length(opt$relindexfiles) == 6){
                opt$indexpath <- paste0(unlist(strsplit(opt$relindexfiles[1], split="/"))[1:length(unlist(strsplit(opt$relindexfiles[1], split="/")))-1], collapse="/")
                indexfile1 <- unlist(strsplit(opt$relindexfiles[1], split="/"))[length(unlist(strsplit(opt$relindexfiles[1], split="/")))]
                opt$host <- unlist(strsplit(indexfile1, split="\\."))[1]
                opt$hostspecies <- opt$host
                flog.info(paste0("Pre-built bowtie2 indices for the ", opt$host, " genome were found at ", opt$indexpath))
            } else {
                flog.info(paste0("Pre-built bowtie2 indices for the host genome were not found or look weird. Check the path supplied. Aborting."))
                opt$abort <- TRUE
            }
        } else {
            flog.info("Host species genome has been supplied as an NCBI species taxid or organism name. Checking Tax ID and best genome to download.")
            flog.info("Trying as Tax ID first.")
            host_assemblies <- get_genomes_NCBI(organisms = "vertebrate_mammalian", outputdir = NULL, species_taxid = opt$host, nobs = TRUE, fileformat = "fasta", simulate = TRUE)
            host_assemblies <- subset(host_assemblies, genome_rep == "Full")
            if(nrow(host_assemblies) == 0){
                flog.info("Not Tax ID. Trying as organism_name.")
                host_assemblies <- get_genomes_NCBI(organisms = "vertebrate_mammalian", outputdir = NULL, organism_name = opt$host, nobs = TRUE, fileformat = "fasta", simulate = TRUE)
            }
            #Only consider full genomes
            host_assemblies <- subset(host_assemblies, genome_rep == "Full")
            if(nrow(host_assemblies) > 0){
                flog.info("Deteriming best host genome assembly to download.")
                #Give priority to reference and then representative genomes, if available.
                if(any(c("reference genome", "representative genome") %in% host_assemblies$refseq_category)){
                    refseqcats <- (c("reference genome", "representative genome") %in% host_assemblies$refseq_category)
                    host_assemblies <- subset(host_assemblies, refseq_category == c("reference genome", "representative genome")[min(which(refseqcats==TRUE))])
                }
                #Go with the most recent
                host_assemblies <- host_assemblies[order(host_assemblies$seq_rel_date, decreasing=TRUE), ][1 , ]
                opt$host_accession_number <- host_assemblies$assembly_accession
                flog.info(paste("Will build index with", host_assemblies$organism_name, "genome under accession number", opt$host_accession_number, "of", host_assemblies$seq_rel_date))
                opt$host <- paste(host_assemblies$organism_name, opt$host_accession_number, sep="_")
                opt$hostspecies <- host_assemblies$organism_name
            } else {
                flog.info("Could not find the Tax ID or organism name to download. Aborting now.")
                opt$abort <- TRUE
            }
        }
    }

    #If testing dependencies and leave is true then it is time to abort.
    if(opt$testdependencies == TRUE){
        opt$abort <- TRUE
    }

    return(opt)
}
