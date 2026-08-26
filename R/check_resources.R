#' check_resources(opt = opt, application = "WGS")
#'
#' JAMSalpha function
#' @export

check_resources <- function(opt = opt, applications_to_check = c("base", "reads", "assembly", "kraken", "host", "annotation", "16S", "blast")){
    opt$abort <- FALSE

    deplist <- list(base = c("pigz", "wget"), reads = c("sratoolkit", "trimmomatic", "seqtk"), host = c("bowtie2"), assembly = c("megahit", "spades", "samtools", "metabat2", "checkm2"), annotation = c("prokka", "convert2bed", "bedtools"), kraken = "kraken2")

    dependencies_to_check <- unname(unlist(deplist[applications_to_check]))

    #Map from dependency name to the actual executable that is called on the command line
    exec_for_dep <- function(dep){
        switch(dep,
            "sratoolkit" = "fasterq-dump",
            "spades"     = "spades.py",
            dep)
    }

    #' get_dep_version
    #' Try to extract a version string for a given executable.
    #' Different tools expose their version through different flags and print it
    #' in different formats, so we try a sequence of common flags and then apply
    #' a light, tool-aware parse of whatever gets printed.
    get_dep_version <- function(dep, cmd){
        #Tool-specific version flags. Order matters: first flag that yields output wins.
        flag_lookup <- list(
            "pigz"         = c("--version"),
            "wget"         = c("--version"),
            "fasterq-dump" = c("--version"),
            "trimmomatic"  = c("-version", "--version"),
            "seqtk"        = c("--version"), # seqtk prints version to stderr with no valid flag; handled below
            "bowtie2"      = c("--version"),
            "megahit"      = c("--version", "-v"),
            "spades.py"    = c("--version", "-v"),
            "samtools"     = c("--version"),
            "metabat2"     = c("--help"), # prints "MetaBAT 2 (2.x ...)" banner; parsed below
            "checkm2"      = c("--version"),
            "prokka"       = c("--version"),
            "convert2bed"  = c("--version"),
            "bedtools"     = c("--version"),
            "kraken2"      = c("--version")
        )

        raw <- character(0)

        #seqtk does not accept any version flag; passing one prints an error to
        #stderr ("unrecognized command '--version'. Abort!") which would otherwise
        #be mistaken for version output. Invoke it bare instead: run with no args
        #it prints a proper "Version: 1.x" banner to stderr.
        if (cmd == "seqtk"){
            raw <- tryCatch(
                suppressWarnings(system("seqtk 2>&1", intern = TRUE)),
                error = function(e) character(0)
            )
            raw <- raw[nzchar(raw)]
        } else {
            flags <- flag_lookup[[cmd]]
            if (is.null(flags)){
                flags <- c("--version", "-version", "-v", "--help")
            }
            for (flg in flags){
                #Capture both stdout and stderr, since many bioinformatics tools
                #print their version banner to stderr.
                out <- tryCatch(
                    suppressWarnings(system(str_c(cmd, " ", flg, " 2>&1"), intern = TRUE, ignore.stderr = FALSE)),
                    error = function(e) character(0)
                )
                out <- out[nzchar(out)]
                if (length(out) > 0){
                    raw <- out
                    break
                }
            }
        }

        if (length(raw) == 0){
            return("UNDETERMINED")
        }

        #Collapse to a single searchable blob but keep the first few lines for parsing.
        blob <- paste(head(raw, 10), collapse = " ")

        #Prefer a "Version: x.y.z" style token if present (samtools, seqtk, bowtie2 etc.)
        ver <- NA_character_
        m <- regmatches(blob, regexpr("[Vv]ersion[:]?[[:space:]]*[Vv]?[0-9]+\\.[0-9][0-9A-Za-z\\._\\-]*", blob))
        if (length(m) > 0 && nzchar(m)){
            #Strip the leading "Version" label, keep the numeric part.
            ver <- sub("^[Vv]ersion[:]?[[:space:]]*", "", m)
        } else {
            #Otherwise grab the first thing that looks like a dotted version number
            #(requires a dot, so a stray standalone integer such as the "2" in
            #"MetaBAT 2" is not mistaken for the version).
            m <- regmatches(blob, regexpr("[0-9]+\\.[0-9][0-9A-Za-z\\._\\-]*", blob))
            if (length(m) > 0 && nzchar(m)){
                ver <- m
            }
        }

        if (is.na(ver) || !nzchar(ver)){
            #As a last resort keep the first non-empty line, trimmed, so we at
            #least record something informative.
            ver <- trimws(raw[1])
        }

        #Tidy trailing punctuation
        ver <- sub("^[vV]", "", ver)
        ver <- gsub("[[:space:],;]+$", "", ver)

        return(ver)
    }

    #Check for non-R software dependencies and record versions
    missingdep <- FALSE
    version_records <- list()
    for (dep in dependencies_to_check) {
        cmd <- exec_for_dep(dep)
        if (system(str_c("which ", cmd), ignore.stdout = TRUE) == 1) { # not found
            flog.info(str_c("You are missing ", dep))
            missingdep <- TRUE
            version_records[[dep]] <- data.frame(Dependency = dep, Executable = cmd, Found = FALSE, Version = NA_character_, stringsAsFactors = FALSE)
        } else {
            depver <- get_dep_version(dep, cmd)
            flog.info(str_c("Found ", dep, " (", cmd, ") version ", depver))
            version_records[[dep]] <- data.frame(Dependency = dep, Executable = cmd, Found = TRUE, Version = depver, stringsAsFactors = FALSE)
        }
    }

    #Bequeath a tidy data frame of dependency versions to opt
    if (length(version_records) > 0){
        opt$dependency_versions <- do.call(rbind, version_records)
        rownames(opt$dependency_versions) <- NULL
    } else {
        opt$dependency_versions <- data.frame(Dependency = character(0), Executable = character(0), Found = logical(0), Version = character(0), stringsAsFactors = FALSE)
    }

    if (missingdep == TRUE) {
        flog.info("Please install missing dependencies before using JAMSalpha. Aborting.")
        opt$abort <- TRUE
    }

    if ("kraken" %in% applications_to_check){
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
                #Pass path to opt$opt$workingkrakendb for backwards compatibility
                opt$workingkrakendb <- opt$krakendb
                if (file.exists(file.path(opt$krakendb, "JAMSKdb.ver"))){
                    opt$JAMS_Kdb_Version <- as.character(read.table(file.path(opt$krakendb, "JAMSKdb.ver"), header=FALSE, quote=NULL)[1,1])
                } else {
                    opt$JAMS_Kdb_Version <- "UNDETERMINED"
                }
                flog.info(paste("The JAMS kraken2 database supplied is version", opt$JAMS_Kdb_Version, "and has", round((krakendbsize/1000000000), 1), "Gb."))
            }
        }
        flog.info("Checking for presence of CheckM database")
        checkmfiles <- list.files(opt$dbdir, recursive=TRUE, full.names=TRUE, include.dirs=FALSE, pattern="*.dmnd")
        if (length(checkmfiles) > 0){
            #opt$CheckMdb <- paste0(unlist(strsplit(checkmfiles[1], split="/"))[1:length(unlist(strsplit(checkmfiles[1], split="/")))-1], collapse="/")
            opt$CheckMdb <- checkmfiles[1]
            flog.info(paste("Found CheckM2 database at", opt$CheckMdb))
        } else {
            flog.warn("CheckM database not found. Will not evaluate genome completeness using CheckM")
            #Ensure NULL if db is absent
            opt$CheckMdb <- NULL
        }
    }

    if ("blast" %in% applications_to_check){
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
    }

    if ("host" %in% applications_to_check){
        #Set host species to none if isolates or if using contigs.
        if (opt$analysis == "isolate"){
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
                host_assemblies <- get_genomes_NCBI(organisms = NULL, outputdir = NULL, species_taxid = opt$host, nobs = TRUE, fileformat = "fasta", simulate = TRUE)
                host_assemblies <- subset(host_assemblies, genome_rep == "Full")
                if (nrow(host_assemblies) == 0){
                    flog.info("Not Tax ID. Trying as organism_name.")
                    host_assemblies <- get_genomes_NCBI(organisms = NULL, outputdir = NULL, organism_name = opt$host, nobs = TRUE, fileformat = "fasta", simulate = TRUE)
                }
                #Only consider full genomes
                host_assemblies <- subset(host_assemblies, genome_rep == "Full")
                if (nrow(host_assemblies) > 0){
                    flog.info("Deteriming best host genome assembly to download.")
                    #Give priority to reference and then representative genomes, if available.
                    if (any(c("reference genome", "representative genome") %in% host_assemblies$refseq_category)){
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
    }

    #If testing dependencies and leave is true then it is time to abort.
    if (opt$testdependencies == TRUE){
        opt$abort <- TRUE
    } else {
        #Detect, create and set temporary folder
        #Test slurm/Biowulf specific situation
        if (is.null(opt$tmpdir) & (file.exists(file.path("/lscratch", as.character(Sys.getenv("SLURM_JOB_ID")))))){
            opt$tempdir <- file.path("/lscratch", as.character(Sys.getenv("SLURM_JOB_ID")))
        }

        if (!is.null(opt$tempdir)){
            tmpdirok <- FALSE
            if (file.exists(opt$tempdir)){
                tmpdirok <- dir.create(file.path(opt$tempdir, opt$prefix), showWarnings = TRUE, recursive = TRUE)
                if (tmpdirok){
                    opt$tempdir <- file.path(opt$tempdir, opt$prefix)
                    flog.info(paste("Temporary folder set to", opt$tempdir))
                } else {
                    flog.warn(paste("Unable to set temporary folder to", file.path(opt$tempdir, opt$prefix)))
                    opt$tempdir <- NULL
                }
            } else {
                flog.info(paste("Temporary folder", opt$tempdir, "does not exist and will not be used."))
                opt$tempdir <- NULL
            }
        }
    }

    return(opt)
}