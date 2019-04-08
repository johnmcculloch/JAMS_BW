#' check_resources(opt=opt)
#'
#' JAMSalpha function
#' @export

check_resources<-function(opt=opt){
    opt$abort<-FALSE

    #Check for non-R software dependencies
    deps <- c("pigz", "wget", "sratoolkit", "samtools", "trimmomatic", "bowtie2", "megahit", "spades", "kraken2", "convert2bed", "bedtools")
    missingdep = FALSE
    for (dep in deps) {
        cmd = dep
        if (cmd=="sratoolkit") cmd = "fastq-dump"
        if (cmd=="spades") cmd = "spades.py"
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
    k2dbfiles<-list.files(opt$dbdir, recursive=TRUE, full.names=TRUE, include.dirs=FALSE, pattern="*.k2d")
    if(length(k2dbfiles) != 3){
        flog.info(paste0("A kraken2 database was not found in", opt$dbdir,". Check your database structure. Aborting."))
        opt$abort<-TRUE        
    }
    opt$krakendb<-paste0(unlist(strsplit(k2dbfiles[1], split="/"))[1:length(unlist(strsplit(k2dbfiles[1], split="/")))-1], collapse="/")
    if(!(dir.exists(opt$krakendb))){
        flog.info("JAMS kraken db directory not found. Aborting.")
        opt$abort<-TRUE
    } else {
        flog.info(paste("JAMS kraken2 database found at", opt$krakendb))
        krakendbsize<-file.size(file.path(opt$krakendb, "hash.k2d"))
        if(is.na(krakendbsize)){
            flog.info("Kraken2 database supplied seems to be corrupt or not valid. Check your database structure. Aborting.")
            opt$abort<-TRUE
        } else if (krakendbsize > opt$totmembytes) {
            flog.info("It seems there is not enough RAM memory to use the kraken2 database supplied. Try again on a system with more RAM. Aborting.")
            opt$abort<-TRUE
        } else {
            if(file.exists(file.path(opt$krakendb, "JAMSKdb.ver"))){
                opt$JAMS_Kdb_Version<-as.character(read.table(file.path(opt$krakendb, "JAMSKdb.ver"), header=FALSE, quote=NULL)[1,1])
            } else {
                opt$JAMS_Kdb_Version<-"UNDETERMINED"
            }
            flog.info(paste("The JAMS kraken2 database supplied is version", opt$JAMS_Kdb_Version, "and has", round((krakendbsize/1000000000), 1), "Gb."))
        }
    }

    #Check for host bowtie indices
    if(opt$host %in% c("mouse", "human")){
        flog.info("Checking for the presence of pre-compiled host genome bowtie2 indices.")
        #Set host prefix
        if(opt$host == "mouse"){
           opt$hostspecies="Mmusculus"
        } else if (opt$host == "human"){
           opt$hostspecies="Hsapiens"
        }
        indexfiles<-list.files(opt$dbdir, recursive=TRUE, full.names=TRUE, include.dirs=FALSE, pattern="*.bt2$")
        opt$relindexfiles<-indexfiles[grep(opt$hostspecies, indexfiles)]
        if(length(opt$relindexfiles) == 6){
            opt$indexpath<-paste0(unlist(strsplit(opt$relindexfiles[1], split="/"))[1:length(unlist(strsplit(opt$relindexfiles[1], split="/")))-1], collapse="/")
            flog.info(paste0("Pre-built bowtie2 indices for the ", opt$host, " genome were found at ", opt$indexpath))
        } else {
            flog.info(paste0("Pre-built bowtie2 indices for the ", opt$host, " genome were not found. Check your database structure. Aborting."))
            opt$abort<-TRUE
        }
    }

    return(opt)
}
