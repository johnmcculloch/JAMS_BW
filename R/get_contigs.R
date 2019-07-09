#' get_contigs
#'
#' JAMSalpha function
#' @export

get_contigs <- function(opt = NULL){
    setwd(opt$workdir)

    #Set working kraken db path for increased speed.
    #if(opt$workdir != opt$sampledir){
    #    dir.create("krakendb", showWarnings = FALSE, recursive = FALSE)
    #    flog.info("Copying krakendb to temporary file for speedier taxonomic classification.")
    #    #copy the whole database to tempfile for speed, if applicable.
    #    system2('cp', args=c("-R", file.path(opt$krakendb,"*"), "krakendb/"))
    #    opt$workingkrakendb<-file.path(opt$workdir, "krakendb")
    #} else {
        opt$workingkrakendb <- opt$krakendb
    #}

    #find out if contigs are ready or need to be assembled from reads

    if (!(is.null(opt$contigs))){
        #copy contigs
        flog.info("Input sequence supplied are contigs.")
        opt <- prepare_contigs_for_JAMS(opt = opt, fastafile = opt$contigs)

        #Write NHcontigs to the system
        NHcontigsfastaout <- paste(opt$prefix, "contigs.fasta", sep="_")
        write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$workdir, NHcontigsfastaout))

        #Bank assembly files to project directory
        if(opt$workdir != opt$sampledir){
            #Write contigs from opt directly to system
            write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$sampledir, NHcontigsfastaout))
        }

    } else if (!(is.null(opt$assemblyaccession))){
        flog.info(paste("Will download contigs from NCBI GenBank under accession number", opt$assemblyaccession[1]))
        assemblytargetname <- opt$assemblyaccession[1]
        contigsaccession <- get_genomes_NCBI(organisms = "bacteria", outputdir = opt$workdir, assembly_accession = assemblytargetname, ntop = 1, fileformat = "fasta", simulate = FALSE)
        if (nrow(contigsaccession) < 1){
            flog.info(paste("A valid entry with accession number", assemblytargetname, "was not found in GenBank. Only bacterial accession numbers are accepted for the time being. Aborting now."))
            opt$abort <- TRUE
        }
        gztn <- paste(assemblytargetname, "fasta", "gz", sep = ".")
        gbfn <- paste(assemblytargetname, "fasta", sep = ".")
        file.rename(contigsaccession$destfn, gztn)
        system(paste("gunzip", gztn))

        #Set opt$contigs to NULL in case the user has tried to set that as an input too
        opt <- prepare_contigs_for_JAMS(opt = opt, fastafile = gbfn)

        #Write NHcontigs to the system
        NHcontigsfastaout <- paste(opt$prefix, "contigs.fasta", sep="_")
        write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$workdir, NHcontigsfastaout))

        #Bank assembly files to project directory
        if(opt$workdir != opt$sampledir){
            #Write contigs from opt directly to system
            write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$sampledir, NHcontigsfastaout))
        }

    } else {
        flog.info("Input sequences are reads.")
        flog.info("Will assemble reads into contigs.")
        opt <- assemble_contigs(opt=opt)
    }

    flog.info("Contigs have been obtained.")

    return(opt)
}
