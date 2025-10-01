#' get_contigs
#'
#' JAMSalpha function
#' @export

get_contigs <- function(opt = NULL){
    setwd(opt$workdir)

    #Set working kraken db path for increased speed.
    if (opt$workdir != opt$sampledir){
        dir.create("krakendb", showWarnings = FALSE, recursive = FALSE)
        flog.info("Copying krakendb to temporary file for speedier taxonomic classification.")
        #copy the whole database to tempfile for speed, if applicable.
        system2('cp', args=c("-R", file.path(opt$krakendb,"*"), "krakendb/"))
        opt$workingkrakendb <- file.path(opt$workdir, "krakendb")
    } else {
        opt$workingkrakendb <- opt$krakendb
    }

    #find out if contigs are ready or need to be assembled from reads
    if (!(is.null(opt$contigsfasta))){
        #copy contigs
        flog.info("Input sequence supplied are contigs.")
        opt <- prepare_contigs_for_JAMS(opt = opt, fastafile = opt$contigsfasta)

        #Write NHcontigs to the system
        NHcontigsfastaout <- paste(opt$prefix, "contigs.fasta", sep="_")
        write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$workdir, NHcontigsfastaout))

        #Bank assembly files to project directory
        if(opt$workdir != opt$sampledir){
            #Write contigs from opt directly to system
            write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$sampledir, NHcontigsfastaout))
        }

    } else if (!(is.null(opt$assemblyaccession))){

        GCA_accessions <- unlist(strsplit(opt$assemblyaccession, split = ","))
        flog.info(paste("Will download contigs from NCBI GenBank under accession number(s)", paste0(GCA_accessions, collapse = ", ")))
        GCA_accession_list <- lapply(GCA_accessions, function(x) { retrieve_GCA_accession_URL(assemblytargetname = x) } )
        names(GCA_accession_list) <- GCA_accessions
        GCA_accession_list <- GCA_accession_list[which(sapply(1:length(GCA_accession_list), function (x) { !is.null(GCA_accession_list[[x]]) }))]

        if (length(GCA_accession_list) == 0){

            #Everything went t!t$ up.
            flog.warn("FATAL ERROR. None of the NCBI GenBank accessions could be retrieved. Check accession numbers and that you have esearch installed on your system. Aborting now.")
            IO_jams_workspace_image(opt = opt, operation = "save", verbose = FALSE)
            q()

        } else {

            #Get list of GCA survivors
            GCA_accessions <- names(GCA_accession_list)
            flog.info(paste("Found", paste0(GCA_accessions, collapse = ", ")))
            #Download, extract and concatenate into a single multifasta file.
            GCA_accessions_df <- data.frame(Assembly_accession = unlist(GCA_accession_list[GCA_accessions]), GCAfn = paste(GCA_accessions, "fasta.gz", sep = "_"))
            for (an in 1:nrow(GCA_accessions_df)){
                #Download
                tryCatch(download.file(GCA_accessions_df[an, "Assembly_accession"], destfile = GCA_accessions_df[an, "GCAfn"] ), error = function(e) print(paste("Unable to download", GCA_accessions_df[an, "Assembly_accession"])))
                #Concatenate
                system(paste("cat", GCA_accessions_df[an, "GCAfn"], ">>", paste0(opt$prefix, "_contigs.fasta.gz")))
            }
            #Expand
            system(paste("unpigz", paste0(opt$prefix, "_contigs.fasta.gz")))
        }

        #Set opt$contigs to NULL in case the user has tried to set that as an input too
        opt <- prepare_contigs_for_JAMS(opt = opt, fastafile = paste0(opt$prefix, "_contigs.fasta"))

        #Write NHcontigs to the system
        NHcontigsfastaout <- paste(opt$prefix, "contigs.fasta", sep="_")
        write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$workdir, NHcontigsfastaout))

        #Bank assembly files to project directory
        if (opt$workdir != opt$sampledir){
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
