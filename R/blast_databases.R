#' blast_databases
#'
#' JAMSalpha function
#' @export

blast_databases <- function(opt = NULL, blastanalyses = NULL, QcovThreshold = 75, HitcovThreshold = 75, PidentThreshold = 75){

    #Set working directory to workdir
    setwd(opt$workdir)

    #Write predicted proteome to serve as query.
    flog.info("Writing predicted proteome to serve as query.")
    write.fasta(sequences = opt$proteome, names = names(opt$proteome), nbchar = 80, file.out = "proteome.faa")
    flog.info("Writing predicted genes to serve as query.")
    write.fasta(sequences = opt$genes, names = names(opt$genes), nbchar = 80, file.out = "genes.fna")

    if (is.null(blastanalyses)){
        blastanalyses <- opt$blastanalyses
    }

    blastanalyses <- blastanalyses[(blastanalyses %in% opt$blastanalyses)]

    for (blastanalysis in blastanalyses){
        flog.info(paste("Blasting genes against", blastanalysis, "database"))

        #Copy database to tmpdir each time for speed.
        blastdbori <- file.path(opt$blastpath, blastanalysis)
        if (file.exists(blastdbori)){
            dir.create("tmpblastdb")
            file.copy(file.path(blastdbori, list.files(path=blastdbori, pattern="sequen")), "tmpblastdb")
            blastdb <- file.path("tmpblastdb", "sequences")

            #Load lookup table if it exists
            if (file.exists(file.path(blastdbori, "lookup.tsv"))){
                lookup <- read.table(file = file.path(blastdbori, "lookup.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
            } else {
                #Try the newer rda lookups
                lookuprdaobject <- paste(blastanalysis, "lookup", sep = "_")
                rdalookup <- file.path(blastdbori, paste(lookuprdaobject, "rda", sep = "."))
                if (file.exists(rdalookup)){
                    load(rdalookup)
                    lookup <- get(lookuprdaobject)
                } else {
                    lookup <- NULL
                }
            }

            #Find out if db is protein or DNA
            if (length(grep("nhr", list.files(blastdbori))) > 0){
                dbtype <- "nucl"
                blastargs <- c("-num_threads", opt$threads, "-query", "genes.fna", "-db", blastdb, "-outfmt", "'6 qseqid qstart qend qlen sseqid sstart send slen evalue length pident gaps'", "-evalue", "1E-20", "-max_target_seqs", "10000", "-culling_limit", "1", ">", "blast.out")
                blastprog <- "blastn"
            } else {
                dbtype <- "prot"
                blastargs <- c("-num_threads", opt$threads, "-query", "proteome.faa", "-db", blastdb, "-outfmt", "'6 qseqid qstart qend qlen sseqid sstart send slen evalue length pident gaps'", "-evalue", "1E-20", "-max_target_seqs", "10000", "-culling_limit", "1", ">", "blast.out")
                blastprog <- "blastp"
            }
            #Blast query sequences against database
            system2(blastprog, args = blastargs, stdout = TRUE, stderr = TRUE)
            numlinblast <- as.numeric(system2('cat', args=c("blast.out", "|", "wc", "-l"), stdout = TRUE, stderr = FALSE))
            if (numlinblast > 0){
                blastout <- fread(file = "blast.out", header = FALSE, sep = "\t", colClasses = c("character", "numeric", "numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric"))
                colnames(blastout) <- c("Feature", "qstart", "qend", "qlen", "Accession", "sstart", "send", "slen", "evalue", "length", "Pident", "gaps")
                file.remove("blast.out")

                #Cleanup
                blastout$Qcov <- round(((blastout$qend - blastout$qstart)/blastout$qlen) * 100, 2)
                blastout$Hitcov <- round(((blastout$send - blastout$sstart)/blastout$slen) * 100, 2)
                blastout <- blastout[, c("Feature", "Accession", "Qcov", "Hitcov", "Pident", "gaps", "evalue")]
                blastout <- subset(blastout, Qcov > QcovThreshold)
                blastout <- subset(blastout, Hitcov > HitcovThreshold)
                blastout <- subset(blastout, Pident > PidentThreshold)
                if (nrow(blastout) > 0){
                    blastout <- left_join(blastout, opt$featuredata, by = "Feature")
                    blastout$LKT <- opt$contigsdata$LKT[match(blastout$Contig, opt$contigsdata$Contig)]
                    if (!is.null(lookup)){
                        lookup$Feature <- NULL
                        blastout <- left_join(blastout, lookup, by = "Accession")
                        blastout[is.na(blastout[])] <- "none"
                    }
                    nextopt <- (length(names(opt)) + 1)
                    opt[[nextopt]] <- as.data.frame(blastout)
                    names(opt)[nextopt] <- blastanalysis
               } else {
                    flog.info(paste("No significant hits found for", blastanalysis))
                }
            } else {
                flog.info(paste("No significant hits found for", blastanalysis))
            }
            #Delete temporary Blast database
            unlink("tmpblastdb", recursive = TRUE)
        } else {
            flog.info(paste(blastanalysis, "database does not exist. Skipping this blast analysis."))
        }
    }
    file.remove("proteome.faa")
    file.remove("genes.fna")

    return(opt)
}
