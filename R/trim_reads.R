#' trim_reads
#'
#' JAMSalpha function
#' @export

trim_reads <- function(opt = NULL, discardleftoverSE = FALSE, qual = 18, autodetect_phred_offest = FALSE){

    #Check if reads are in tmpdir
    if (opt$workdir != opt$sampledir){
        readsworkdir <- file.path(opt$workdir, "reads")
    } else {
        readsworkdir <- opt$readsdir
    }
    setwd(readsworkdir)

    if (!(all(file.exists(file.path(readsworkdir, opt$rawreads))))){
        flog.info("Could not find raw reads to trim. Aborting now")
        q()
    }

    opt$phredoffset <- 33

    if (opt$seqtype %in% c("illuminape", "illuminamp")){

        #Write copy of Illumina adapters to system
        data(seqadapters)
        write.fasta(sequences = seqadapters, names = names(seqadapters), nbchar = 80, file.out = file.path(readsworkdir, "adapter.fasta"))
        adapters <- file.path(readsworkdir, "adapter.fasta")

        flog.info("Sequences are Illumina reads. Clipping adapters and quality trimming fastq reads.")
        #Trimmomatic options
        sliding=4
        crop=0
        minlen=36

        trimmcommand <- "trimmomatic"
        #Deal with JAVA memory (safer in newer versions of trimmomatic)
        mem_gb <- floor(opt$totmembytes * 0.9 / 1024^3)
        java_opts <- paste0("-Xmx", mem_gb, "G")

        if (autodetect_phred_offest == TRUE){
            scores <- system2('head', args = c("-1000", opt$rawreads[1]), stdout = TRUE)
            scores <- scores[seq(4, 1000, by = 4)]
            atomized_scores <- unlist(strsplit(scores, split = ""))
            atomized_scores_values <- sapply(atomized_scores, function(x) { utf8ToInt(x) } )
            atomized_scores_values <- atomized_scores_values - 33
            scores_quantiles <- quantile(atomized_scores_values)
            if (unname(scores_quantiles["25%"]) > 40){
                opt$phredoffset <- 64
            } else {
                opt$phredoffset <- 33
            }
            flog.info(paste("Looks like the fastqs have a Phred offset of", opt$phredoffset))
        }

        if (opt$libstructure == "pairedend"){
            libstruct <- "PE"
            input.read1 <- opt$rawreads[1]
            input.read2 <- opt$rawreads[2]
            output.trim1 <- paste(opt$prefix, libstruct, "R1_trim.fastq.gz", sep="_")
            output.trim2 <- paste(opt$prefix, libstruct, "R2_trim.fastq.gz", sep="_")
            output.trim_unpaired1 <- paste(opt$prefix, libstruct, "R1_trim_unpaired.fastq.gz", sep="_")
            output.trim_unpaired2 <- paste(opt$prefix, libstruct, "R2_trim_unpaired.fastq.gz", sep="_")
            #filter reads
            argstorun <- paste(libstruct, "-threads", opt$threads, input.read1, input.read2, output.trim1, output.trim_unpaired1, output.trim2, output.trim_unpaired2, paste0("ILLUMINACLIP:", adapters, ":2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:", sliding, ":", qual), paste0("-phred", as.character(opt$phredoffset)), paste0("HEADCROP:",crop), paste0("MINLEN:", minlen), sep=" ")
            system2(command = trimmcommand, args = argstorun, env = c(JAVA_OPTS = java_opts))

            #Since JAMS > 2.0, use of unpaired trimmed reads is deprecated.
            #merge unpaired
            #output.trimSE <- paste(opt$prefix, "SE_trim.fastq.gz", sep="_")
            #commandtorun <- paste("cat *_unpaired.fastq.gz >>", output.trimSE)
            #system(commandtorun)
            #delete unpaired1 and unpaired2
            file.remove(c(output.trim_unpaired1, output.trim_unpaired2))

            #if (all(c(discardleftoverSE==TRUE), (file.info(output.trimSE)$size > 1))){
                #If leftover single unpaied, only add if smaller than 250 Mb. Otherwise sequencing errors may seep into the contigs and assembly becomes much more complicated.
                #SEsize <- file.size(output.trimSE)
                #if (SEsize < 250000000){
                    #opt$trimreads <- c(opt$trimreads, output.trimSE)
                #} else {
                    #flog.info(paste("There are leftover unpaired reads from trimming (when only one pair survives trimming), but as the file size is large, these will be ignored for read assembly otherwise sequencing errors may seep into contigs."))
                #}
            #} else {
                #opt$trimreads<-c(opt$trimreads, output.trimSE)
            #}

            #add output of what was trimmed to opt
            opt$trimreads <- c(output.trim1, output.trim2)

        } else if (opt$libstructure == "singleend"){
            libstruct <- "SE"
            input.read1 <- opt$rawreads[1]
            output.trimSE <- paste(opt$prefix, libstruct, "trim.fastq.gz", sep="_")

            #filter reads
            argstorun <- paste(libstruct, "-threads", opt$threads, input.read1, output.trimSE, paste0("ILLUMINACLIP:", adapters, ":2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:", sliding, ":", qual), paste0("-phred", as.character(opt$phredoffset)), paste0("HEADCROP:",crop), paste0("MINLEN:", minlen), sep=" ")
            system2(command = trimmcommand, args = argstorun, env = c(JAVA_OPTS = java_opts))

            #add info of what was acieved to opt
            opt$trimreads <- output.trimSE
        }
        file.remove("adapter.fasta")

    } else if (opt$seqtype == "iontorrent") {
        if (opt$libstructure == "singleend"){
            flog.info("Sequences are single-end Ion Torrent reads. Using IonHammer (SPAdes) for fastq read correction.")
            libstruct <- "SE"
            input.read1 <- opt$rawreads[1]
            output.trimSE <- paste(opt$prefix, libstruct, "trim.fastq.gz", sep="_")

            #filter reads
            commandtorun <- paste("spades.py", "-t", opt$threads, "-s", input.read1, "-o", "corrreads", "--only-error-correction", "--disable-gzip-output", "1>", "corrreads.log", "2>", "corrreads.err")
            system(commandtorun)
            corrfile <- list.files(path=file.path(readsworkdir, "corrreads", "corrected"), pattern=".fastq.gz")[1]
            file.copy(file.path(readsworkdir, "corrreads", "corrected", corrfile), file.path(readsworkdir, output.trimSE))
            unlink(file.path(readsworkdir, "corrreads", "corrected"))
            #add info of what was acieved to opt
            opt$trimreads <- output.trimSE
         } else {
            flog.info("Sequences are paired-end Ion Torrent reads. Defaulting the assembler to SPAdes, which will correct the reads with IonHammer")
            opt$assembler <- "spades"
            libstruct <- "PE"
            input.read1 <- opt$rawreads[1]
            input.read2 <- opt$rawreads[2]
            output.trim1 <- paste(opt$prefix, libstruct, "R1_trim.fastq.gz", sep="_")
            output.trim2 <- paste(opt$prefix, libstruct, "R2_trim.fastq.gz", sep="_")
            file.copy(input.read1, output.trim1)
            file.copy(input.read2, output.trim2)
            opt$trimreads <- c(output.trim1, output.trim2)
         }

    } #End conditionals for sequencing  platform type

    #Disregard unpaired reads from trimming if not assembling an isolate.
    if (opt$analysis != "isolate"){
        if (length(opt$trimreads) > 2){
            opt$trimreads <- opt$trimreads[1:2]
        }
    }

    #Collect information on trimmed reads and purge raw reads to conserve space
    file.remove(opt$rawreads)
    flog.info("Computing trimmed reads")
    opt$fastqstats <- rbind(opt$fastqstats, countfastq_files(fastqfiles = opt$trimreads, threads = opt$threads, type_tag = "Trimmed"))

    if (opt$workdir != opt$sampledir){
        #Bank trimmed reads to project directory
        file.copy(opt$trimreads, opt$readsdir)
    }

    return(opt)
}
