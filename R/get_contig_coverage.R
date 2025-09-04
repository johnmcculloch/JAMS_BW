#' get_contig_coverage(opt = NULL, markduplicates = FALSE, align_as_unpaired_reads = TRUE)
#'
#' JAMSalpha function
#' @export

get_contig_coverage <- function(opt = NULL, markduplicates = FALSE, align_as_unpaired_reads = TRUE){
    setwd(opt$workdir)

    #Only get coverage from reads if contigs were actually assembled from reads.
    if ("readsdir" %in% names(opt)){
        #Check if reads are in tmpdir
        if(opt$workdir != opt$sampledir){
            readsworkdir <- file.path(opt$workdir, "reads")
        } else {
            readsworkdir <- opt$readsdir
        }
        setwd(readsworkdir)

        #Choose input reads, either NAHS or if not present, trimmed reads.
        if ((opt$host == "none") || (opt$analysis %in% c("isolate", "isolaternaseq"))){
            inputreads <- opt$trimreads
        } else {
            inputreads <- opt$nahsreads
        }
        #Adjust input reads to absolute path
        inputreads <- file.path(readsworkdir, inputreads)

        #Check the files exist or die.
        if (!(all(file.exists(inputreads)))){
            flog.info("Could not find input reads which were used for assembly. Aborting now.")
            q()
        }

        flog.info("Creating directory to calculate read coverage.")
        opt$covdir <- file.path(opt$workdir, paste(opt$prefix, "coverage", sep="_"))
        dir.create(opt$covdir, recursive = TRUE)
        setwd(opt$covdir)

        #Write a copy of contigs to be analyzed to the system
        flog.info("Writing contigs for read alignment.")
        write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = "contigscov.fa")
        #Copy annotation data
        file.copy(opt$bedfile, opt$covdir)

        #Set general options
        bowtiethreads <- opt$threads

        #build index
        flog.info("Building index.")
        system2('bowtie2-build', args = c("-f", "contigscov.fa", "contigscov"), stdout = TRUE, stderr = TRUE)

        #Build command from options present
        bowtiecommonargs <- c("-q", "--no-unal", "--threads", bowtiethreads, "-x", "contigscov")
        if (opt$seqtype == "illuminamp"){
            bowtiecommonargs <- c(bowtiecommonargs, "--rf", "--maxins", "10000")
        }
        bowtiereadargs <- list("-U", c("-1", "-2"), c("-1", "-2", "-U"))

        if (align_as_unpaired_reads){
            flog.info("Aligning reads back to contigs as unpaired reads.")
            inputreads <- paste0(inputreads, collapse = ",")
        }

        bowtiereadinput <- as.vector(rbind(bowtiereadargs[[length(inputreads)]], inputreads))

        bowtiereadoutput <- c("--un", (paste(opt$prefix, "NAss_SE.fastq", sep="_")), "--un-conc", (paste(opt$prefix, "NAss_PE.fastq", sep="_")))[1:(2 * min(length(inputreads), 2))]

        samfile <- paste(paste(opt$prefix, "contigs", sep = "_"), "sam", sep = ".")
        bowtiesam <- c("-S", samfile)
        bowtieargs <- c(bowtiecommonargs, bowtiereadinput, bowtiereadoutput, bowtiesam)

        flog.info("Aligning reads back to contigs.")
        opt$bowtie_cov_output <- system2('bowtie2', args = bowtieargs, stdout = TRUE, stderr = TRUE)
        flog.info(paste("Alignment of reads to contigs is complete with", opt$bowtie_cov_output[length(opt$bowtie_cov_output)]))
        opt$NAssfastqs <- list.files(pattern = "NAss_")

        #Just count NAssfastqs and get rid of them
        flog.info("Computing Non-Assembled (NAss) reads stats")
        opt$fastqstats <- rbind(opt$fastqstats, countfastq_files(fastqfiles = opt$NAssfastqs, threads = opt$threads, type_tag = "NAss"))
        file.remove(opt$NAssfastqs)

        #index
        flog.info("Indexing and sorting alignment.")
        system2('samtools', args = c("faidx", "contigscov.fa"), stdout = TRUE, stderr = TRUE)
        #convert to bam & sort
        bamfile <- paste(paste(opt$prefix, "contigs", sep="_"), "bam", sep=".")
        sortedbamfile <- paste(paste(opt$prefix, "contigs", sep="_"), "sorted", "bam", sep=".")
        sortedmarkdupbamfile <- paste(paste(opt$prefix, "contigs", sep="_"), "sorted", "markdup", "bam", sep=".")
        system2('samtools', args = c("view", "--threads", opt$threads, "-q", "10", "-F", "4", "-bt", "contigscov.fa.fai", samfile, ">", bamfile), stdout = TRUE, stderr = TRUE)
        system(paste("samtools", "sort", "--threads", opt$threads, bamfile, ">", sortedbamfile))
        system(paste("samtools", "index", sortedbamfile))

        #Marking duplicates with picard is now deprecated
        #if (markduplicates == TRUE){
        if (FALSE) {
            flog.info("Removing duplicate reads from alignment.")
            sortedmarkdupbamfile <- paste(paste(opt$prefix, "contigs", sep = "_"), "sorted", "markdup", "bam", sep = ".")
            sortedmarkdupsortedbamfile <- paste(paste(opt$prefix, "contigs", sep = "_"), "sorted", "markdup", "sorted", "bam", sep = ".")
            sortedmarkdupsortedflagstatfile <- paste(paste(opt$prefix, "contigs", sep = "_"), "sorted", "markdup", "sorted", "flagstat", sep = ".")
            sortedmarkdupsortedstatsfile <- paste(paste(opt$prefix, "contigs", sep = "_"), "sorted", "markdup", "sorted", "bam", "stats", sep = ".")

            picardcmd <- paste0("picard MarkDuplicates INPUT=", sortedbamfile, " OUTPUT=", sortedmarkdupbamfile, " METRICS_FILE=", paste(sortedmarkdupbamfile, "metrics", sep="."), " AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE READ_NAME_REGEX=null")
            system(picardcmd)
            system(paste("samtools", "sort", "--threads", opt$threads, sortedmarkdupbamfile, ">", sortedmarkdupsortedbamfile))
            system(paste("samtools", "index", sortedmarkdupsortedbamfile))
            system(paste("samtools", "flagstat", sortedmarkdupsortedbamfile, ">", sortedmarkdupsortedflagstatfile))
            system(paste("samtools", "stats", sortedmarkdupsortedbamfile, ">", sortedmarkdupsortedstatsfile))
            coveragebamfile <- sortedmarkdupsortedbamfile
        } else {
            coveragebamfile <- sortedbamfile
        }

        covdflist <- list()

        #Calculating contig base coverage
        flog.info("Calculating base coverage depth in contigs.")
        coverageargs <- c("-ibam", coveragebamfile, "-d", ">", "perbasecoverage.map")
        system2('genomeCoverageBed', args = coverageargs, stdout = TRUE, stderr = TRUE)
        covdflist$contigcoverage <- fread(file = "perbasecoverage.map", sep = "\t", header = FALSE, nThread = opt$threads)
        colnames(covdflist$contigcoverage) <- c("Feature", "Position", "Depth")
        covdflist$contigcoverage <- covdflist$contigcoverage[ , c("Feature", "Depth")]

        #Calculating feature coverage
        flog.info("Calculating coverage depth for features.")
        sortedbed <- paste0(opt$prefix, ".sorted.markdup.sorted.bed")
        bedcmd <- paste0("bamToBed -i ", coveragebamfile, " | sort --temporary-directory=", opt$workdir, " -k1,1 -k2,2n > ", sortedbed)
        system(bedcmd)
        bedfeatargs <- c("coverage", "-hist", "-sorted", "-b", sortedbed, "-a", opt$bedfile, ">", "featuredep.tsv")
        system2('bedtools', args = bedfeatargs, stdout = TRUE, stderr = FALSE)
        covdflist$featuredep <- fread(file = "featuredep.tsv", sep = "\t", stringsAsFactors = FALSE, fill = TRUE, nThread = opt$threads)
        colnames(covdflist$featuredep) <- c("Contig", "Start", "End", "Feature", "MapQual", "Strand", "Annotby", "FeatType", "Spin", "Annot", "Depth", "PartNumBasesofFeature", "LengthofFeat", "PctofAatdep")
        covdflist$featuredep <- subset(covdflist$featuredep, FeatType %in% c("CDS", "tRNA",  "rRNA",  "tmRNA"))
        covdflist$featuredep$PartDepthofFeature <- (covdflist$featuredep$Depth * covdflist$featuredep$PartNumBasesofFeature)
        covdflist$featuredep <- covdflist$featuredep[ , c("Feature", "PartDepthofFeature")]
        colnames(covdflist$featuredep) <- c("Feature", "Depth")

        #Aggregate in parallel
        sum_up_bases <- function(df = NULL){
            aggregate_df <- df[ , list(NumBases = sum(Depth)), by = "Feature"]

            return(aggregate_df)
        }

        if (opt$threads > 3){
            aggregate_df_list <- mclapply(names(covdflist), function(x) { sum_up_bases(df = covdflist[[x]])}, mc.cores = 4)
        } else {
            aggregate_df_list <- lapply(names(covdflist), function(x) { sum_up_bases(df = covdflist[[x]])})
        }

        names(aggregate_df_list) <- names(covdflist)
        colnames(aggregate_df_list$contigcoverage) <- c("Contig", "NumBases")
        opt$contigsdata <- left_join(opt$contigsdata, as.data.frame(aggregate_df_list$contigcoverage), by = "Contig")
        opt$contigsdata$NumBases <- as.numeric(opt$contigsdata$NumBases)
        opt$featuredata <- left_join(opt$featuredata, as.data.frame(aggregate_df_list$featuredep), by = "Feature")
        opt$featuredata$NumBases <- as.numeric(opt$featuredata$NumBases)

    } else {

        flog.info("Short reads were not used for contig assembly, so sequencing depth will be set at 1X for all contigs.")
        #if there are no reads, assume depth of 1, ie just repeat length of contig as NumBases.
        opt$contigsdata$NumBases <- opt$contigsdata$Length
        opt$featuredata$NumBases <- opt$featuredata$LengthDNA

    } #End conditional get doses from short read alignment

    #Enforce correct row names to ensure downstream efficiency
    rownames(opt$contigsdata) <- opt$contigsdata$Contig 
    rownames(opt$featuredata) <- opt$featuredata$Feature 

    return(opt)
}
