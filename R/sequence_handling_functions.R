#' countfastq(fastqfile)
#' Wrapper for counting number of reads and number of bases of a fastq file in the system.
#' Returns a vector with read counts and base counts.
#' Caveat: Input fastq file must have exactly four lines per sequence. If your fastq files do not fit this criterion, you are totally bonkers.
#' @export

countfastq <- function(fastqfile){

    countargs <- c(fastqfile, "|", "paste", "-", "-", "-", "-", "|", "cut", "-f", "2", "|", "wc", "-lc")
    fastqstats <- system2('cat', args = countargs, stdout = TRUE, stderr = FALSE)
    fastqstats <- unlist(strsplit(fastqstats, split = " "))
    fastqstats <- fastqstats[which(fastqstats != "")]
    fastqstats <- as.numeric(fastqstats)

    return(fastqstats)
}


#' countfastq_files(fastqfiles = NULL, threads = NULL)
#' Wrapper for applying the countfastq() function to a vector of filenames using multiple threads.
#' Returns a dataframe with read counts and base counts for each input fastq file.
#' Caveat: Input fastq file must have exactly four lines per sequence. If your fastq files do not fit this criterion, you are totally bonkers.
#' @export

countfastq_files <- function(fastqfiles = NULL, threads = NULL){
    #fastqstatslist <- lapply(1:length(fastqfiles), function (x) { countfastq(fastqfiles[x]) })
    fastqstatslist <- mclapply(1:length(fastqfiles), function (x) { countfastq(fastqfiles[x]) }, mc.cores = threads)
    readcounts <- sapply(1:length(fastqstatslist), function (x) { fastqstatslist[[x]][1] })
    basecounts <- sapply(1:length(fastqstatslist), function (x) { fastqstatslist[[x]][2] })
    fastqstatsdf <- data.frame(Reads = fastqfiles, Count = readcounts, Bases = basecounts, stringsAsFactors = FALSE)
    fastqstatsdf$Readlength <- round((fastqstatsdf$Bases / fastqstatsdf$Count), 0)
    fastqstatsdf$Bases <- format(fastqstatsdf$Bases, scientific = FALSE)
    fastqstatsdf$Count <- format(fastqstatsdf$Count, scientific = FALSE)
    fastqstatsdf$Readlength <- format(fastqstatsdf$Readlength, scientific = FALSE)
    fastqstatsdf$Reads <- as.character(fastqstatsdf$Reads)
    fastqstatsdf$Count <- as.numeric(fastqstatsdf$Count)
    fastqstatsdf$Bases <- as.numeric(fastqstatsdf$Bases)
    fastqstatsdf$Readlength <- as.numeric(fastqstatsdf$Readlength)

    return(fastqstatsdf)
}


#' rename_sequences_consecutively(sequence=NULL, sequence=NULL, headerprefix="ctg")
#'
#' Renames sequences numerically consecutively
#' @export

rename_sequences_consecutively <- function(sequence = NULL, headerprefix = "ctg"){
    sequence_names <- paste(headerprefix, formatC(1:length(sequence), width = nchar(length(sequence)), flag = "0"), sep = "_")

    #Fix attribute when loaded using sequinr load.fasta with set.attributes being TRUE and ensure CAPS
    for(t in 1:length(sequence)){
        tmps <- sequence[[t]]
        attr(tmps, "name") <- sequence_names[t]
        tmps <- toupper(tmps)
        sequence[[t]] <- tmps
    }
    names(sequence) <- sequence_names

    return(sequence)
}


#' filter_sequence_by_length(sequence=NULL, minlength=0, maxlength=Inf)
#'
#' Filters sequences by size
#' @export

filter_sequence_by_length <- function(sequence = NULL, minlength = 0, maxlength = Inf){
    sequence_lengths <- lapply(1:length(sequence), function(x){ length(sequence[[x]]) })
    seqsIwant <- which(sequence_lengths > minlength & sequence_lengths < maxlength)
    sequence <- sequence[seqsIwant]

    return(sequence)
}

#' filter_sequence_by_name(input_sequences = NULL, sequencenames = NULL, keep = TRUE)
#'
#' Filters sequinr sequences by name either keeping or discarding specified sequence names
#' @export

filter_sequence_by_name <- function(input_sequences = NULL, sequencenames = NULL, keep = TRUE){
    input_seqnames <- names(input_sequences)

    if (!(is.null(sequencenames))){
        #Filter by read name
        if (keep == TRUE){
            seqtoKeep <- which(input_seqnames %in% sequencenames)

        } else {
            seqtoKeep <- which(!(input_seqnames %in% sequencenames))
            #Check that all sequences you asked to eliminate existed in the input file.
            if (!(all(sequencenames %in% input_seqnames))){
                flog.info("WARNING: Not all sequences you requested to exclude were available in input file.")
            }
        }
        output_sequences <- input_sequences[seqtoKeep]
    } else {
        output_sequences <- input_sequences
    }

    output_sequences <- make_sequences_caps(output_sequences)

    return(output_sequences)
}


#' make_sequences_caps(sequence = NULL)
#'
#' Transforms sequinr sequences to uppercase
#' @export

make_sequences_caps <- function(sequence = NULL){
    seqnames <- names(sequence)
    sequence <- lapply(1:length(sequence), function(x){ toupper(sequence[[x]]) })
    names(sequence) <- seqnames

    return(sequence)
}


#' get_sequence_stats
#'
#' JAMSalpha function
#' @export

get_sequence_stats<-function(sequences=NULL){
    sequencenames<-names(sequences)
    sequencelengths<-getLength(sequences)
    sequenceGCs<-sapply(1:length(sequences), function(x){ GC(sequences[[x]]) })

    seqstats<-data.frame(Sequence=sequencenames, Length=sequencelengths, GC=sequenceGCs)

    return(seqstats)
}
