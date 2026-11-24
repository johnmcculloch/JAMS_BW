#' describe_JAMSalpha_run(opt = NULL)
#'
#' Internal helper for generate_JAMSalpha_methods() and tabulate_JAMSalpha_methods().
#' Distils a JAMSalpha opt (live at the end of a run, restored from a sample's workspace image,
#' or rebuilt from a .jams file with rebuild_JAMSalpha_opt_from_jams()) into a flat list of
#' methods-relevant facts. Both the prose and the table read from this single list, so they can
#' never disagree. Anything that cannot be determined from opt is left NA/NULL, never guessed;
#' the sentence or row depending on it is then omitted.
#' @keywords internal

describe_JAMSalpha_run <- function(opt = NULL){

    if (is.null(opt) || !is.list(opt)){
        stop("opt must be a JAMSalpha opt list.")
    }

    # ---- helpers -------------------------------------------------------------
    #Always use [[ ]] (exact matching). opt$x partially matches, so e.g. opt$methods would
    #silently return opt$methods_table on a run where only the latter exists.
    g   <- function(x){ opt[[x]] }
    nz  <- function(x){ !is.null(x) && is.atomic(x) && length(x) > 0 && !all(is.na(x)) && !all(as.character(x) == "", na.rm = TRUE) }
    chr <- function(x){ if (nz(x)) as.character(x)[1] else NA_character_ }
    num <- function(x){ if (nz(x)) suppressWarnings(as.numeric(x)[1]) else NA_real_ }
    isdf <- function(x){ is.data.frame(x) && nrow(x) > 0 }
    bacarch <- c("d__2157_Archaea", "d__Archaea", "d__2_Bacteria", "d__Bacteria")
    quality_order <- c("HQ", "MHQ", "MQ", "LQ", "Contaminated")
    tally_quality <- function(df){
        if (!isdf(df) || !("Quality" %in% colnames(df))) return(NULL)
        tq <- table(factor(df$Quality, levels = quality_order))
        setNames(as.integer(tq), names(tq))
    }

    f <- list()
    f$rebuilt_from_jams <- isTRUE(g("rebuilt_from_jams"))
    #The procedures described are those of the installed JAMS code. Flag runs made with another version,
    #and do not read "object absent" as "zero" for them (the object may simply not have existed then).
    f$JAMS_version <- if (nz(g("verstr"))) sub("^JAMS v", "", chr(g("verstr"))) else NA_character_
    f$installed_JAMS_version <- tryCatch(as.character(packageVersion("JAMS")), error = function(e) NA_character_)
    f$version_matches_installed <- !is.na(f$JAMS_version) && identical(f$JAMS_version, f$installed_JAMS_version)
    absent_means_zero <- !f$rebuilt_from_jams && f$version_matches_installed

    # ---- identity and input --------------------------------------------------
    f$prefix            <- chr(g("prefix"))
    f$analysis          <- chr(g("analysis"))
    f$seqtype           <- chr(g("seqtype"))
    f$libstructure      <- chr(g("libstructure"))
    f$readorigin        <- chr(g("readorigin"))
    f$sraaccession      <- chr(g("sraaccession"))
    f$assemblyaccession <- if (nz(g("assemblyaccession"))) paste0(unlist(strsplit(as.character(g("assemblyaccession")), split = ",")), collapse = ", ") else NA_character_
    f$process           <- chr(g("process"))
    f$from_reads <- nz(g("readsdir")) || identical(f$process, "Assemble_from_reads") || isdf(g("readstats"))

    if (f$from_reads){
        f$input_type <- "reads"
    } else if (!is.na(f$assemblyaccession)){
        f$input_type <- "assemblyaccession"
    } else if (nz(g("contigsfasta")) || nz(g("contigs")) || identical(f$process, "Contigs_supplied")){
        f$input_type <- "contigs"
    } else {
        f$input_type <- NA_character_
    }

    # ---- host ----------------------------------------------------------------
    f$host                  <- chr(g("host"))
    f$hostspecies           <- chr(g("hostspecies"))
    f$host_accession_number <- chr(g("host_accession_number"))
    f$host_filtered <- f$from_reads && !is.na(f$host) && (f$host != "none")
    #assemble_contigs() and get_contig_coverage() use trimmed reads when there is no host, or for
    #isolates/isolate RNA-seq, and host-depleted (NAHS) reads otherwise.
    f$assembly_input <- if (is.na(f$host) || (f$host == "none") || (f$analysis %in% c("isolate", "isolaternaseq"))) "Trimmed" else "NAHS"

    # ---- read statistics -----------------------------------------------------
    fq <- g("fastqstats")
    readrows <- function(type, files){
        if (!isdf(fq)) return(NULL)
        if ("Read_Type" %in% colnames(fq)){
            r <- fq[which(fq$Read_Type == type), , drop = FALSE]
        } else if (nz(files)){
            r <- fq[which(fq$Reads %in% files), , drop = FALSE]
        } else {
            return(NULL)
        }
        if (nrow(r) == 0) NULL else r
    }
    rawrows  <- readrows("Raw",     g("rawreads"))
    trimrows <- readrows("Trimmed", g("trimreads"))
    nahsrows <- readrows("NAHS",    g("nahsreads"))

    f$raw_reads <- f$raw_readlength <- NA_real_
    if (!is.null(rawrows)){
        f$raw_reads      <- sum(as.numeric(rawrows$Count), na.rm = TRUE)
        rb <- sum(as.numeric(rawrows$Bases), na.rm = TRUE)
        f$raw_readlength <- if (f$raw_reads > 0) round(rb / f$raw_reads) else NA_real_
    }
    f$raw_pairs <- if (identical(f$libstructure, "pairedend") && !is.na(f$raw_reads)) f$raw_reads / 2 else NA_real_

    rs <- g("readstats")
    rsv <- function(type){
        if (isdf(rs) && all(c("Read_type", "Base_counts") %in% colnames(rs)) && (type %in% rs$Read_type)){
            as.numeric(rs$Base_counts[match(type, rs$Read_type)])
        } else {
            NA_real_
        }
    }
    f$bases_raw        <- rsv("Raw")
    f$bases_trimmed    <- rsv("Trimmed")
    f$bases_nonhost    <- rsv("NonHost")
    f$bases_subsampled <- rsv("Subsampled")
    f$bases_assembled  <- rsv("Assembled")
    #Fall back on fastqstats where readstats is missing
    if (is.na(f$bases_raw) && !is.null(rawrows))      f$bases_raw     <- sum(as.numeric(rawrows$Bases),  na.rm = TRUE)
    if (is.na(f$bases_trimmed) && !is.null(trimrows)) f$bases_trimmed <- sum(as.numeric(trimrows$Bases), na.rm = TRUE)
    if (is.na(f$bases_nonhost) && !is.null(nahsrows)) f$bases_nonhost <- sum(as.numeric(nahsrows$Bases), na.rm = TRUE)

    f$pct_bases_trimmed_kept <- if (!is.na(f$bases_raw) && !is.na(f$bases_trimmed) && f$bases_raw > 0) 100 * f$bases_trimmed / f$bases_raw else NA_real_

    f$pct_host <- NA_real_
    if (f$host_filtered){
        if (!is.na(f$bases_trimmed) && !is.na(f$bases_nonhost) && f$bases_trimmed > 0){
            f$pct_host <- 100 * (1 - (f$bases_nonhost / f$bases_trimmed))
            f$pct_host_basis <- "bases"
        } else if (nz(g("bowtie_host_output"))){
            bho <- grep("overall alignment rate", as.character(g("bowtie_host_output")), value = TRUE)
            if (length(bho) > 0){
                f$pct_host <- suppressWarnings(as.numeric(sub("%.*$", "", trimws(tail(bho, 1)))))
                f$pct_host_basis <- "reads"
            }
        }
    }

    # ---- subsampling ---------------------------------------------------------
    f$maxbases_input <- num(g("maxbases_input"))
    tbs <- g("totbasesbeforesubsampling")
    f$subsampled <- isdf(tbs) || !is.na(f$bases_subsampled)
    #n.b. the readstats "Subsampled" (and therefore "Assembled") values cannot be trusted for subsampled runs:
    #compute_readcounts() takes them from fastqstats rows counted BEFORE seqtk overwrote the read files, so
    #"Subsampled" equals the pre-subsampling count. Use the subsampling target set by assemble_contigs()
    #instead (an approximation, flagged as such) and never the readstats value.
    f$bases_subsampled <- if (isdf(tbs) && ("MaxBasesAllowed" %in% colnames(tbs))) sum(as.numeric(tbs$MaxBasesAllowed), na.rm = TRUE) else NA_real_

    #Bases offered to the assembler (and later aligned back to contigs)
    if (f$subsampled){
        f$bases_assembly_input <- f$bases_subsampled
    } else if (f$assembly_input == "NAHS"){
        f$bases_assembly_input <- f$bases_nonhost
    } else {
        f$bases_assembly_input <- f$bases_trimmed
    }
    #Assembled bases = assembly input minus bases in reads that did not align back to the contigs (NAss).
    #Recomputed here rather than read from readstats, which is wrong for subsampled runs (above) and for
    #isolate RNA-seq with a host (compute_readcounts() subtracts from NonHost although trimmed reads were used).
    nassrows <- readrows("NAss", g("NAssfastqs"))
    f$bases_unassembled <- if (!is.null(nassrows)) sum(as.numeric(nassrows$Bases), na.rm = TRUE) else NA_real_
    if (!is.na(f$bases_unassembled) && !is.na(f$bases_assembly_input)){
        f$bases_assembled <- f$bases_assembly_input - f$bases_unassembled
    } else if (f$subsampled || (identical(f$analysis, "isolaternaseq") && f$host_filtered)){
        f$bases_assembled <- NA_real_
    }
    f$pct_assembled <- if (!is.na(f$bases_assembled) && !is.na(f$bases_assembly_input) && f$bases_assembly_input > 0) 100 * f$bases_assembled / f$bases_assembly_input else NA_real_

    # ---- assembly ------------------------------------------------------------
    f$contig_minlength <- 500
    f$assembler <- f$assembler_command <- NA_character_
    f$estcoverage <- NA_real_
    f$megahit <- NULL
    if (f$input_type %in% "reads"){
        #n.b. opt$assembler is set to "megahit" by default even when contigs are supplied, so only
        #read it for runs that actually assembled reads.
        f$assembler         <- chr(g("assembler"))
        f$assembler_command <- chr(g("assembler_command"))
        f$estcoverage       <- num(g("estcoverage"))
        f$deepseq           <- if (nz(g("deepseq"))) isTRUE(as.logical(g("deepseq"))) else NA

        if (identical(f$assembler, "megahit")){
            grab <- function(flag){
                if (is.na(f$assembler_command)) return(NA_real_)
                m <- regmatches(f$assembler_command, regexec(paste0("--", flag, " ([0-9]+)"), f$assembler_command))[[1]]
                if (length(m) == 2) as.numeric(m[2]) else NA_real_
            }
            mh <- list(kmin = grab("k-min"), kmax = grab("k-max"), kstep = grab("k-step"), mincount = grab("min-count"), source = "recorded_command")
            #Runs predating opt$assembler_command: re-derive the values exactly as assemble_contigs() does.
            if (is.na(mh$kmin) && !is.na(f$estcoverage)){
                if (f$estcoverage <= 5){
                    mh$kmin <- 21; mh$kstep <- 8;  mh$mincount <- 1
                } else if (f$estcoverage <= 40){
                    mh$kmin <- 23; mh$kstep <- 10; mh$mincount <- 2
                } else {
                    mh$kmin <- 27; mh$kstep <- 12; mh$mincount <- 3
                }
                mh$source <- "reconstructed_from_estcoverage"
            }
            if (is.na(mh$kmax)){
                inrows <- if (f$assembly_input == "NAHS") nahsrows else trimrows
                if (!is.null(inrows)){
                    inrows <- inrows[order(as.numeric(inrows$Count), decreasing = TRUE), , drop = FALSE]
                    estreadlength <- round(mean(as.numeric(inrows$Readlength[1:min(nrow(inrows), 2)])), 0)
                    if (!is.na(estreadlength)){
                        mh$kmax <- if (as.logical(estreadlength %% 2)) min(estreadlength, 255) else min((estreadlength - 1), 255)
                        mh$source <- if (identical(mh$source, "recorded_command")) "kmax_reconstructed_from_fastqstats" else paste(mh$source, "kmax_reconstructed_from_fastqstats", sep = "; ")
                    }
                }
            }
            if (all(is.na(unlist(mh[c("kmin", "kmax", "kstep", "mincount")])))) mh$source <- NA_character_
            f$megahit <- mh
        }
    }

    # ---- contigs -------------------------------------------------------------
    cd <- g("contigsdata")
    f$n_contigs <- f$contigs_total_length <- f$contigs_N50 <- f$contigs_max <- NA_real_
    if (isdf(cd) && ("Length" %in% colnames(cd))){
        L <- sort(as.numeric(cd$Length), decreasing = TRUE)
        L <- L[!is.na(L)]
        if (length(L) > 0){
            f$n_contigs            <- length(L)
            f$contigs_total_length <- sum(L)
            f$contigs_N50          <- L[which(cumsum(L) >= (sum(L) / 2))[1]]
            f$contigs_max          <- L[1]
        }
    }
    #prepare_contigs_for_JAMS() only creates opt$hostcontigsdata when it removes something, so on a
    #live/restored opt absence means zero. On a rebuilt opt (not banked to .jams) it means unknown.
    if (isdf(g("hostcontigsdata"))){
        f$n_metazoan_contigs <- nrow(g("hostcontigsdata"))
    } else {
        f$n_metazoan_contigs <- if (absent_means_zero) 0 else NA_real_
    }

    # ---- taxonomic classification --------------------------------------------
    f$kraken_confidence <- num(g("krakenconfidencescore"))
    f$kdb_version <- chr(g("JAMS_Kdb_Version"))
    if (identical(f$kdb_version, "UNDETERMINED")) f$kdb_version <- NA_character_
    f$checkm_db <- basename_safe(chr(g("CheckMdb")))

    # ---- denoising -----------------------------------------------------------
    f$denoise_applied <- identical(f$analysis, "metagenome") && f$from_reads
    if (f$denoise_applied){
        if (isdf(g("unvetable_contigs_before_censoring"))){
            f$n_censored_contigs <- nrow(g("unvetable_contigs_before_censoring"))
        } else {
            f$n_censored_contigs <- if (absent_means_zero) 0 else NA_real_
        }
    } else {
        f$n_censored_contigs <- NA_real_
    }

    # ---- annotation ----------------------------------------------------------
    fd <- g("featuredata")
    f$n_CDS <- f$n_rRNA <- f$n_tRNA <- NA_real_
    if (isdf(fd) && ("FeatType" %in% colnames(fd))){
        f$n_CDS  <- sum(fd$FeatType == "CDS",  na.rm = TRUE)
        f$n_rRNA <- sum(fd$FeatType == "rRNA", na.rm = TRUE)
        f$n_tRNA <- sum(fd$FeatType == "tRNA", na.rm = TRUE)
    }

    # ---- taxonomic abundance tables / bins -----------------------------------
    ab <- g("abundances")
    tx <- if (is.list(ab)) ab[["taxonomic"]] else NULL
    f$tax_spaces <- if (is.list(tx) && length(tx) > 0) names(tx)[vapply(names(tx), function(x){ isdf(tx[[x]]) }, logical(1))] else character(0)

    f$n_contig_LKTs <- if ("Contig_LKT" %in% f$tax_spaces) nrow(tx[["Contig_LKT"]]) else NA_real_

    f$mags_present <- "MB2bin" %in% f$tax_spaces
    f$binning_applied <- f$mags_present || (identical(f$analysis, "metagenome") && f$from_reads)
    f$n_mags <- if (f$mags_present) nrow(tx[["MB2bin"]]) else NA_real_
    f$mag_quality <- if (f$mags_present) tally_quality(tx[["MB2bin"]]) else NULL

    f$cgb_present <- "ConsolidatedGenomeBin" %in% f$tax_spaces
    f$n_cgb <- f$n_cgb_from_mb2 <- f$pct_depth_mb2 <- NA_real_
    f$cgb_quality <- NULL
    f$isolate_genome <- NULL
    if (f$cgb_present){
        cgb <- tx[["ConsolidatedGenomeBin"]]
        if (identical(f$analysis, "isolate")){
            pick <- function(col){ if (col %in% colnames(cgb)) as.character(cgb[1, col]) else NA_character_ }
            f$isolate_genome <- list(LKT = pick("LKT"), Completeness = suppressWarnings(as.numeric(pick("Completeness"))),
                                     Contamination = suppressWarnings(as.numeric(pick("Contamination"))), Quality = pick("Quality"),
                                     Genome_Size = suppressWarnings(as.numeric(pick("Genome_Size"))))
        } else {
            f$n_cgb <- nrow(cgb)
            f$cgb_quality <- tally_quality(cgb)
            cgbnames <- if ("ConsolidatedGenomeBin" %in% colnames(cgb)) cgb$ConsolidatedGenomeBin else rownames(cgb)
            f$n_cgb_from_mb2 <- sum(grepl("^MB2__", cgbnames))
            if (isdf(cd) && all(c("ConsolidatedGenomeBin", "NumBases") %in% colnames(cd))){
                nb <- as.numeric(cd$NumBases)
                if (sum(nb, na.rm = TRUE) > 0){
                    f$pct_depth_mb2 <- 100 * sum(nb[grepl("^MB2__", cd$ConsolidatedGenomeBin)], na.rm = TRUE) / sum(nb, na.rm = TRUE)
                }
            }
        }
    }
    f$n_nonbacterial_contigs <- NA_real_
    if (identical(f$analysis, "isolate") && isdf(cd) && ("Domain" %in% colnames(cd))){
        f$n_nonbacterial_contigs <- sum(!(cd$Domain %in% bacarch))
    }

    # ---- BLAST ---------------------------------------------------------------
    f$blastanalyses <- if (nz(g("blastanalyses"))) unique(as.character(g("blastanalyses"))) else character(0)
    f$blast_hits <- vapply(f$blastanalyses, function(b){
        x <- opt[[b]]
        if (is.data.frame(x) && ("Feature" %in% colnames(x))) length(unique(x$Feature)) else 0L
    }, integer(1))
    names(f$blast_hits) <- f$blastanalyses

    # ---- functional ----------------------------------------------------------
    fx <- if (is.list(ab)) ab[["functional"]] else NULL
    f$func_spaces <- if (is.list(fx) && length(fx) > 0) names(fx)[vapply(names(fx), function(x){ isdf(fx[[x]]) }, logical(1))] else character(0)
    func_analyses <- unique(unlist(lapply(f$func_spaces, function(x){ if ("Analysis" %in% colnames(fx[[x]])) unique(as.character(fx[[x]]$Analysis)) else NULL })))
    #Default noninterproanalyses of harvest_functions()
    basic_analyses <- c("FeatType", "ECNumber", "Product", "resfinder", "plasmidfinder", "napdos", "serofinderH", "serofinderO", "vfdb", "abricate")
    f$func_analyses <- func_analyses
    f$func_basic_analyses <- intersect(c("FeatType", "Product", "ECNumber"), func_analyses)
    ipro <- g("interproscanoutput")
    if (isdf(ipro) && ("Analysis" %in% colnames(ipro))){
        f$ipro_analyses <- sort(unique(as.character(ipro$Analysis)))
        if ("IproAcc" %in% colnames(ipro)) f$ipro_analyses <- c(f$ipro_analyses, "Interpro")
        if ("GOterms" %in% colnames(ipro)) f$ipro_analyses <- c(f$ipro_analyses, "GO")
    } else {
        f$ipro_analyses <- setdiff(func_analyses, c(basic_analyses, f$blastanalyses))
    }
    f$ipro_ran <- length(f$ipro_analyses) > 0
    f$useipro <- if (nz(g("useipro"))) isTRUE(as.logical(g("useipro"))) else NA

    # ---- software versions ---------------------------------------------------
    dv <- g("dependency_versions")
    depver <- function(dep){
        if (isdf(dv) && all(c("Dependency", "Version") %in% colnames(dv)) && (dep %in% dv$Dependency)){
            r <- dv[match(dep, dv$Dependency), , drop = FALSE]
            v <- as.character(r$Version)
            found <- if ("Found" %in% colnames(dv)) isTRUE(as.logical(r$Found)) else TRUE
            #Only accept something that looks like a version; get_dep_version() can fall back on a usage banner.
            if (found && nz(v) && grepl("^[vV]?[0-9]+\\.[0-9]", v)) return(sub("^[vV]", "", v))
        }
        NA_character_
    }

    used <- c("kraken2", "prokka", "convert2bed")
    if (f$input_type %in% "reads"){
        if (identical(f$readorigin, "sraaccession"))              used <- c(used, "sratoolkit")
        if (f$seqtype %in% c("illuminape", "illuminamp"))         used <- c(used, "trimmomatic")
        if (f$subsampled)                                         used <- c(used, "seqtk")
        used <- c(used, "bowtie2", "samtools", "bedtools")
        if (!is.na(f$assembler))                                  used <- c(used, f$assembler)
        if (identical(f$seqtype, "iontorrent") && identical(f$libstructure, "singleend")) used <- c(used, "spades")
    }
    if (f$binning_applied)                                        used <- c(used, "metabat2")
    if (f$analysis %in% c("metagenome", "isolate") || f$mags_present || f$cgb_present) used <- c(used, "checkm2")
    if (length(f$blastanalyses) > 0)                              used <- c(used, "blast")
    used <- unique(used)

    vers <- setNames(vapply(used, depver, character(1)), used)
    vsrc <- setNames(ifelse(is.na(vers), NA_character_, "dependency_versions"), used)
    #Older runs: the assembler's own --version banner was banked to opt$assemblerversion
    if (!is.na(f$assembler) && is.na(vers[f$assembler]) && nz(g("assemblerversion"))){
        av <- paste(as.character(g("assemblerversion")), collapse = " ")
        m <- regmatches(av, regexpr("[0-9]+\\.[0-9][0-9A-Za-z\\._\\-]*", av))
        if (length(m) == 1 && nzchar(m)){
            vers[f$assembler] <- m
            vsrc[f$assembler] <- "assemblerversion"
        }
    }
    f$versions <- vers
    f$version_source <- vsrc
    f$tools_used <- used

    #f$JAMS_version (set above) is the version of the run itself. Deliberately NOT packageVersion("JAMS"):
    #when describing an old run retrospectively, the installed version is not the one that produced the data.
    f$R_version    <- chr(g("R_version_string"))
    f$threads      <- num(g("threads"))
    fmt_time <- function(x){ if (inherits(x, "POSIXt")) format(x, "%Y-%m-%d %H:%M:%S %Z") else chr(as.character(x)) }
    f$starttime    <- fmt_time(g("starttime"))
    f$endtime      <- fmt_time(g("endtime"))
    f$phredoffset  <- num(g("phredoffset"))

    return(f)
}


#' JAMSalpha_tool_labels
#'
#' Internal. Display names for the dependency keys used in check_resources().
#' @keywords internal
JAMSalpha_tool_labels <- c(sratoolkit = "SRA Toolkit", trimmomatic = "Trimmomatic", seqtk = "seqtk", bowtie2 = "Bowtie2",
                           megahit = "MEGAHIT", spades = "SPAdes", samtools = "SAMtools", metabat2 = "MetaBAT2",
                           checkm2 = "CheckM2", prokka = "Prokka", convert2bed = "BEDOPS convert2bed", bedtools = "BEDTools",
                           kraken2 = "Kraken2", blast = "BLAST+")


#' JAMSalpha_references
#'
#' Internal. Short-form references for the tools and criteria JAMSalpha relies on, keyed as in
#' check_resources(). Please verify against the journal before submission.
#' @keywords internal
JAMSalpha_references <- c(
    JAMS        = "McCulloch JA, et al. JAMS - A framework for the taxonomic and functional exploration of microbiological genomic data. bioRxiv 2023. doi:10.1101/2023.03.03.531026",
    sratoolkit  = "NCBI SRA Toolkit. https://github.com/ncbi/sra-tools",
    trimmomatic = "Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 2014.",
    seqtk       = "Li H. seqtk. https://github.com/lh3/seqtk",
    bowtie2     = "Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods 2012.",
    megahit     = "Li D, et al. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics 2015.",
    spades      = "Bankevich A, et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J Comput Biol 2012.",
    metaspades  = "Nurk S, et al. metaSPAdes: a new versatile metagenomic assembler. Genome Res 2017.",
    rnaspades   = "Bushmanova E, et al. rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. GigaScience 2019.",
    samtools    = "Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience 2021.",
    bedtools    = "Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 2010.",
    convert2bed = "Neph S, et al. BEDOPS: high-performance genomic feature operations. Bioinformatics 2012.",
    kraken2     = "Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol 2019.",
    prokka      = "Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 2014.",
    metabat2    = "Kang DD, et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ 2019.",
    checkm2     = "Chklovski A, et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat Methods 2023.",
    mimag       = "Bowers RM, et al. Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. Nat Biotechnol 2017.",
    blast       = "Camacho C, et al. BLAST+: architecture and applications. BMC Bioinformatics 2009.",
    interproscan = "Jones P, et al. InterProScan 5: genome-scale protein function classification. Bioinformatics 2014."
)


#' rebuild_JAMSalpha_opt_from_jams(list.data = NULL, prefix = NULL)
#'
#' Reconstructs a minimal JAMSalpha opt list for one sample from the objects banked in its .jams
#' file, as returned by load_jamsfiles_from_system(), so that generate_JAMSalpha_methods() and
#' tabulate_JAMSalpha_methods() can describe samples whose workspace image (.RData) is no longer
#' available. Run parameters are recovered from the projinfo table. Anything never banked into a
#' .jams file (e.g. host-contig and denoising records, and the dependency version snapshot for runs
#' predating JAMS v2.3.5) is simply absent, and the sentences that depend on it are omitted. BLAST
#' databases that produced no hits leave no trace in a .jams file and therefore cannot be recovered.
#' If prefix is NULL and list.data holds exactly one sample, that sample is used.
#' @export

rebuild_JAMSalpha_opt_from_jams <- function(list.data = NULL, prefix = NULL){

    if (is.null(list.data) || !is.list(list.data) || length(list.data) == 0){
        stop("list.data must be a non-empty list of JAMS objects, as returned by load_jamsfiles_from_system().")
    }

    projinfo_names <- grep("_projinfo$", names(list.data), value = TRUE)
    if (is.null(prefix)){
        if (length(projinfo_names) == 1){
            prefix <- sub("_projinfo$", "", projinfo_names)
        } else {
            stop(paste("list.data holds", length(projinfo_names), "samples. Specify which one with prefix."))
        }
    }

    #Objects banked by JAMSalpha into .jams files. Anything else of the form <prefix>_<x> is only taken
    #if <x> has no underscore and looks like a BLAST hit table; this stops sample "S1" from swallowing
    #the objects of sample "S1_2" (whose names would leave the suffix "2_projinfo", etc.).
    known_objects <- c("contigsdata", "featuredata", "abundances", "taxonomic_Quality_split_list", "bam_contig_depths",
                       "TNF_contigs", "readstats", "fastqstats", "projinfo", "dependency_versions", "methods", "methods_table")
    stem <- paste0(prefix, "_")
    cand <- names(list.data)[startsWith(names(list.data), stem)]
    suffixes <- substring(cand, nchar(stem) + 1)

    opt <- list()
    blastanalyses <- character(0)
    for (i in seq_along(cand)){
        sfx <- suffixes[i]
        obj <- list.data[[cand[i]]]
        if (sfx %in% known_objects){
            opt[[sfx]] <- obj
        } else if (!grepl("_", sfx) && is.data.frame(obj) && all(c("Feature", "Accession") %in% colnames(obj))){
            opt[[sfx]] <- obj
            blastanalyses <- c(blastanalyses, sfx)
        }
    }

    if (is.null(opt[["projinfo"]])){
        stop(paste("No projinfo object found for sample", prefix, "in list.data."))
    }

    #Parameters banked by gather_run_info(). Multi-element values were unlisted with numeric suffixes
    #appended to their names (e.g. assemblerversion1, assemblerversion2), so fold those back.
    useful_parameters <- c("analysis", "verstr", "JAMS_Kdb_Version", "krakenconfidencescore", "prefix", "threads", "contigs",
                           "assemblyaccession", "readorigin", "seqtype", "libstructure", "R1", "R2", "SE", "sraaccession",
                           "assembler", "assemblerversion", "host", "outdir", "maxbases_input", "useipro")
    pinfo <- as.data.frame(opt[["projinfo"]], stringsAsFactors = FALSE)
    if (all(c("Opt_name", "Run_value") %in% colnames(pinfo))){
        for (rn in seq_len(nrow(pinfo))){
            on <- as.character(pinfo$Opt_name[rn])
            val <- as.character(pinfo$Run_value[rn])
            if (is.na(on)) next
            if (on %in% c("starttime", "endtime")){
                opt[[on]] <- val
                next
            }
            base <- if (on %in% useful_parameters) on else sub("[0-9]+$", "", on)
            if (base %in% useful_parameters){
                opt[[base]] <- c(opt[[base]], val)
            }
        }
        if ("Run_info" %in% colnames(pinfo) && ("Process" %in% pinfo$Run_info)){
            opt[["process"]] <- as.character(pinfo$Run_value[match("Process", pinfo$Run_info)])
        }
    }
    for (nm in c("krakenconfidencescore", "threads", "maxbases_input")){
        if (!is.null(opt[[nm]])) opt[[nm]] <- suppressWarnings(as.numeric(opt[[nm]][1]))
    }
    if (!is.null(opt[["useipro"]])) opt[["useipro"]] <- as.logical(opt[["useipro"]][1])
    if (is.null(opt[["prefix"]])) opt[["prefix"]] <- prefix
    opt[["blastanalyses"]] <- if (length(blastanalyses) > 0) blastanalyses else NULL
    opt[["rebuilt_from_jams"]] <- TRUE

    return(opt)
}


#' generate_JAMSalpha_methods(opt = NULL, list.data = NULL, prefix = NULL,
#'                            include_sample_statistics = TRUE, include_references = TRUE,
#'                            outfile = NULL)
#'
#' Synthesises a prose Methods description of a JAMSalpha run (metagenome, metatranscriptome,
#' isolate or isolate RNA-seq; reads or contigs as input), the JAMSalpha sibling of
#' generate_JAMS16_methods(). Every sentence is derived from what is recorded in opt, and any
#' sentence whose inputs are absent is omitted rather than guessed, so it works equally on a live
#' opt at the end of a run, on an opt restored from an older run's workspace image that never had
#' opt$methods, or on an opt rebuilt from a .jams file. Supply either opt, or list.data (as returned
#' by load_jamsfiles_from_system()) plus prefix, in which case the opt is rebuilt with
#' rebuild_JAMSalpha_opt_from_jams().
#'
#' Tool versions come from opt$dependency_versions (recorded by check_resources()), falling back on
#' opt$assemblerversion for the assembler; the JAMS version is the one recorded in opt$verstr at run
#' time, never the currently installed one.
#' Settings, versions and outcomes come from the run, but the procedures described are those of the
#' installed JAMS code; a warning is logged when the run's JAMS version differs from it.
#'
#' include_sample_statistics = TRUE (default) adds the sample's own figures (read and base counts,
#' host fraction, assembly size and N50, bins obtained, etc.) and names the sample. FALSE gives
#' sample-agnostic text, stating the rules applied rather than the values they produced, which is
#' what you want for a cohort paper. include_references = TRUE appends short-form references for the
#' tools mentioned. Paragraphs are separated by blank lines in $text.
#' @export

generate_JAMSalpha_methods <- function(opt = NULL, list.data = NULL, prefix = NULL,
                                       include_sample_statistics = TRUE, include_references = TRUE,
                                       outfile = NULL){

    if (is.null(opt)){
        if (is.null(list.data)){
            stop("Supply either opt or list.data (with prefix).")
        }
        opt <- rebuild_JAMSalpha_opt_from_jams(list.data = list.data, prefix = prefix)
    }

    f <- describe_JAMSalpha_run(opt = opt)
    stats <- isTRUE(include_sample_statistics)
    if (!f$version_matches_installed){
        flog.warn(paste0("This run was made with ", if (is.na(f$JAMS_version)) "an unrecorded JAMS version" else paste0("JAMS v", f$JAMS_version),
            " but the installed JAMS is v", f$installed_JAMS_version, ". Settings and versions are taken from the run itself, but the procedures ",
            "described are those of the installed code. Check them against the run's version before publishing."))
    }
    #switch() keys must never be NA
    key <- function(x){ if (is.null(x) || length(x) == 0 || is.na(x[1])) "<none>" else as.character(x[1]) }

    # ---- helpers -------------------------------------------------------------
    para <- function(x){ x <- unlist(x); x <- x[!is.na(x) & nzchar(x)]; paste0(x, collapse = " ") }
    fmt_int <- function(x){ format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE) }
    fmt_gbp <- function(x){ paste0(formatC(x / 1e9, format = "f", digits = 2), " Gbp") }
    fmt_mbp <- function(x){ paste0(formatC(x / 1e6, format = "f", digits = 2), " Mbp") }
    fmt_pct <- function(x){ paste0(formatC(x, format = "f", digits = 1), "%") }
    has <- function(x){ !is.null(x) && length(x) > 0 && !all(is.na(x)) }
    #Inline version for a tool, e.g. " v0.39", or "" when unknown
    v <- function(tool){
        vv <- f$versions[tool]
        if (length(vv) == 1 && !is.na(vv)) paste0(" v", sub("^[vV]", "", vv)) else ""
    }
    #Tool name with version in parentheses, e.g. "Trimmomatic (v0.39)" or "Trimmomatic"
    tv <- function(tool, label = JAMSalpha_tool_labels[[tool]]){
        vs <- trimws(v(tool))
        if (nzchar(vs)) paste0(label, " (", vs, ")") else label
    }
    quality_phrase <- function(q){
        if (is.null(q)) return(NA_character_)
        q <- q[q > 0]
        if (length(q) == 0) return(NA_character_)
        labs <- c(HQ = "HQ", MHQ = "MHQ", MQ = "MQ", LQ = "LQ", Contaminated = "contaminated")
        paste0(paste0(q, " ", labs[names(q)]), collapse = ", ")
    }
    #"LKT__s__562_Escherichia_coli" -> "Escherichia coli (LKT__s__562_Escherichia_coli)"
    pretty_lkt <- function(x){
        y <- sub("^LKT__", "", x)
        y <- sub("^[A-Za-z0-9]+__", "", y)
        y <- sub("^[0-9]+_", "", y)
        y <- trimws(gsub("_", " ", y))
        if (nzchar(y) && (y != x)) paste0(y, " (", x, ")") else x
    }
    refs_used <- "JAMS"
    cite <- function(...){ refs_used <<- unique(c(refs_used, c(...))) }

    analysis_label <- switch(key(f$analysis),
        "metagenome"        = "shotgun metagenomic",
        "metatranscriptome" = "metatranscriptomic",
        "isolate"           = "isolate whole-genome",
        "isolaternaseq"     = "isolate RNA-seq",
        "sequencing")
    seq_label <- switch(key(f$seqtype),
        "illuminape" = "Illumina",
        "illuminamp" = "Illumina mate-pair",
        "iontorrent" = "Ion Torrent",
        "pacbio"     = "PacBio",
        NA_character_)
    lib_label <- switch(key(f$libstructure), "pairedend" = "paired-end", "singleend" = "single-end", NA_character_)
    jams_label <- if (!is.na(f$JAMS_version)) paste0("JAMS v", f$JAMS_version) else "JAMS"
    sample_clause <- if (stats && !is.na(f$prefix)) paste0(" for sample ", f$prefix) else ""

    P <- list()   # paragraphs, each a list of sentences added conditionally
    S <- list()   # flat named list of every sentence, for programmatic access

    # ==== 1. Input and read processing =======================================
    R <- list()
    if (identical(f$input_type, "reads")){
        descr <- paste(na.omit(c(lib_label, seq_label)), collapse = " ")
        R$intro <- paste0(trimws(paste(descr, analysis_label)), " reads", sample_clause,
            " were processed with the JAMSalpha pipeline of ", jams_label, " (McCulloch et al. 2023).")
        R$intro <- paste0(toupper(substring(R$intro, 1, 1)), substring(R$intro, 2))

        if (identical(f$readorigin, "sraaccession")){
            cite("sratoolkit")
            R$sra <- paste0("Reads were retrieved from the NCBI Sequence Read Archive",
                if (stats && !is.na(f$sraaccession)) paste0(" (run ", f$sraaccession, ")") else "",
                " with prefetch and fasterq-dump from the ", tv("sratoolkit"), ", skipping technical reads.")
        }

        if (f$seqtype %in% c("illuminape", "illuminamp")){
            cite("trimmomatic")
            phred <- if (!is.na(f$phredoffset)) f$phredoffset else 33
            R$trim <- paste0("Adapter sequences were clipped and reads quality-trimmed with ", tv("trimmomatic"),
                " (ILLUMINACLIP against the JAMS Illumina adapter set with settings 2:30:10, LEADING:15, TRAILING:15, ",
                "SLIDINGWINDOW:4:18, MINLEN:36, Phred+", phred, " quality encoding)",
                if (identical(f$libstructure, "pairedend")) "; reads whose mate did not survive trimming were discarded" else "",
                ".")
        } else if (identical(f$seqtype, "iontorrent")){
            if (identical(f$libstructure, "singleend")){
                cite("spades")
                R$trim <- paste0("Reads were error-corrected with IonHammer as implemented in ", tv("spades"), " (--only-error-correction).")
            } else {
                R$trim <- "Paired-end Ion Torrent reads were not trimmed and were passed directly to SPAdes, run in Ion Torrent mode (--iontorrent)."
            }
        }

        if (f$host_filtered){
            cite("bowtie2")
            #check_resources() renames an NCBI-downloaded host to "<organism>_<GCx accession>"; on an opt rebuilt
            #from a .jams file hostspecies/host_accession_number are gone, so recover them from that name.
            ncbi_host <- regmatches(f$host, regexec("^(.+)_(GC[AF]_[0-9]+\\.?[0-9]*)$", f$host))[[1]]
            host_acc  <- if (!is.na(f$host_accession_number)) f$host_accession_number else if (length(ncbi_host) == 3) ncbi_host[3] else NA_character_
            host_org  <- if (!is.na(f$hostspecies)) f$hostspecies else if (length(ncbi_host) == 3) ncbi_host[2] else f$host
            host_label <- switch(f$host,
                "human" = "human (Homo sapiens)",
                "mouse" = "mouse (Mus musculus)",
                if (!is.na(host_acc)) paste0(host_org, " (NCBI assembly ", host_acc, ", downloaded and indexed with bowtie2-build)")
                else paste0("user-supplied host (", f$host, ")"))
            R$host <- paste0("To remove host-derived sequence, trimmed reads were aligned to the ", host_label,
                " reference genome", if (f$host %in% c("human", "mouse")) " supplied with the JAMS database" else "",
                " with ", tv("bowtie2"), " (--very-fast-local --ignore-quals), and only ",
                if (identical(f$libstructure, "pairedend")) "read pairs that did not align concordantly" else "reads that did not align", " were retained",
                if (identical(f$analysis, "isolaternaseq")) " (for isolate RNA-seq, however, the trimmed rather than the host-depleted reads were used for assembly and abundance estimation)" else "",
                ".")
        }

        if (f$subsampled){
            cite("seqtk")
            R$subsample <- paste0(if (stats) "As the number of bases available for assembly exceeded the cap of " else "Where the number of bases available for assembly exceeded the cap of ",
                if (!is.na(f$maxbases_input)) paste0(fmt_int(f$maxbases_input), " bases") else "the user-defined maximum",
                ", reads were randomly subsampled to approximately that amount with ", tv("seqtk"),
                " (sample, seed 100), keeping read pairs together.")
        } else if (!stats && !is.na(f$maxbases_input)){
            R$subsample <- paste0("Where the number of bases available for assembly exceeded ", fmt_int(f$maxbases_input),
                ", reads were randomly subsampled to approximately that amount with ", tv("seqtk"), " (sample, seed 100).")
            cite("seqtk")
        }

        if (stats){
            bits <- character(0)
            if (!is.na(f$raw_reads)){
                bits <- c(bits, paste0("The sample yielded ", fmt_int(f$raw_reads), " raw reads",
                    if (!is.na(f$raw_pairs)) paste0(" (", fmt_int(f$raw_pairs), " pairs)") else "",
                    if (!is.na(f$raw_readlength)) paste0(" with a mean length of ", f$raw_readlength, " bp") else "",
                    if (!is.na(f$bases_raw)) paste0(", totalling ", fmt_gbp(f$bases_raw)) else "", "."))
            } else if (!is.na(f$bases_raw)){
                bits <- c(bits, paste0("The sample yielded ", fmt_gbp(f$bases_raw), " of raw sequence."))
            }
            if (!is.na(f$pct_bases_trimmed_kept)){
                bits <- c(bits, paste0(fmt_pct(f$pct_bases_trimmed_kept), " of bases (", fmt_gbp(f$bases_trimmed), ") survived quality trimming."))
            }
            if (!is.na(f$pct_host)){
                bits <- c(bits, paste0(fmt_pct(f$pct_host), " of trimmed ", if (identical(f$pct_host_basis, "reads")) "reads" else "bases",
                    " aligned to the host genome",
                    if (!is.na(f$bases_nonhost)) paste0(", leaving ", fmt_gbp(f$bases_nonhost), " of non-host sequence") else "", "."))
            }
            if (f$subsampled && !is.na(f$bases_subsampled)){
                bits <- c(bits, paste0("After subsampling, approximately ", fmt_gbp(f$bases_subsampled), " were used for assembly."))
            }
            if (length(bits) > 0) R$readstats <- paste(bits, collapse = " ")
        }
    } else if (f$input_type %in% c("contigs", "assemblyaccession")){
        src <- if (identical(f$input_type, "assemblyaccession")){
            paste0("retrieved from NCBI GenBank", if (stats) paste0(" (assembly accession ", f$assemblyaccession, ")") else "")
        } else {
            "supplied as a FASTA file"
        }
        R$intro <- paste0("Pre-assembled ", if (analysis_label != "sequencing") paste0(analysis_label, " ") else "",
            "contigs", sample_clause, ", ", src, ", were processed with the JAMSalpha pipeline of ", jams_label, " (McCulloch et al. 2023).")
        R$intro <- paste0(toupper(substring(R$intro, 1, 1)), substring(R$intro, 2))
    }
    P$reads <- R

    # ==== 2. Assembly ========================================================
    A <- list()
    if (identical(f$input_type, "reads") && !is.na(f$assembler)){
        if (identical(f$assembler, "megahit")){
            cite("megahit")
            mh <- f$megahit
            if (stats && !is.null(mh) && !is.na(mh$kmin) && !is.na(mh$kmax)){
                A$assembler <- paste0("Reads were assembled with ", tv("megahit"), " using k-mer sizes from ", mh$kmin, " to ", mh$kmax,
                    if (!is.na(mh$kstep)) paste0(" in steps of ", mh$kstep) else "",
                    if (!is.na(mh$mincount)) paste0(" and a minimum k-mer multiplicity of ", mh$mincount) else "",
                    #estcoverage of 0 is the assemble_contigs() artefact of -b being set but not exceeded; do not print it.
                    if (!is.na(f$estcoverage) && f$estcoverage > 0) paste0(", parameters chosen by JAMS from the estimated sequencing coverage (", f$estcoverage,
                        "x, assuming 50 genomes of 4 Mbp)") else "",
                    ".")
            } else {
                A$assembler <- paste0("Reads were assembled with ", tv("megahit"), ", with k-mer parameters chosen by JAMS from the estimated ",
                    "sequencing coverage (total input bases divided by 50 genomes of 4 Mbp): ",
                    "minimum k-mer size 21, 23 or 27, k-mer step 8, 10 or 12 and minimum k-mer multiplicity 1, 2 or 3 for estimated coverages of ",
                    "at most 5x, 5-40x and above 40x, respectively, with the maximum k-mer size set to the largest odd number not exceeding the read length (at most 255).")
            }
        } else if (identical(f$assembler, "spades")){
            cite("spades")
            mode <- switch(key(f$analysis),
                "metagenome"        = { cite("metaspades"); " in metagenomic mode (metaSPAdes, --meta)" },
                "metatranscriptome" = { cite("rnaspades");  " in transcriptomic mode (rnaSPAdes, --rna)" },
                "isolaternaseq"     = { cite("rnaspades");  " in transcriptomic mode (rnaSPAdes, --rna)" },
                "isolate"           = " in isolate mode (--isolate --cov-cutoff auto)",
                "")
            extra <- character(0)
            if (identical(f$seqtype, "iontorrent")) extra <- c(extra, "Ion Torrent mode (--iontorrent)")
            if (identical(f$seqtype, "illuminamp")) extra <- c(extra, "reads supplied as a high-quality mate-pair library")
            A$assembler <- paste0("Reads were assembled with ", tv("spades"), mode,
                if (length(extra) > 0) paste0(", with ", paste(extra, collapse = " and ")) else "", ".")
        }
    }
    if (!is.na(f$input_type)){
        A$minlength <- paste0("Contigs shorter than ", f$contig_minlength, " bp were discarded and the remainder renamed consecutively.")
    }
    if (stats && !is.na(f$n_contigs)){
        A$stats <- paste0("The final contig set comprised ", fmt_int(f$n_contigs), " contigs",
            " (total length ", fmt_mbp(f$contigs_total_length), "; N50 ", fmt_int(f$contigs_N50), " bp; longest contig ", fmt_int(f$contigs_max), " bp)",
            if (identical(f$input_type, "reads") && !is.na(f$pct_assembled)) paste0(", and ", fmt_pct(f$pct_assembled), " of the bases used for assembly aligned back to them") else "",
            ".")
    }
    P$assembly <- A

    # ==== 3. Contig classification, annotation and abundance =================
    TX <- list()
    if (!is.na(f$input_type)){
        cite("kraken2")
        TX$kraken <- paste0("Each contig was taxonomically classified with ", tv("kraken2"), " against the JAMS Kraken2 database",
            if (!is.na(f$kdb_version)) paste0(" (version ", f$kdb_version, ")") else "",
            if (!is.na(f$kraken_confidence)) paste0(" with a confidence threshold of ", format(f$kraken_confidence, scientific = FALSE)) else "",
            ", and its NCBI taxid was expanded into a full lineage and a Last Known Taxon (LKT), the most specific taxonomic rank resolved, using the JAMS taxonomy table.")
        TX$metazoa <- paste0("Contigs classified as Metazoa (e.g. residual host sequence) were removed",
            if (stats && has(f$n_metazoan_contigs)) paste0(" (", fmt_int(f$n_metazoan_contigs), " contig", if (f$n_metazoan_contigs != 1) "s" else "", ")") else "",
            ".")
    }
    #Annotation precedes read mapping in the pipeline (annotate_contigs before get_contig_coverage)
    if (!is.na(f$input_type)){
        cite("prokka", "convert2bed")
        TX$prokka <- paste0("Genes were predicted and annotated on the retained contigs with ", tv("prokka"),
            " (--addgenes, genetic code 11", if (identical(f$analysis, "metagenome")) ", --metagenome" else "", "; tbl2asn step skipped), ",
            "and the resulting GFF file was converted to BED format with convert2bed (BEDOPS", v("convert2bed"), ").")
        if (stats && !is.na(f$n_CDS)){
            TX$prokkastats <- paste0("This yielded ", fmt_int(f$n_CDS), " protein-coding sequences, ", fmt_int(f$n_rRNA), " rRNA and ", fmt_int(f$n_tRNA), " tRNA genes.")
        }
    }
    if (identical(f$input_type, "reads")){
        cite("bowtie2", "samtools", "bedtools")
        TX$coverage <- paste0("To quantify abundance, the reads used for assembly were aligned back to the contigs as unpaired reads with ", tv("bowtie2"),
            "; alignments with a mapping quality below 10 were discarded with ", tv("samtools"),
            ", and per-base sequencing depth over contigs and genes was computed with ", tv("bedtools"), " (genomecov and coverage). ",
            "The abundance of each contig and gene was expressed as the number of sequenced bases aligned to it (the sum of its per-base depths), ",
            "and relative abundance as parts per million (PPM) of all bases aligned to the sample's contigs.")
    } else if (f$input_type %in% c("contigs", "assemblyaccession")){
        TX$coverage <- "As no reads were available, the abundance of each contig and gene was set to its length, i.e. a uniform depth of 1x."
    }
    if (f$denoise_applied){
        TX$denoise <- paste0("To limit spurious taxonomic assignments, the classification of contigs shorter than 2,000 bp was reset to unclassified when their species ",
            "was not represented among contigs of at least 2,000 bp and the short contigs of that species together amounted to less than 10% of the median genome size ",
            "recorded for its taxid or to less than 50 PPM of the sample's sequencing depth",
            if (stats && has(f$n_censored_contigs)) paste0(" (", fmt_int(f$n_censored_contigs), " contig", if (f$n_censored_contigs != 1) "s" else "", " affected)") else "",
            ".")
    }
    P$contigs <- TX

    # ==== 4. Genome completeness, binning and consolidation ==================
    B <- list()
    checkm_db_clause <- if (!is.na(f$checkm_db)) paste0(" with the ", f$checkm_db, " database") else ""
    if (identical(f$analysis, "metagenome")){
        cite("checkm2")
        B$lkt <- paste0("The contigs assigned to each LKT were then evaluated as a putative genome: completeness and contamination were estimated with ",
            tv("checkm2"), checkm_db_clause, " (predict) for bacterial and archaeal LKTs whose contigs summed to at least 5 kb, and otherwise as the summed contig length ",
            "relative to the median genome size of the taxon, with completeness above 100% capped at 100% and the excess added to contamination.")
    }
    if (!identical(f$analysis, "isolate") && (f$mags_present || f$cgb_present || identical(f$analysis, "metagenome"))){
        cite("mimag")
        B$quality <- paste0("Genomes and bins were rated with modified MIMAG criteria (Bowers et al. 2017): high quality (HQ), completeness above 90% and contamination below 5%; ",
            "medium-high quality (MHQ), above 70% and below 10%; medium quality (MQ), above 50% and below 10%; low quality (LQ), all others; ",
            "bins with contamination above 20% were classed as contaminated.")
    }
    if (f$binning_applied){
        cite("metabat2")
        B$mags <- paste0("Contigs were binned into metagenome-assembled genomes (MAGs) with ", tv("metabat2"),
            " (default settings, including its 2,500 bp minimum contig length) using contig depths computed with jgi_summarize_bam_contig_depths from the read alignments above. ",
            "Each bin was assessed with CheckM2 and classified as a single entity by running Kraken2 (confidence 0) on its concatenated sequence; ",
            "bins were named after their bin number and LKT.")
        if (stats && has(f$n_mags)){
            qp <- quality_phrase(f$mag_quality)
            B$magstats <- paste0("MetaBAT2 produced ", fmt_int(f$n_mags), " bin", if (f$n_mags != 1) "s" else "",
                if (!is.na(qp)) paste0(" (", qp, ")") else "", ".")
        }
    }
    if (f$cgb_present && !identical(f$analysis, "isolate")){
        B$cgb <- paste0("Contigs were then consolidated into a single non-redundant set of genome bins per sample (consolidated genome bins, CGBs). ",
            "HQ and MHQ MetaBAT2 bins were retained first; where such a bin was resolved to species or strain level and was below 99% complete (considering one bin per species), ",
            "contigs not already in an HQ or MHQ bin that Kraken2 assigned to the same species were added to it if the augmented bin still rated HQ or MHQ on re-evaluation with CheckM2. ",
            "The remaining contigs were grouped by LKT, and LKT completeness and contamination estimates were summed within each taxon at successively higher taxonomic levels ",
            "(infraspecies, species, genus, family, order, class, phylum and kingdom; summed completeness above 100% counted as contamination). ",
            "For each quality tier in turn (HQ, MHQ, MQ, LQ), each LKT's contigs were consolidated at the most specific of these levels at which the aggregate reached that tier ",
            "(or at the LKT itself if no higher level did), and every candidate bin was then verified with CheckM2 for bacteria and archaea (bins of at least 5 kb) ",
            "or estimated from contig length relative to the expected genome size for other taxa, being accepted only if it met the tier. ",
            "A final round consolidated the remaining contigs by LKT irrespective of quality, and any contigs left after that were also consolidated by their LKT.")
        if (stats && has(f$n_cgb)){
            qp <- quality_phrase(f$cgb_quality)
            B$cgbstats <- paste0("This yielded ", fmt_int(f$n_cgb), " CGBs", if (!is.na(qp)) paste0(" (", qp, ")") else "",
                if (has(f$n_cgb_from_mb2) && f$n_cgb_from_mb2 == 0) ", none of which derived from MetaBAT2 bins" else
                if (has(f$n_cgb_from_mb2)) paste0(", ", fmt_int(f$n_cgb_from_mb2), " of which derived from MetaBAT2 bins",
                    if (has(f$pct_depth_mb2)) paste0(" and together accounted for ", fmt_pct(f$pct_depth_mb2), " of the sample's sequencing depth") else "") else "",
                ".")
        }
    }
    if (identical(f$analysis, "isolate") && (f$cgb_present || !is.na(f$input_type))){
        cite("checkm2", "mimag")
        B$isolate <- paste0("The isolate genome was evaluated as a single entity: completeness and contamination were estimated with ", tv("checkm2"), checkm_db_clause,
            " (predict) on all contigs classified as bacterial or archaeal, the concatenated sequence was classified with Kraken2 (confidence 0) to assign a single LKT to the genome, ",
            "and genome quality was rated with modified MIMAG criteria (Bowers et al. 2017; HQ, completeness above 90% and contamination below 5%; MHQ, above 70% and below 10%; ",
            "MQ, above 50% and below 10%; LQ, all others; contaminated, contamination above 20%).")
        if (stats && !is.null(f$isolate_genome)){
            ig <- f$isolate_genome
            bits <- character(0)
            if (!is.na(ig$LKT)) bits <- c(bits, paste0("The genome was classified as ", pretty_lkt(ig$LKT)))
            cc <- character(0)
            if (!is.na(ig$Completeness))  cc <- c(cc, paste0(formatC(ig$Completeness, format = "f", digits = 2), "% completeness"))
            if (!is.na(ig$Contamination)) cc <- c(cc, paste0(formatC(ig$Contamination, format = "f", digits = 2), "% contamination"))
            if (length(cc) > 0) bits <- c(bits, paste0(if (length(bits) > 0) ", with " else "The genome had ", paste(cc, collapse = " and ")))
            if (!is.na(ig$Quality) && length(bits) > 0) bits <- c(bits, paste0(" (", ig$Quality, ")"))
            if (length(bits) > 0) B$isolatestats <- paste0(paste0(bits, collapse = ""), ".")
            if (has(f$n_nonbacterial_contigs) && f$n_nonbacterial_contigs > 0){
                B$isolatestats <- trimws(paste0(if (is.null(B$isolatestats)) "" else B$isolatestats, " ", fmt_int(f$n_nonbacterial_contigs), " contig", if (f$n_nonbacterial_contigs != 1) "s were" else " was",
                    " classified as neither bacterial nor archaeal."))
            }
        }
    }
    P$binning <- B

    # ==== 5. Functional annotation and profiling ============================
    Fn <- list()
    if (length(f$blastanalyses) > 0){
        cite("blast")
        dbnames <- c(resfinder = "ResFinder", plasmidfinder = "PlasmidFinder", vfdb = "VFDB", abricate = "ABRicate", napdos = "NaPDoS",
                     serofinderH = "SerotypeFinder (H antigen)", serofinderO = "SerotypeFinder (O antigen)")
        dbl <- ifelse(f$blastanalyses %in% names(dbnames), dbnames[f$blastanalyses], f$blastanalyses)
        dbl <- unname(dbl)
        dbtxt <- if (length(dbl) > 1) paste0(paste(dbl[-length(dbl)], collapse = ", "), " and ", dbl[length(dbl)]) else dbl
        Fn$blast <- paste0("Predicted genes were searched against the ", dbtxt, " database", if (length(dbl) > 1) "s" else "",
            " bundled with the JAMS database using ", tv("blast"), " (blastn for nucleotide and blastp, with the predicted proteins, for protein databases; ",
            "e-value at most 1e-20, culling limit 1), retaining hits with query coverage, subject coverage and identity each above 75%.")
        if (stats){
            hb <- f$blast_hits[f$blast_hits > 0]
            if (length(hb) > 0){
                hl <- ifelse(names(hb) %in% names(dbnames), dbnames[names(hb)], names(hb))
                Fn$blaststats <- paste0("Genes with hits: ", paste0(unname(hl), " ", hb, collapse = "; "), ".")
            }
        }
    }
    if (f$ipro_ran){
        cite("interproscan")
        Fn$ipro <- paste0("Predicted proteins, excluding those with internal stop codons, were further annotated with InterProScan (with GO term and pathway lookups)",
            if (stats) paste0(", yielding ", paste(f$ipro_analyses, collapse = ", "), " signatures") else "", ".")
    }
    if (length(f$func_spaces) > 0){
        sig <- character(0)
        if ("Product" %in% f$func_basic_analyses)  sig <- c(sig, "Prokka product descriptions")
        if ("ECNumber" %in% f$func_basic_analyses) sig <- c(sig, "Enzyme Commission (EC) numbers")
        if ("FeatType" %in% f$func_basic_analyses) sig <- c(sig, "feature types")
        if (length(f$blastanalyses) > 0)           sig <- c(sig, "BLAST database hits")
        if (f$ipro_ran)                            sig <- c(sig, "InterProScan signatures, InterPro entries and GO terms")
        space_labels <- c(Contig_LKT = "contig LKT", MB2bin = "MetaBAT2 bin", ConsolidatedGenomeBin = "consolidated genome bin")
        sl <- unname(ifelse(f$func_spaces %in% names(space_labels), space_labels[f$func_spaces], f$func_spaces))
        sigtxt <- if (length(sig) > 1) paste0(paste(sig[-length(sig)], collapse = ", "), " and ", sig[length(sig)]) else sig
        Fn$functional <- paste0("Functional abundance tables were built by summing the abundance of all genes carrying each functional signature",
            if (length(sig) > 0) paste0(" (", sigtxt, ")") else "",
            ", stratified by the taxonomic unit to which each gene's contig was assigned (", paste(sl, collapse = ", "), ").")
    }
    P$functional <- Fn

    # ==== 6. Software and databases ==========================================
    W <- list()
    vk <- f$versions[!is.na(f$versions)]
    if (!is.na(f$JAMS_version) || length(vk) > 0){
        parts <- character(0)
        if (!is.na(f$JAMS_version)) parts <- c(parts, paste0("JAMS v", f$JAMS_version, if (!is.na(f$R_version)) paste0(" under ", f$R_version) else ""))
        if (length(vk) > 0) parts <- c(parts, paste0(JAMSalpha_tool_labels[names(vk)], " v", sub("^[vV]", "", vk)))
        dbs <- character(0)
        if (!is.na(f$kdb_version)) dbs <- c(dbs, paste0("JAMS Kraken2 database version ", f$kdb_version))
        if (!is.na(f$checkm_db) && ("checkm2" %in% f$tools_used)) dbs <- c(dbs, paste0("CheckM2 database ", f$checkm_db))
        W$versions <- paste0("Software versions: ", paste(parts, collapse = ", "), ".",
            if (length(dbs) > 0) paste0(" Reference databases: ", paste(dbs, collapse = "; "), ".") else "")
    }
    P$software <- W

    # ---- assemble ------------------------------------------------------------
    P <- P[vapply(P, function(x){ length(x) > 0 }, logical(1))]
    paragraphs <- vapply(P, para, character(1))
    paragraphs <- paragraphs[nzchar(paragraphs)]
    methods_text <- paste(paragraphs, collapse = "\n\n")

    for (pn in names(P)){
        for (sn in names(P[[pn]])){
            S[[paste(pn, sn, sep = "_")]] <- P[[pn]][[sn]]
        }
    }

    refs <- unname(JAMSalpha_references[intersect(names(JAMSalpha_references), refs_used)])
    if (isTRUE(include_references) && length(refs) > 0){
        methods_text <- paste0(methods_text, "\n\nReferences\n", paste(refs, collapse = "\n"))
    }

    if (!is.null(outfile)){
        writeLines(methods_text, con = outfile)
        flog.info(paste("Wrote methods narrative to", outfile))
    }

    return(invisible(list(text = methods_text, paragraphs = paragraphs, sentences = S, versions = f$versions,
                          references = refs, rebuilt_from_jams = f$rebuilt_from_jams)))
}


#' tabulate_JAMSalpha_methods(opt = NULL, list.data = NULL, prefix = NULL, outfile = NULL)
#'
#' Builds a tidy Category / Parameter / Value data frame of every methods-relevant choice and
#' outcome of a JAMSalpha run, for machine-readable provenance alongside the prose of
#' generate_JAMSalpha_methods() (the JAMSalpha sibling of tabulate_JAMS16_methods()). Rows whose
#' underlying value is absent are omitted, so the table never carries empty or misleading cells.
#' Rows whose value was re-derived rather than recorded say so in a companion *_source row.
#' Supply either opt, or list.data plus prefix. Optionally writes the table to TSV.
#' @export

tabulate_JAMSalpha_methods <- function(opt = NULL, list.data = NULL, prefix = NULL, outfile = NULL){

    if (is.null(opt)){
        if (is.null(list.data)){
            stop("Supply either opt or list.data (with prefix).")
        }
        opt <- rebuild_JAMSalpha_opt_from_jams(list.data = list.data, prefix = prefix)
    }

    f <- describe_JAMSalpha_run(opt = opt)

    rows <- list()
    #Small accumulator: only lands a row when the value is non-NULL / non-NA / non-empty.
    add <- function(parameter, value, category = NA_character_){
        if (is.null(value)) return(invisible(NULL))
        if (is.numeric(value)){
            value <- ifelse(is.na(value), NA_character_, format(value, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
        }
        value <- as.character(value)
        if (length(value) == 0 || all(is.na(value)) || all(value == "", na.rm = TRUE)) return(invisible(NULL))
        rows[[length(rows) + 1]] <<- data.frame(
            Category  = category,
            Parameter = parameter,
            Value     = paste0(value[!is.na(value)], collapse = ", "),
            stringsAsFactors = FALSE)
    }
    yn <- function(x){ if (is.na(x)) NA_character_ else if (isTRUE(x)) "yes" else "no" }
    rnd <- function(x, d = 2){ if (is.null(x) || is.na(x)) NA else round(x, d) }

    # ---- Project / input -----------------------------------------------------
    add("Sample_name",        f$prefix,            "Project")
    add("Analysis_type",      f$analysis,          "Project")
    add("Input_type",         f$input_type,        "Input")
    add("Reads_source",       f$readorigin,        "Input")
    add("SRA_accession",      f$sraaccession,      "Input")
    add("GenBank_assembly_accession", f$assemblyaccession, "Input")
    add("Sequencing_platform", f$seqtype,          "Input")
    add("Library_structure",  f$libstructure,      "Input")

    # ---- Reads ---------------------------------------------------------------
    if (identical(f$input_type, "reads")){
        add("Raw_reads",              f$raw_reads,       "Read_depth")
        add("Raw_read_pairs",         f$raw_pairs,       "Read_depth")
        add("Median_raw_read_length_bp", f$raw_readlength, "Read_depth")
        add("Raw_bases",              f$bases_raw,       "Read_depth")
        add("Trimmed_bases",          f$bases_trimmed,   "Read_depth")
        add("Pct_bases_surviving_trimming", rnd(f$pct_bases_trimmed_kept), "Read_depth")
        if (f$seqtype %in% c("illuminape", "illuminamp")){
            #Trimmomatic constants, hard-wired in trim_reads(); recording them makes the row self-contained.
            add("Trimming_tool",       "Trimmomatic",   "Read_QC")
            add("ILLUMINACLIP",        "JAMS seqadapters:2:30:10", "Read_QC")
            add("LEADING",             "15",            "Read_QC")
            add("TRAILING",            "15",            "Read_QC")
            add("SLIDINGWINDOW",       "4:18",          "Read_QC")
            add("MINLEN",              "36",            "Read_QC")
            add("Phred_offset",        if (!is.na(f$phredoffset)) f$phredoffset else 33, "Read_QC")
            add("Unpaired_trimmed_reads_used", "no",    "Read_QC")
        } else if (identical(f$seqtype, "iontorrent")){
            add("Read_correction", if (identical(f$libstructure, "singleend")) "IonHammer (spades.py --only-error-correction)" else "none (SPAdes --iontorrent)", "Read_QC")
        }
        add("Host_species",           f$host,            "Host")
        if (f$host_filtered){
            add("Host_genome_label",      f$hostspecies,           "Host")
            add("Host_genome_accession",  f$host_accession_number, "Host")
            add("Host_aligner_settings",  "bowtie2 --very-fast-local --ignore-quals", "Host")
            add("Pct_host",               rnd(f$pct_host),         "Host")
            add("Pct_host_basis",         f$pct_host_basis,        "Host")
            add("Non_host_bases",         f$bases_nonhost,         "Host")
        }
        add("Assembly_input_reads",   if (f$assembly_input == "NAHS") "host-depleted (NAHS)" else "trimmed", "Assembly")
        add("Max_bases_for_assembly", f$maxbases_input,  "Subsampling")
        add("Subsampling_applied",    yn(f$subsampled),  "Subsampling")
        if (f$subsampled){
            add("Subsampling_tool",   "seqtk sample -s100", "Subsampling")
            add("Subsampled_bases",   f$bases_subsampled,   "Subsampling")
        }
    }

    # ---- Assembly ------------------------------------------------------------
    add("Assembler",                  f$assembler,        "Assembly")
    add("Assembler_command",          f$assembler_command, "Assembly")
    add("Estimated_coverage_x",       f$estcoverage,      "Assembly")
    if (!is.null(f$deepseq)) add("Deep_sequencing_flag", yn(f$deepseq), "Assembly")
    if (!is.null(f$megahit)){
        add("MEGAHIT_k_min",          f$megahit$kmin,     "Assembly")
        add("MEGAHIT_k_max",          f$megahit$kmax,     "Assembly")
        add("MEGAHIT_k_step",         f$megahit$kstep,    "Assembly")
        add("MEGAHIT_min_count",      f$megahit$mincount, "Assembly")
        add("MEGAHIT_parameters_source", f$megahit$source, "Assembly")
    }
    if (identical(f$assembler, "spades")){
        add("SPAdes_mode", switch(if (is.na(f$analysis)) "<none>" else f$analysis, "metagenome" = "--meta", "metatranscriptome" = "--rna", "isolaternaseq" = "--rna",
                                  "isolate" = "--isolate --cov-cutoff auto", "default"), "Assembly")
    }
    if (!is.na(f$input_type)) add("Min_contig_length_bp", f$contig_minlength, "Assembly")
    add("Num_contigs",                f$n_contigs,             "Assembly")
    add("Total_contig_length_bp",     f$contigs_total_length,  "Assembly")
    add("Contig_N50_bp",              f$contigs_N50,           "Assembly")
    add("Max_contig_length_bp",       f$contigs_max,           "Assembly")
    add("Assembled_bases",            f$bases_assembled,       "Assembly")
    add("Pct_input_bases_assembled",  rnd(f$pct_assembled),    "Assembly")

    # ---- Taxonomy ------------------------------------------------------------
    if (!is.na(f$input_type)){
        add("Contig_classifier",      "kraken2",               "Taxonomy")
        if (!is.na(f$kraken_confidence)) add("Kraken2_confidence", format(f$kraken_confidence, scientific = FALSE), "Taxonomy")
        add("JAMS_Kraken2_db_version", f$kdb_version,          "Taxonomy")
        add("Metazoan_contigs_removed", f$n_metazoan_contigs,  "Taxonomy")
    }
    add("Num_contig_LKTs",            f$n_contig_LKTs,         "Taxonomy")

    # ---- Abundance -----------------------------------------------------------
    if (identical(f$input_type, "reads")){
        add("Abundance_measure",      "NumBases (sum of per-base depth); relative abundance in PPM", "Abundance")
        add("Coverage_aligner_settings", "bowtie2 (default end-to-end, --no-unal), reads aligned as unpaired", "Abundance")
        add("Min_mapping_quality",    "10",                    "Abundance")
    } else if (!is.na(f$input_type)){
        add("Abundance_measure",      "contig/gene length (uniform 1x depth; no reads)", "Abundance")
    }

    # ---- Denoising -----------------------------------------------------------
    if (f$denoise_applied){
        add("Denoising_applied",               "yes",   "Denoising")
        add("Denoising_short_contig_threshold_bp", "2000", "Denoising")
        add("Denoising_min_est_completeness_pct", "10", "Denoising")
        add("Denoising_min_PPM",               "50",    "Denoising")
        add("Contigs_taxonomy_censored",       f$n_censored_contigs, "Denoising")
    }

    # ---- Completeness / binning / consolidation ------------------------------
    if (f$analysis %in% c("metagenome", "isolate")){
        add("Completeness_tool",      "CheckM2 predict",       "Completeness")
        add("CheckM2_database",       f$checkm_db,             "Completeness")
        add("CheckM2_min_bin_length_bp", "5000",               "Completeness")
        if (identical(f$analysis, "metagenome")) add("Completeness_fallback", "contig length sum / median taxid genome size (excess over 100% added to contamination)", "Completeness")
        add("Quality_HQ",             "completeness > 90 and contamination < 5",  "Completeness")
        add("Quality_MHQ",            "completeness > 70 and contamination < 10", "Completeness")
        add("Quality_MQ",             "completeness > 50 and contamination < 10", "Completeness")
        add("Quality_LQ",             "all other bins",                           "Completeness")
        add("Quality_Contaminated",   "contamination > 20",                       "Completeness")
    }
    if (f$binning_applied){
        add("Binning_tool",           "MetaBAT2 (default settings)", "Binning")
        add("Num_MAGs",               f$n_mags,                "Binning")
        if (!is.null(f$mag_quality)) for (q in names(f$mag_quality)) add(paste0("MAGs_", q), f$mag_quality[[q]], "Binning")
    }
    if (f$cgb_present && !identical(f$analysis, "isolate")){
        add("Num_CGBs",               f$n_cgb,                 "Consolidation")
        if (!is.null(f$cgb_quality)) for (q in names(f$cgb_quality)) add(paste0("CGBs_", q), f$cgb_quality[[q]], "Consolidation")
        add("CGBs_from_MetaBAT2",     f$n_cgb_from_mb2,        "Consolidation")
        add("Pct_depth_in_MetaBAT2_CGBs", rnd(f$pct_depth_mb2), "Consolidation")
    }
    if (!is.null(f$isolate_genome)){
        add("Isolate_LKT",            f$isolate_genome$LKT,           "Isolate")
        add("Isolate_completeness",   f$isolate_genome$Completeness,  "Isolate")
        add("Isolate_contamination",  f$isolate_genome$Contamination, "Isolate")
        add("Isolate_quality",        f$isolate_genome$Quality,       "Isolate")
        add("Isolate_genome_size_bp", f$isolate_genome$Genome_Size,   "Isolate")
        add("Isolate_nonbacterial_contigs", f$n_nonbacterial_contigs, "Isolate")
    }

    # ---- Annotation / functional ---------------------------------------------
    if (!is.na(f$input_type)){
        add("Annotation_tool",        paste0("prokka --addgenes --gcode 11", if (identical(f$analysis, "metagenome")) " --metagenome" else "", " (tbl2asn skipped)"), "Annotation")
    }
    add("Num_CDS",                    f$n_CDS,                 "Annotation")
    add("Num_rRNA",                   f$n_rRNA,                "Annotation")
    add("Num_tRNA",                   f$n_tRNA,                "Annotation")
    if (length(f$blastanalyses) > 0){
        add("BLAST_databases",        f$blastanalyses,         "BLAST")
        add("BLAST_settings",         "evalue 1e-20, culling_limit 1, max_target_seqs 10000", "BLAST")
        add("BLAST_hit_filters",      "query coverage > 75, subject coverage > 75, identity > 75", "BLAST")
        for (b in names(f$blast_hits)) add(paste0("Genes_with_hits_", b), f$blast_hits[[b]], "BLAST")
    }
    add("InterProScan_requested",     yn(f$useipro),           "InterProScan")
    if (f$ipro_ran) add("InterProScan_analyses", f$ipro_analyses, "InterProScan")
    add("Functional_analyses",        f$func_analyses,         "Functional")
    add("Taxonomic_spaces",           f$tax_spaces,            "Functional")
    add("Functional_stratification",  f$func_spaces,           "Functional")

    # ---- Software ------------------------------------------------------------
    add("Version_JAMS",               f$JAMS_version,          "Software")
    add("Procedures_described_per_JAMS_version", f$installed_JAMS_version, "Software")
    add("Run_version_matches_installed", yn(f$version_matches_installed), "Software")
    add("Version_R",                  f$R_version,             "Software")
    for (tool in names(f$versions)){
        add(paste0("Version_", tool), f$versions[[tool]],      "Software")
        if (!is.na(f$version_source[[tool]]) && f$version_source[[tool]] != "dependency_versions"){
            add(paste0("Version_", tool, "_source"), f$version_source[[tool]], "Software")
        }
    }
    missing_v <- names(f$versions)[is.na(f$versions)]
    if (length(missing_v) > 0) add("Versions_not_recorded", missing_v, "Software")

    # ---- Provenance ----------------------------------------------------------
    add("Threads",                    f$threads,               "Provenance")
    add("Run_start",                  f$starttime,             "Provenance")
    add("Run_end",                    f$endtime,               "Provenance")
    add("Description_source", if (f$rebuilt_from_jams) "rebuilt_from_jams_file" else "opt", "Provenance")
    add("Table_generated",            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "Provenance")

    # ---- assemble ------------------------------------------------------------
    methods_table <- if (length(rows) > 0) do.call(rbind, rows) else
        data.frame(Category = character(0), Parameter = character(0), Value = character(0))
    rownames(methods_table) <- NULL

    if (!is.null(outfile)){
        write.table(methods_table, file = outfile, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
        flog.info(paste("Wrote methods parameter table to", outfile))
    }

    return(invisible(methods_table))
}
