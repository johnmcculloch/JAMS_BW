#' check_for_banked_samples(search_path = NULL, query_samples = NULL, mode = c("inventory", "query"), recursive = TRUE, flag_duplicates = TRUE)
#'
#' Scans a directory tree for banked JAMS samples and reports what is available.
#'
#' Unlike the previous version, which required the EXACT path to a single JAMSarchive
#' (with jamsfiles/ and JAMStarballs/ immediately inside), this version searches
#' recursively from any upstream path. It can therefore span multiple JAMSarchive
#' directories belonging to different projects, even when these sit at different
#' depths. A "banked sample" is defined as a prefix for which BOTH a .jams file and a
#' JAMSalpha run tarball (_JAMS.tar.gz) exist within the same JAMSarchive.
#'
#' Two modes of use:
#'   - "inventory" (default): returns a data frame of every banked sample found
#'     downstream of search_path, with the archive path, whether reads were banked,
#'     and (where obtainable) the .jams file modification date and owner.
#'   - "query": pass a vector of sample prefixes to query_samples and receive back a
#'     data frame (Sample, Banked_in_archive TRUE/FALSE, plus archive path and details
#'     for those found). Useful for vetting a JAMSbeta sample list before launching.
#'
#' @param search_path Absolute (or relative) path to a directory to search downward from. May contain several JAMSarchive directories at any depth.
#' @param query_samples Optional character vector of sample prefixes to check for. If supplied, mode is forced to "query".
#' @param mode Either "inventory" or "query". Ignored (set to "query") if query_samples is given. Default "inventory".
#' @param recursive Logical. Search subdirectories recursively. Default TRUE.
#' @param flag_duplicates Logical. If TRUE (default), warn when the same prefix is banked in more than one archive downstream of search_path.
#' @param return_vector Logical. In "query" mode only, if TRUE return a named logical vector (names = query_samples) instead of a data frame. Default FALSE.
#' @export

check_for_banked_samples <- function(search_path = NULL, query_samples = NULL, mode = "inventory", recursive = TRUE, flag_duplicates = TRUE, return_vector = FALSE){

    #Vet input path
    if (is.null(search_path)){
        stop("You must supply a search_path to scan for banked JAMS samples.")
    }
    search_path <- fixrelpath(search_path)
    if (!dir.exists(search_path)){
        stop(paste("The search_path supplied does not exist:", search_path))
    }

    #If query_samples is supplied, we are in query mode regardless of what was asked.
    if (!is.null(query_samples)){
        mode <- "query"
        query_samples <- unique(as.character(query_samples))
    }

    if (!(mode %in% c("inventory", "query"))){
        stop("mode must be either \"inventory\" or \"query\".")
    }

    #############################################################################
    ## 1. Find every .jams file and every JAMSalpha tarball downstream.
    ## We do NOT assume a fixed jamsfiles/ or JAMStarballs/ layout. Instead we
    ## locate the files wherever they are, and treat the directory that CONTAINS
    ## each .jams file as a stand-in for its JAMSarchive. This tolerates archives
    ## at different depths and multiple projects under one search_path.
    #############################################################################
    flog.info(paste("Scanning for banked JAMS samples under", search_path))

    jamsfiles <- list.files(path = search_path, pattern = "\\.jams$", full.names = TRUE, recursive = recursive, include.dirs = FALSE)
    tarballs  <- list.files(path = search_path, pattern = "_JAMS\\.tar\\.gz$", full.names = TRUE, recursive = recursive, include.dirs = FALSE)

    if (length(jamsfiles) < 1){
        flog.warn(paste("No .jams files were found under", search_path))
        #Return an empty-but-correctly-shaped object so downstream code does not break.
        if (mode == "query"){
            if (return_vector){
                res <- setNames(rep(FALSE, length(query_samples)), query_samples)
                return(res)
            } else {
                return(data.frame(Sample = query_samples, Banked_in_archive = FALSE, Archive_path = NA_character_, Reads_banked = NA, Jams_file_date = as.POSIXct(NA), Banked_by = NA_character_, stringsAsFactors = FALSE))
            }
        } else {
            return(data.frame(Sample = character(0), Archive_path = character(0), Reads_banked = logical(0), Jams_file_date = as.POSIXct(character(0)), Banked_by = character(0), stringsAsFactors = FALSE))
        }
    }

    #Helper: derive the prefix from a .jams filename (strip directory and .jams suffix)
    jams_prefix <- function(fp){ gsub("\\.jams$", "", basename(fp)) }
    #Helper: derive the prefix from a _JAMS.tar.gz filename
    tarball_prefix <- function(fp){ gsub("_JAMS\\.tar\\.gz$", "", basename(fp)) }

    #Build a data frame of every .jams file found, with its prefix and containing folder.
    jdf <- data.frame(Jams_fullpath = jamsfiles, Sample = sapply(jamsfiles, jams_prefix, USE.NAMES = FALSE), Jams_dir = dirname(jamsfiles), stringsAsFactors = FALSE)

    #############################################################################
    ## 2. Resolve the JAMSarchive root for each .jams file.
    ## JAMSbankit organises an archive as <JAMSarchive>/{jamsfiles,JAMStarballs,
    ## reports,reads,...}. If the .jams file sits inside a folder literally called
    ## "jamsfiles", the archive is its PARENT; otherwise we take the containing
    ## directory itself as the archive root. This keeps both canonical and
    ## flattened layouts working.
    #############################################################################
    resolve_archive_root <- function(jams_dir){
        if (basename(jams_dir) == "jamsfiles"){
            return(dirname(jams_dir))
        } else {
            return(jams_dir)
        }
    }
    jdf$Archive_path <- sapply(jdf$Jams_dir, resolve_archive_root, USE.NAMES = FALSE)

    #############################################################################
    ## 3. Match each .jams to a JAMSalpha tarball to confirm it is truly "banked".
    ## A sample counts as banked only if BOTH exist within the SAME archive root.
    ## We resolve each tarball's archive root the same way (parent of a
    ## JAMStarballs/ folder, else its own directory), then match on
    ## Archive_path + prefix.
    #############################################################################
    if (length(tarballs) > 0){
        tdf <- data.frame(Tarball_fullpath = tarballs, Sample = sapply(tarballs, tarball_prefix, USE.NAMES = FALSE), Tar_dir = dirname(tarballs), stringsAsFactors = FALSE)
        resolve_archive_root_tar <- function(tar_dir){
            if (basename(tar_dir) == "JAMStarballs"){
                return(dirname(tar_dir))
            } else {
                return(tar_dir)
            }
        }
        tdf$Archive_path <- sapply(tdf$Tar_dir, resolve_archive_root_tar, USE.NAMES = FALSE)
        tdf$Archive_sample_key <- paste(tdf$Archive_path, tdf$Sample, sep = "§")
        tarball_keys <- unique(tdf$Archive_sample_key)
    } else {
        flog.warn("No _JAMS.tar.gz tarballs were found; no sample can be confirmed as fully banked (a .jams without its tarball is incomplete).")
        tarball_keys <- character(0)
    }

    jdf$Archive_sample_key <- paste(jdf$Archive_path, jdf$Sample, sep = "§")
    jdf$Has_tarball <- jdf$Archive_sample_key %in% tarball_keys

    #Keep only samples that have both a .jams and a tarball in the same archive.
    banked <- jdf[jdf$Has_tarball, , drop = FALSE]

    if (nrow(banked) < 1){
        flog.warn("Found .jams file(s) but none had a matching JAMSalpha tarball in the same archive. Nothing is confirmed as fully banked.")
    }

    #############################################################################
    ## 4. Enrich with reads-banked status, file date and owner.
    ## Reads: JAMSbankit places banked reads in a reads/ subfolder of the archive,
    ## named <prefix>_*reads*.fastq.gz (Microbiota or QCtrimmed). We look for any
    ## such file for the prefix. Date and owner come from file.info() on the .jams.
    #############################################################################
    reads_present <- function(archive_path, prefix){
        reads_dir <- file.path(archive_path, "reads")
        if (!dir.exists(reads_dir)){
            #Fall back: look anywhere in the archive for this prefix's banked reads.
            hits <- list.files(path = archive_path, pattern = paste0("^", prefix, "_.*reads.*\\.fastq\\.gz$"), full.names = FALSE, recursive = TRUE)
        } else {
            hits <- list.files(path = reads_dir, pattern = paste0("^", prefix, "_.*reads.*\\.fastq\\.gz$"), full.names = FALSE, recursive = FALSE)
        }
        return(length(hits) > 0)
    }

    if (nrow(banked) > 0){
        finfo <- file.info(banked$Jams_fullpath)
        banked$Jams_file_date <- finfo$mtime
        #file.info() may or may not populate uname depending on platform; guard it.
        if ("uname" %in% colnames(finfo)){
            banked$Banked_by <- finfo$uname
        } else {
            banked$Banked_by <- NA_character_
        }
        banked$Banked_by[is.na(banked$Banked_by)] <- "unknown"
        banked$Reads_banked <- mapply(reads_present, banked$Archive_path, banked$Sample, USE.NAMES = FALSE)
    }

    #############################################################################
    ## 5. Flag duplicate prefixes across archives.
    ## A prefix appearing in more than one archive is bad practice (double banking
    ## or cross-project name clashes) and is surfaced, not silently collapsed.
    #############################################################################
    if (flag_duplicates && nrow(banked) > 0){
        dup_prefixes <- unique(banked$Sample[duplicated(banked$Sample)])
        if (length(dup_prefixes) > 0){
            flog.warn(paste("The following prefix(es) are banked in MORE THAN ONE archive downstream of the search path. This is bad practice (double-banking or cross-project name collisions) and you should disambiguate:", paste0(dup_prefixes, collapse = ", ")))
            for (dp in dup_prefixes){
                where <- banked$Archive_path[banked$Sample == dp]
                flog.warn(paste0("   '", dp, "' found in: ", paste0(where, collapse = " ; ")))
            }
        }
    }

    #Assemble the tidy inventory data frame.
    inventory <- banked[ , c("Sample", "Archive_path", "Reads_banked", "Jams_file_date", "Banked_by"), drop = FALSE]
    inventory <- inventory[order(inventory$Sample, inventory$Archive_path), ]
    rownames(inventory) <- NULL

    #############################################################################
    ## 6. Return according to mode.
    #############################################################################
    if (mode == "inventory"){
        flog.info(paste("Found", length(unique(inventory$Sample)), "unique banked sample prefix(es) across", length(unique(inventory$Archive_path)), "archive(s)."))
        return(inventory)
    }

    #Query mode: report, for each requested sample, whether it is banked anywhere.
    banked_set <- unique(inventory$Sample)
    is_banked <- query_samples %in% banked_set

    if (return_vector){
        res <- setNames(is_banked, query_samples)
        missing_q <- query_samples[!is_banked]
        if (length(missing_q) > 0){
            flog.warn(paste("The following queried sample(s) are NOT banked under the search path:", paste0(missing_q, collapse = ", ")))
        }
        return(res)
    }

    #Data frame form: one row per queried sample. For banked ones, attach the
    #first archive where found (and its details); for duplicates, note the count.
    query_df <- data.frame(Sample = query_samples, Banked_in_archive = is_banked, stringsAsFactors = FALSE)

    #Join details for the banked ones (taking the first archive per sample, and a count).
    first_hit <- inventory[!duplicated(inventory$Sample), , drop = FALSE]
    rownames(first_hit) <- first_hit$Sample
    archive_counts <- table(inventory$Sample)

    query_df$Archive_path   <- ifelse(is_banked, first_hit[query_df$Sample, "Archive_path"], NA_character_)
    query_df$Num_archives   <- ifelse(is_banked, as.integer(archive_counts[query_df$Sample]), 0L)
    query_df$Reads_banked   <- ifelse(is_banked, first_hit[query_df$Sample, "Reads_banked"], NA)
    query_df$Jams_file_date <- as.POSIXct(ifelse(is_banked, first_hit[query_df$Sample, "Jams_file_date"], NA), origin = "1970-01-01")
    query_df$Banked_by      <- ifelse(is_banked, first_hit[query_df$Sample, "Banked_by"], NA_character_)

    missing_q <- query_samples[!is_banked]
    if (length(missing_q) > 0){
        flog.warn(paste("The following queried sample(s) are NOT banked under the search path:", paste0(missing_q, collapse = ", ")))
    }
    rownames(query_df) <- NULL

    return(query_df)
}