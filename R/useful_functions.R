#' log2PPMtoPct(log2PPM)
#'
#' Returns a rounded percentage given a log2 transformed PPM value
#' @export

log2PPMtoPct <- function(log2PPM= NULL, signifdigits = 1){
    PPM <- ((2 ^ log2PPM) - 1)
    Pct <- round((PPM / 10000), signifdigits)

    return(Pct)
}

#' Pct2log2PPM(Pct)
#'
#' Returns a log2 transformed PPM value given a percentage
#' @export

Pct2log2PPM <- function(Pct){
    PPM <- Pct * 10000
    log2PPM <- log2(PPM + 1)

    return(log2PPM)
}

#' is.redundant(vec)
#'
#' Returns logical for a vector being redundant, i.e. has at least one element repeated in it.
#' @export
is.redundant <- function(vec){
    propunique <- length(unique(vec)) / length(vec)
    if (propunique < 1){
        redundant <- TRUE
    } else {
        redundant <- FALSE
    }

    return(redundant)
}


#' is.useful(vec)
#'
#' Returns logical for a vector having more than one class in it
#' @export

is.useful <- function(vec){
    numunique <- length(unique(vec))
    if (numunique < 2){
        useful <- FALSE
    } else {
        useful <- TRUE
    }

    return(useful)
}


#' convert_matrix_log2(mat = NULL, transformation = NULL)
#'
#' Wrapper for transforming a matrix to and from log2 formats
#' @export

convert_matrix_log2 <- function(mat = NULL, transformation = NULL){

    if (!(transformation %in% c("to_log2", "from_log2"))){
        stop("Please choose either \"to_log2\", or \"from_log2\" as transformation to apply to matrix.")
    }

    if (class(mat)[1] != "matrix"){
        stop("Object to transform must be a matrix.")
    }

    if (transformation == "to_log2"){
        tfun <- function (x) {
            fx <- log2(x + 1)
            return(fx)
        }
    } else if (transformation == "from_log2"){
        tfun <- function (x) {
            fx <- ((2 ^ x) - 1)
            return(fx)
        }
    }

    transmat <- base::apply(mat, MARGIN = 2, FUN = tfun)

    return(transmat)

}


#' convert_matrix_PA(mat = NULL, threshPA = 0)
#'
#' Wrapper for transforming a matrix to and from log2 formats
#' @export

convert_matrix_PA <- function(mat = NULL, threshPA = 0){

    if (class(mat)[1] != "matrix"){
        stop("Object to transform must be a matrix.")
    }

    PAvec <- NULL
    for (rownum in 1:nrow(mat)){
        PAvec <- mat[rownum, ]
        if (threshPA != 0){
            PAvec[which(PAvec < threshPA)] <- 0
            PAvec[which(PAvec >= threshPA)] <- 1
        } else {
            PAvec[which(PAvec > 0)] <- 1
        }
        mat[rownum, ] <- PAvec
    }

    return(mat)

}


#' detectHardwareResources()
#' Generic function for getting the number of CPUs and the total memory available to JAMS.
#' Works on laptops (macOS Intel/Apple silicon, Linux, Windows), on stand-alone servers and on
#' HPC clusters. Base R only (no benchmarkme), so it also works on Windows.
#'
#' Resolution order for CPUs: JAMS_THREADS env var (manual override) > Slurm
#' (SLURM_CPUS_PER_TASK, SLURM_CPUS_ON_NODE, sacct) > PBS/Torque (PBS_NUM_PPN, NCPUS) >
#' SGE/UGE (NSLOTS) > LSF (LSB_DJOB_NUMPROC) > parallel::detectCores().
#' Resolution order for memory: JAMS_MEMORY_GB env var (manual override) > Slurm
#' (SLURM_MEM_PER_NODE, SLURM_MEM_PER_CPU x CPUs) > Linux cgroup limit (containers, Slurm with
#' cgroups, Kubernetes) > physical RAM. A scheduler/cgroup value is never allowed to exceed physical RAM.
#'
#' Returns a named numeric vector with available CPUs ("threads") and memory in bytes ("memory").
#' @export

detectHardwareResources <- function(){

    env_num <- function(var){
        val <- Sys.getenv(var, unset = "")
        if (!nzchar(val)){
            return(NA_real_)
        }
        #Slurm can report things like "16(x2)"; keep the leading integer.
        suppressWarnings(as.numeric(sub("^([0-9.]+).*$", "\\1", val)))
    }

    read_first_line <- function(path){
        if (!file.exists(path)){
            return(NA_character_)
        }
        out <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character(0))
        if (length(out) == 0){
            return(NA_character_)
        }
        out
    }

    ## CPUs ##
    detectBatchCPUs <- function(){
        ncores <- env_num("JAMS_THREADS")

        if (is.na(ncores) && nzchar(Sys.getenv("SLURM_JOB_ID"))){
            ncores <- env_num("SLURM_CPUS_PER_TASK")
            if (is.na(ncores)){
                ncores <- env_num("SLURM_CPUS_ON_NODE")
            }
            if (is.na(ncores) && nzchar(Sys.which("sacct"))){
                #Last resort for Slurm: ask the accounting database.
                sacctraw <- tryCatch(suppressWarnings(system2("sacct", args = c("-j", Sys.getenv("SLURM_JOB_ID"), "-X", "-P", "-n", "-o", "AllocCPUS"), stdout = TRUE, stderr = FALSE, timeout = 20)), error = function(e) character(0))
                if (length(sacctraw) > 0){
                    ncores <- suppressWarnings(as.numeric(trimws(sacctraw[1])))
                }
            }
        }

        #Other common schedulers
        for (schedvar in c("PBS_NUM_PPN", "NCPUS", "NSLOTS", "LSB_DJOB_NUMPROC")){
            if (is.na(ncores)){
                ncores <- env_num(schedvar)
            }
        }

        if (is.na(ncores)){
            ncores <- parallel::detectCores(logical = TRUE)
        }

        if (is.na(ncores) || (ncores < 1)){
            warning("Could not determine how many CPUs you have. Defaulting to 1. Set the JAMS_THREADS environment variable to override.")
            ncores <- 1
        }

        return(as.numeric(floor(ncores)))
    }

    ## Memory ##
    #Physical RAM of the machine, in bytes, using only base R and OS-native tools.
    physical_ram_bytes <- function(){
        sysname <- Sys.info()[["sysname"]]
        totmembytes <- NA_real_

        if (sysname == "Linux"){
            if (file.exists("/proc/meminfo")){
                meminfo <- readLines("/proc/meminfo", warn = FALSE)
                memtotal <- grep("^MemTotal:", meminfo, value = TRUE)
                if (length(memtotal) == 1){
                    totmembytes <- as.numeric(gsub("[^0-9]", "", memtotal)) * 1024
                }
            }
        } else if (sysname == "Darwin"){
            out <- tryCatch(suppressWarnings(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
            if (length(out) > 0){
                totmembytes <- suppressWarnings(as.numeric(out[1]))
            }
        } else if (.Platform$OS.type == "windows"){
            #wmic is gone from recent Windows 11 builds (this is what broke benchmarkme), so ask PowerShell first.
            out <- tryCatch(suppressWarnings(system2("powershell", c("-NoProfile", "-NonInteractive", "-Command", shQuote("(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory", type = "cmd")), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
            out <- suppressWarnings(as.numeric(trimws(out)))
            out <- out[!is.na(out)]
            if (length(out) > 0){
                totmembytes <- out[1]
            } else {
                out <- tryCatch(suppressWarnings(system2("wmic", c("ComputerSystem", "get", "TotalPhysicalMemory", "/value"), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
                out <- grep("TotalPhysicalMemory=", out, value = TRUE)
                if (length(out) > 0){
                    totmembytes <- suppressWarnings(as.numeric(gsub("[^0-9]", "", out[1])))
                }
            }
        } else {
            #Other unixes (e.g. FreeBSD)
            out <- tryCatch(suppressWarnings(system2("sysctl", c("-n", "hw.physmem"), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
            if (length(out) > 0){
                totmembytes <- suppressWarnings(as.numeric(out[1]))
            }
        }

        return(totmembytes)
    }

    #Memory limit imposed by Linux control groups (Docker/Apptainer, Kubernetes, Slurm with cgroups), in bytes.
    cgroup_limit_bytes <- function(){
        if (Sys.info()[["sysname"]] != "Linux" || !file.exists("/proc/self/cgroup")){
            return(NA_real_)
        }
        cg <- tryCatch(readLines("/proc/self/cgroup", warn = FALSE), error = function(e) character(0))
        limits <- NA_real_
        #cgroup v2 (unified hierarchy): a line like "0::/path"
        v2 <- grep("^0::", cg, value = TRUE)
        if (length(v2) > 0){
            cgpath <- sub("^0::", "", v2[1])
            #Walk up the hierarchy: the effective limit is the smallest one on the way to the root.
            repeat {
                lim <- read_first_line(file.path("/sys/fs/cgroup", cgpath, "memory.max"))
                if (!is.na(lim) && (lim != "max")){
                    limits <- c(limits, suppressWarnings(as.numeric(lim)))
                }
                if (cgpath %in% c("/", "", ".")){
                    break
                }
                cgpath <- dirname(cgpath)
            }
        }
        #cgroup v1: a line like "N:memory:/path"
        v1 <- grep(":memory:", cg, value = TRUE, fixed = TRUE)
        if (length(v1) > 0){
            cgpath <- sub("^[0-9]+:memory:", "", v1[1])
            for (cand in unique(c(file.path("/sys/fs/cgroup/memory", cgpath, "memory.limit_in_bytes"), "/sys/fs/cgroup/memory/memory.limit_in_bytes"))){
                lim <- suppressWarnings(as.numeric(read_first_line(cand)))
                limits <- c(limits, lim)
            }
        }
        limits <- limits[!is.na(limits) & (limits > 0)]
        #cgroup v1 reports "unlimited" as a huge number (~9.2e18); ignore anything above 1 EiB.
        limits <- limits[limits < 2^60]
        if (length(limits) == 0){
            return(NA_real_)
        }

        return(min(limits))
    }

    detectAvailRAM <- function(){
        physmem <- physical_ram_bytes()

        #Manual override
        override_gb <- env_num("JAMS_MEMORY_GB")
        if (!is.na(override_gb)){
            return(override_gb * 1e9)
        }

        candidates <- physmem
        if (nzchar(Sys.getenv("SLURM_JOB_ID"))){
            #Slurm reports memory in megabytes.
            mempernode <- env_num("SLURM_MEM_PER_NODE")
            if (!is.na(mempernode)){
                candidates <- c(candidates, mempernode * 1e6)
            } else {
                mempercpu <- env_num("SLURM_MEM_PER_CPU")
                if (!is.na(mempercpu)){
                    slurmcpus <- env_num("SLURM_CPUS_PER_TASK")
                    if (is.na(slurmcpus)){
                        slurmcpus <- env_num("SLURM_CPUS_ON_NODE")
                    }
                    if (is.na(slurmcpus)){
                        slurmcpus <- 1
                    }
                    candidates <- c(candidates, mempercpu * slurmcpus * 1e6)
                }
            }
        }
        candidates <- c(candidates, cgroup_limit_bytes())
        candidates <- candidates[!is.na(candidates) & (candidates > 0)]

        if (length(candidates) == 0){
            warning("Could not determine how much RAM you have. Set the JAMS_MEMORY_GB environment variable to tell JAMS.")
            return(NA_real_)
        }

        return(as.numeric(min(candidates)))
    }

    hardwareRes <- c(threads = detectBatchCPUs(), memory = detectAvailRAM())

    return(hardwareRes)
}

#' whoopsieplot(msg = NULL)
#' Shuts down the device and gives a message on error, useful for when making reports.
#' @export
whoopsieplot <- function(msg = "trying to do this."){
    flog.info(paste("Whoops, something went wrong while", msg))
    dev.off()
}


#' filetype(path)
#' Wrapper for returning class of system file.
#' @export
filetype <- function(path){
    f <- file(path)
    ext <- summary(f)$class
    close.connection(f)
    ext
}


#' trim_whitespace_from_df(df)
#' Trims leading and trailing whitespace from a dataframe
#' @export

trim_whitespace_from_df <- function(df = NULL) {
    # Is it a data frame, or its ilk?
    if (!inherits(df, "data.frame")) {
        stop("Input must be a data.frame or inherit from one (e.g., tibble).")
    }

    # Is the data frame empty?
    if (ncol(df) == 0) {
        return(df)
    }

    # Ok, now we can fix each column.
    for (colm in seq_len(ncol(df))) {
        # Process character or factor columns
        if (is.character(df[[colm]])) {
            df[[colm]] <- trimws(df[[colm]], whitespace = "[\\s\\h\\v]")
        } else if (is.factor(df[[colm]])) {
            levels(df[[colm]]) <- trimws(levels(df[[colm]]), whitespace = "[\\s\\h\\v]")
        }
    }

    return(df)
}


#' replace_NAs_with_character(df = NULL, replacement = "N_A", coerce_all_to_character = TRUE, replace_empty_strings = TRUE)
#'
#' Conforms a data frame to the JAMS metadata convention: every cell is of class character,
#' and every missing cell holds a single explicit sentinel (by default "N_A"). This is what
#' allows omission to be controlled downstream solely through the class_to_ignore argument,
#' including for numeric variables, whose values travel as "231", "56", ... and are coerced
#' back with can_be_made_numeric() only once the ignored classes have been dropped.
#'
#' @param df A data frame, or anything inheriting from one (tibble, data.table).
#' @param replacement Character scalar written into every missing cell. Default "N_A". If a
#'   vector is passed (class_to_ignore is allowed to hold several values elsewhere in JAMS),
#'   only its first element is used, as there can be only one canonical sentinel.
#' @param coerce_all_to_character Logical. If TRUE (default), EVERY column is returned as
#'   character, not just the ones which contained NAs. Set to FALSE for the historical
#'   behaviour of only touching columns which hold at least one NA.
#' @param replace_empty_strings Logical. If TRUE (default), character cells which are empty or
#'   contain only whitespace are treated as missing and also receive the replacement. Note that
#'   with coerce_all_to_character = FALSE this can only apply to columns which are already
#'   character.
#'
#' @export

replace_NAs_with_character <- function(df = NULL, replacement = "N_A", coerce_all_to_character = TRUE, replace_empty_strings = TRUE){
    #Accept anything inheriting from data.frame, not just a bare data.frame. A data.table in
    #particular turns up whenever metadata was read with fread() without data.table = FALSE,
    #and its `[<-` semantics differ from a data.frame's, so flatten it before touching anything.
    if (!inherits(df, "data.frame")){
        stop("Input must be a data.frame or inherit from one (e.g., a tibble or a data.table).")
    }
    if (!identical(class(df), "data.frame")){
        original_rownames <- rownames(df)
        df <- as.data.frame(df, stringsAsFactors = FALSE)
        rownames(df) <- original_rownames
    }

    #There can be only one sentinel. Callers routinely hand their class_to_ignore straight to
    #this argument, and that is allowed to be a vector of several values, which would otherwise
    #be recycled across the missing cells of a column.
    replacement <- as.character(replacement)
    replacement <- replacement[!is.na(replacement)]
    if (length(replacement) < 1){
        stop("replacement must contain at least one non-missing value.")
    }
    replacement <- replacement[1]

    if (ncol(df) == 0){
        return(df)
    }

    df[] <- lapply(df, function(x){
        #Record missingness BEFORE any coercion. as.character(NaN) is the string "NaN", which
        #is no longer NA, so a NaN would otherwise survive as a phantom class.
        na_mask <- is.na(x)

        #Dates, factors, integers and logicals must all come out as plain character. A Date
        #column is the nastiest of these, because it is numeric underneath, so it passes
        #can_be_made_numeric() and is then compared against class_to_ignore by
        #define_kinds_of_variables(), where `!=` dispatches to `!=.Date` and dies.
        if (coerce_all_to_character || any(na_mask)){
            x <- as.character(x)
        }

        if (replace_empty_strings && is.character(x)){
            #fread() writes an empty tsv field as "" rather than NA, so a blank cell would
            #otherwise survive as a class which class_to_ignore never catches. %in% is used
            #rather than == so that the already-NA positions do not propagate NA into the mask.
            na_mask <- na_mask | (trimws(x, whitespace = "[\\s\\h\\v]") %in% "")
        }

        #Guarded so that, when coerce_all_to_character is FALSE, a clean numeric column is
        #handed back with its type intact rather than being silently turned into character.
        if (any(na_mask)){
            x[na_mask] <- replacement
        }

        return(x)
    })

    return(df)
}


#' enforce_JAMS_style_metadata(pheno = NULL, phenolabels = NULL, sample_column = "Sample", class_to_ignore = "N_A", fix_column_names = TRUE)
#'
#' Checks a metadata table (phenotable) and conforms it to JAMS conventions, so that it will not cause trouble in downstream JAMS functions. It is meant to be run once, right after reading metadata in and before building or updating SummarizedExperiment objects, and it reports everything it changes.
#'
#' In order, it:
#' (1) checks that pheno is a data.frame (tibbles and data.tables are flattened to a plain data.frame);
#' (2) makes column names safe: leading/trailing whitespace is removed; "%" becomes "pct" and "+" becomes "plus"; every other character that is not an ASCII letter, digit, underscore or dot (spaces, "-", "!@#$^&()", "/", accented letters, etc.) becomes "_"; runs of "_" are collapsed; leading/trailing "_" and "." are removed; names starting with a digit get an "X" prefix; R reserved words (TRUE, NA, if, function...) get a "_var" suffix; empty names become "Unnamed_column_<n>"; and names left duplicated are made unique by appending _1, _2, etc. (e.g. "IL-6 (pg/mL)" becomes "IL_6_pg_mL", "Olink_IL_1 alpha_Serum" becomes "Olink_IL_1_alpha_Serum");
#' (3) finds the Sample column. If there is no column named exactly sample_column, a single column whose name differs from it only in case (e.g. "sample", "SAMPLE") is renamed to it; otherwise it stops;
#' (4) trims leading and trailing whitespace from every cell (trim_whitespace_from_df);
#' (5) enforces the JAMS all-character convention, with every NA or empty cell set to class_to_ignore[1] (replace_NAs_with_character);
#' (6) vets the sample names. Rows where every cell is empty are dropped quietly (common at the bottom of spreadsheets); rows with a missing sample name are dropped with a warning; duplicated sample names make it stop, listing them, as there is no safe way to choose between them (or, with on_duplicate_samples = "drop", every row of a duplicated sample is dropped with a warning). Sample names are NEVER altered, as they must match the JAMS file prefixes; a warning is issued if they contain characters other than letters, digits, "_", "." or "-".
#' (7) warns if any cell contains the section sign character, which JAMS uses internally to separate subset tiers.
#'
#' Everything is reported through futile.logger, one line per renamed column, so when it is called from JAMSbeta or JAMS16 (via load_metadata_from_file) the full record of changes ends up in the project's _JAMS.log file.
#'
#' @param pheno data.frame (or tibble/data.table) with the metadata, one row per sample.
#' @param phenolabels Optional data.frame with columns Var_label and Var_type. If given, its Var_label entries are renamed in step with the pheno columns. A label is matched to its column by the column's original name or, failing that, by the dotted form older JAMS versions produced on reading (e.g. "IL.6..pg.mL." for "IL-6 (pg/mL)"), so existing phenolabels keep working. The Sample column is ensured to be the only one with Var_type "Sample". See also reconcile_phenolabels.
#' @param sample_column Name of the column holding the sample names (prefixes). Default is "Sample". If this name is itself changed by the column name fixing, the fixed name is used.
#' @param class_to_ignore String (or vector, of which the first element is used for filling) used for missing data. Default is "N_A".
#' @param fix_column_names Logical. If TRUE (the default), unsafe column names are fixed. If FALSE, they are only reported with a warning and left as they are.
#' @param reconcile_phenolabels Logical, only used if phenolabels is given. If TRUE (the default), phenolabels entries with no matching pheno column are dropped, pheno columns missing from phenolabels are added with an inferred Var_type, and phenolabels is put in pheno column order. If FALSE, phenolabels is only renamed (and its Sample entry fixed), and nothing is added, dropped or reordered; this is what load_metadata_from_file uses, because JAMSbeta treats phenolabels as the list of columns to keep.
#' @param on_duplicate_samples Either "stop" (the default) or "drop". What to do when a sample name occurs in more than one row.
#' @return If phenolabels is NULL, the fixed phenotable data.frame. If phenolabels is given, a list with elements phenotable and phenolabels. In both cases the phenotable carries an attribute "JAMS_column_name_changes", a data.frame with columns Original and New listing every column that was renamed (zero rows if none), retrievable with attr(x, "JAMS_column_name_changes"). It also carries the attributes "JAMS_sample_column" (the name of the Sample column after fixing) and "JAMS_duplicated_samples_dropped" (sample names dropped with on_duplicate_samples = "drop").
#' @examples
#' \dontrun{
#' pheno <- enforce_JAMS_style_metadata(pheno = pheno)
#' attr(pheno, "JAMS_column_name_changes")
#' md <- enforce_JAMS_style_metadata(pheno = pheno, phenolabels = phenolabels)
#' }
#' @export

enforce_JAMS_style_metadata <- function(pheno = NULL, phenolabels = NULL, sample_column = "Sample", class_to_ignore = "N_A", fix_column_names = TRUE, reconcile_phenolabels = TRUE, on_duplicate_samples = "stop"){

    #### 1. Is it a data frame? ####
    if (is.null(pheno) || !inherits(pheno, "data.frame")){
        stop("pheno must be a data.frame (or a tibble or data.table).")
    }
    if (!(is.character(on_duplicate_samples) && (length(on_duplicate_samples) == 1) && (on_duplicate_samples %in% c("stop", "drop")))){
        stop("on_duplicate_samples must be either \"stop\" or \"drop\".")
    }
    #Flatten tibble/data.table; optional = TRUE so that as.data.frame() does not itself mangle the names with make.names()
    pheno <- as.data.frame(pheno, optional = TRUE)
    if ((ncol(pheno) == 0) || (nrow(pheno) == 0)){
        stop("pheno has no columns or no rows.")
    }
    flog.info(paste0("enforce_JAMS_style_metadata: checking metadata with ", nrow(pheno), " row(s) and ", ncol(pheno), " column(s) for conformity with JAMS style."))
    nrow_in <- nrow(pheno)
    if (!(is.character(class_to_ignore) && (length(class_to_ignore) >= 1) && !is.na(class_to_ignore[1]))){
        stop("class_to_ignore must be a character string (or vector of strings), e.g. \"N_A\".")
    }
    if (!(is.character(sample_column) && (length(sample_column) == 1) && !is.na(sample_column))){
        stop("sample_column must be a single column name.")
    }

    #### 2. Column names ####
    orignames <- colnames(pheno)
    if (is.null(orignames)){
        orignames <- rep("", ncol(pheno))
    }
    orignames[is.na(orignames)] <- ""

    sanitize_one_name <- function(x, i){
        #To ASCII first: every non-ASCII byte becomes "_", which also sidesteps any invalid-encoding errors in gsub()
        y <- iconv(enc2utf8(x), from = "UTF-8", to = "ASCII", sub = "_")
        if (is.na(y)){
            y <- ""
        }
        y <- trimws(y)
        y <- gsub("%", "pct", y, fixed = TRUE)
        y <- gsub("+", "plus", y, fixed = TRUE)
        y <- gsub("[^A-Za-z0-9_.]", "_", y, perl = TRUE)
        y <- gsub("_+", "_", y, perl = TRUE)
        y <- gsub("^[_.]+|[_.]+$", "", y, perl = TRUE)
        if (y == ""){
            y <- paste0("Unnamed_column_", i)
        }
        if (grepl("^[0-9]", y)){
            y <- paste0("X", y)
        }
        #What survives is syntactic except for reserved words (TRUE, NA, if, function, ...), which make.names() alters
        if (make.names(y) != y){
            y <- paste0(y, "_var")
        }
        return(y)
    }

    proposed <- vapply(seq_along(orignames), function(i){ sanitize_one_name(orignames[i], i) }, character(1))
    if (anyDuplicated(proposed)){
        dupnames <- unique(proposed[duplicated(proposed)])
        flog.warn(paste0("enforce_JAMS_style_metadata: after making them safe, these column names occur more than once, and will be made unique with a numeric suffix (_1, _2, ...): ", paste0(dupnames, collapse = ", "), ". Please check which column is which."))
        proposed <- make.unique(proposed, sep = "_")
    }

    if (fix_column_names){
        colnames(pheno) <- proposed
    } else {
        unsafe <- which(orignames != proposed)
        if (length(unsafe) > 0){
            flog.warn(paste0("enforce_JAMS_style_metadata: ", length(unsafe), " column name(s) are unsafe for JAMS but fix_column_names = FALSE, so they were left as they are: ", paste0("\"", orignames[unsafe], "\" (suggested: \"", proposed[unsafe], "\")", collapse = "; ")))
        }
        colnames(pheno) <- orignames
    }
    #Finds which ORIGINAL column a name refers to: the name as written, else the dotted forms older JAMS versions
    #produced on reading (read.table's make.names, with or without unique = TRUE; read.xlsx's spaces-to-dots).
    #Used both for the Sample column and for phenolabels, so the two always agree.
    match_to_original_column <- function(labels){
        tn <- trimws(orignames)
        idx <- match(labels, tn)
        for (alt in list(make.names(tn), make.names(tn, unique = TRUE), gsub(" ", ".", tn, fixed = TRUE))){
            idx[is.na(idx)] <- match(labels[is.na(idx)], alt)
        }
        return(idx)
    }

    #If the requested Sample column name was itself unsafe (or given in an old dotted form), follow it to its current name
    if (!(sample_column %in% colnames(pheno))){
        m <- match_to_original_column(sample_column)
        if (!is.na(m)){
            sample_column <- colnames(pheno)[m]
        }
    }

    #### 3. Sample column ####
    if (!(sample_column %in% colnames(pheno))){
        candidates <- which(tolower(colnames(pheno)) == tolower(sample_column))
        if (length(candidates) == 1){
            flog.warn(paste0("enforce_JAMS_style_metadata: no column called \"", sample_column, "\"; renaming column \"", colnames(pheno)[candidates], "\" to \"", sample_column, "\"."))
            colnames(pheno)[candidates] <- sample_column
        } else if (length(candidates) > 1){
            stop(paste0("There is no column called \"", sample_column, "\", and several columns differ from it only in case (", paste0(colnames(pheno)[candidates], collapse = ", "), "). Rename the right one to \"", sample_column, "\"."))
        } else {
            hints <- grep("sample", colnames(pheno), ignore.case = TRUE, value = TRUE)
            stop(paste0("No \"", sample_column, "\" column found. JAMS needs a column holding the sample names (the prefixes used for the JAMS files). ", if (length(hints) > 0) paste0("Possible candidates: ", paste0(hints, collapse = ", "), ". ") else "", "Rename the right column to \"", sample_column, "\", or pass its name with the sample_column argument."))
        }
    }

    #Record every rename (from steps 2 and 3)
    changed <- which(orignames != colnames(pheno))
    namechanges <- data.frame(Original = orignames[changed], New = colnames(pheno)[changed], stringsAsFactors = FALSE)
    if (nrow(namechanges) > 0){
        flog.warn(paste0("enforce_JAMS_style_metadata: ", nrow(namechanges), " of ", ncol(pheno), " metadata column name(s) were changed to conform to JAMS style (no spaces, hyphens or other special characters). The original -> new names follow, one per line. Use the NEW names in all JAMS functions (compareby, subsetby, colourby, etc.)."))
        for (nc in seq_len(nrow(namechanges))){
            flog.warn(paste0("enforce_JAMS_style_metadata: metadata column renamed: \"", namechanges$Original[nc], "\" -> \"", namechanges$New[nc], "\""))
        }
    } else {
        flog.info("enforce_JAMS_style_metadata: all metadata column names already conform to JAMS style; none were changed.")
    }

    #### 4 and 5. Trim whitespace, then enforce all-character with explicit missing-data sentinel ####
    pheno <- trim_whitespace_from_df(df = pheno)
    pheno <- replace_NAs_with_character(df = pheno, replacement = class_to_ignore[1])

    #### 6. Sample names ####
    is_missing <- function(v){ v %in% c(class_to_ignore, "") }
    all_empty_rows <- which(apply(pheno, 1, function(r){ all(is_missing(r)) }))
    if (length(all_empty_rows) > 0){
        flog.info(paste0("enforce_JAMS_style_metadata: dropping ", length(all_empty_rows), " completely empty row(s)."))
        pheno <- pheno[-all_empty_rows, , drop = FALSE]
    }
    no_name_rows <- which(is_missing(pheno[[sample_column]]))
    if (length(no_name_rows) > 0){
        flog.warn(paste0("enforce_JAMS_style_metadata: dropping ", length(no_name_rows), " row(s) with no sample name in column \"", sample_column, "\" (row number(s) ", paste0(utils::head(no_name_rows, 20), collapse = ", "), if (length(no_name_rows) > 20) ", ..." else "", " after removing empty rows)."))
        pheno <- pheno[-no_name_rows, , drop = FALSE]
    }
    if (nrow(pheno) == 0){
        stop("No rows with a sample name are left in pheno.")
    }
    sampnames <- pheno[[sample_column]]
    dupsamps <- unique(sampnames[duplicated(sampnames)])
    if (length(dupsamps) > 0){
        if (on_duplicate_samples == "drop"){
            flog.warn(paste0("enforce_JAMS_style_metadata: these sample names appear more than once in column \"", sample_column, "\", and ALL their rows have been dropped: ", paste0(dupsamps, collapse = ", "), ". If you want these samples, fix the metadata so each has exactly one row."))
            pheno <- pheno[!(pheno[[sample_column]] %in% dupsamps), , drop = FALSE]
            if (nrow(pheno) == 0){
                stop("No rows are left in pheno after dropping duplicated samples.")
            }
            sampnames <- pheno[[sample_column]]
        } else {
            stop(paste0("These sample names appear more than once in column \"", sample_column, "\": ", paste0(dupsamps, collapse = ", "), ". Each sample must have exactly one row. Fix the metadata and try again."))
        }
    }
    oddsamps <- sampnames[grepl("[^A-Za-z0-9_.-]", sampnames, perl = TRUE)]
    if (length(oddsamps) > 0){
        flog.warn(paste0("enforce_JAMS_style_metadata: ", length(oddsamps), " sample name(s) contain characters other than letters, digits, \"_\", \".\" or \"-\" (e.g. ", paste0("\"", utils::head(oddsamps, 5), "\"", collapse = ", "), "). They have NOT been changed, because sample names must match the JAMS file prefixes, but such names are likely to cause trouble."))
    }
    rownames(pheno) <- NULL

    #### 7. Reserved JAMS separator ####
    section_cols <- colnames(pheno)[vapply(pheno, function(v){ any(grepl("\u00a7", v, fixed = TRUE)) }, logical(1))]
    if (length(section_cols) > 0){
        flog.warn(paste0("enforce_JAMS_style_metadata: the character \"\u00a7\" occurs in column(s) ", paste0(section_cols, collapse = ", "), ". JAMS uses it internally to separate subset tiers, so subsetting by these columns will misbehave. Please replace it."))
    }

    attr(pheno, "JAMS_column_name_changes") <- namechanges
    attr(pheno, "JAMS_duplicated_samples_dropped") <- if ((length(dupsamps) > 0) && (on_duplicate_samples == "drop")) dupsamps else character(0)
    attr(pheno, "JAMS_sample_column") <- sample_column
    flog.info(paste0("enforce_JAMS_style_metadata: done. Metadata now has ", nrow(pheno), " row(s) (", nrow_in - nrow(pheno), " dropped) and ", ncol(pheno), " column(s) (", nrow(namechanges), " renamed). Sample column is \"", sample_column, "\"."))

    #### Optional: keep phenolabels in step ####
    if (!is.null(phenolabels)){
        if (!inherits(phenolabels, "data.frame") || !all(c("Var_label", "Var_type") %in% colnames(phenolabels))){
            stop("phenolabels must be a data.frame with columns Var_label and Var_type.")
        }
        phenolabels <- as.data.frame(phenolabels, optional = TRUE)
        phenolabels[] <- lapply(phenolabels, as.character)
        #Blank rows (e.g. from an Excel sheet read with skipEmptyRows = FALSE) carry no information; drop them
        blank_labels <- is.na(phenolabels$Var_label) | (phenolabels$Var_label %in% c(class_to_ignore, ""))
        if (any(blank_labels)){
            flog.info(paste0("enforce_JAMS_style_metadata: dropping ", sum(blank_labels), " phenolabels row(s) with an empty Var_label."))
            phenolabels <- phenolabels[!blank_labels, , drop = FALSE]
        }
        phenolabels <- trim_whitespace_from_df(df = phenolabels)

        #Rename labels in step with the columns: first by the original column name, then by the dotted forms
        #that older JAMS versions produced on reading (read.table check.names / read.xlsx sep.names), so that
        #phenolabels written against those dotted names keep working.
        old_labels <- phenolabels$Var_label
        idx <- match_to_original_column(phenolabels$Var_label)
        phenolabels$Var_label[!is.na(idx)] <- colnames(pheno)[idx[!is.na(idx)]]
        relabelled <- which(old_labels != phenolabels$Var_label)
        for (rl in relabelled){
            flog.info(paste0("enforce_JAMS_style_metadata: phenolabels Var_label renamed: \"", old_labels[rl], "\" -> \"", phenolabels$Var_label[rl], "\""))
        }

        if (reconcile_phenolabels){
            orphan_labels <- phenolabels$Var_label[!(phenolabels$Var_label %in% colnames(pheno))]
            if (length(orphan_labels) > 0){
                flog.warn(paste0("enforce_JAMS_style_metadata: phenolabels entries with no matching pheno column were dropped: ", paste0(orphan_labels, collapse = ", ")))
                phenolabels <- phenolabels[phenolabels$Var_label %in% colnames(pheno), , drop = FALSE]
            }
            phenolabels <- phenolabels[!duplicated(phenolabels$Var_label), , drop = FALSE]

            unlabelled <- colnames(pheno)[!(colnames(pheno) %in% phenolabels$Var_label)]
            if (length(unlabelled) > 0){
                newtypes <- vapply(unlabelled, function(x){ infer_column_type(phenotable = pheno, colm = x, class_to_ignore = class_to_ignore) }, character(1))
                flog.info(paste0("enforce_JAMS_style_metadata: added ", length(unlabelled), " pheno column(s) missing from phenolabels, with inferred Var_type: ", paste0(unlabelled, " (", newtypes, ")", collapse = ", ")))
                extra <- phenolabels[rep(NA_integer_, length(unlabelled)), , drop = FALSE]
                extra$Var_label <- unlabelled
                extra$Var_type <- unname(newtypes)
                phenolabels <- rbind(phenolabels, extra)
            }
        } else {
            orphan_labels <- phenolabels$Var_label[!(phenolabels$Var_label %in% colnames(pheno))]
            if (length(orphan_labels) > 0){
                flog.warn(paste0("enforce_JAMS_style_metadata: these phenolabels entries do not match any metadata column (left as they are): ", paste0(orphan_labels, collapse = ", ")))
            }
            duplabels <- unique(phenolabels$Var_label[duplicated(phenolabels$Var_label)])
            if (length(duplabels) > 0){
                flog.warn(paste0("enforce_JAMS_style_metadata: these phenolabels entries occur more than once (after renaming); only the first of each is kept: ", paste0(duplabels, collapse = ", ")))
                phenolabels <- phenolabels[!duplicated(phenolabels$Var_label), , drop = FALSE]
            }
        }

        #Exactly one Sample, and it is sample_column
        other_samples <- which((phenolabels$Var_type %in% "Sample") & !(phenolabels$Var_label %in% sample_column))
        if (length(other_samples) > 0){
            flog.warn(paste0("enforce_JAMS_style_metadata: column(s) ", paste0(phenolabels$Var_label[other_samples], collapse = ", "), " were declared with Var_type \"Sample\" in phenolabels; changed to \"discrete\", as \"", sample_column, "\" is the Sample column."))
            phenolabels$Var_type[other_samples] <- "discrete"
        }
        phenolabels$Var_type[phenolabels$Var_label %in% sample_column] <- "Sample"

        if (reconcile_phenolabels){
            #Same order as the pheno columns
            phenolabels <- phenolabels[match(colnames(pheno), phenolabels$Var_label), , drop = FALSE]
        }
        rownames(phenolabels) <- NULL

        return(list(phenotable = pheno, phenolabels = phenolabels))
    }

    return(pheno)
}


#' exporttabletsv(dataobj = NULL, project = NULL, basename = NULL, row.names = TRUE, col.names = TRUE, path = NULL)
#' Generic function to export a dataframe to tsv format
#' @export

exporttabletsv <- function(dataobj = NULL, project = NULL, basename = NULL, row.names = TRUE, col.names = TRUE, path = NULL){
    flname <- paste(paste("JAMS", project, basename, sep = "_"), "tsv", sep = ".")
    if (!is.null(path)){
        flname <- paste(path, flname, sep = "/")
    }
    write.table(dataobj, file = flname, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = row.names, col.names = col.names)
}

#' quiet(function(x))
#' Suppresses output messages
#' By Hadley Wickham
#' http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
#'
#' @export

quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}


#' can_be_made_numeric(x, cats_to_ignore = NULL)
#' Returns logical for if a vector can be coerced into numeric.
#'
#' @export

can_be_made_numeric <- function(x, cats_to_ignore = NULL){

    if (!is.null(cats_to_ignore)){
        x <- x[which(!(x %in% cats_to_ignore))]
    }
    numtest <- length(which(is.na(suppressWarnings(as.numeric(x))) == TRUE))
    if (numtest == 0){
        testresult <- TRUE
    } else {
        testresult <- FALSE
    }

    return(testresult)
}


#' IO_jams_workspace_image(opt = NULL, workspaceimage = NULL, threads = 8, operation = c("save", "load"), verbose = FALSE, verify = TRUE)
#' Safe way of either loading or saving an R workspace image (the global environment). If argument workspaceimage is null, workspace image file will be searched for in opt (opt$projimage). If that is also NULL, saving or loading is aborted. If opt is passed, number of CPUs will be set to opt$threads, trumping the threads argument.
#'
#' Saving is multi-threaded whenever the pigz program is on the PATH (no R package needed): the image is streamed through "pigz -p <threads>", which writes an ordinary gzip-compressed .RData file that base R load() reads on any OS, laptops included. Without pigz, a single-threaded base R save is used. Either way the image is first written to a temporary file next to the target and only renamed over the previous image once it is complete (and, with verify = TRUE, has passed "pigz -t"), so a crash, a full disk or a job hitting its walltime mid-save never destroys the previous image.
#' @param verify If TRUE (default), test the integrity of a pigz-compressed image before it replaces the previous one.
#'
#' @export

IO_jams_workspace_image <- function(opt = NULL, workspaceimage = NULL, threads = 8, operation = NULL, verbose = FALSE, verify = TRUE){
    if (is.null(workspaceimage)){
        workspaceimage <- opt$projimage
    }

    if (is.null(workspaceimage)){
        stop("Workspace image not found or specified.")
    }

    if (!is.null(opt)){
        if (!is.null(opt$threads)){
            threads <- opt$threads
        }
        RAMbytesavail <- opt$totmembytes
    } else {
        RAMbytesavail <- NULL
    }
    threads <- suppressWarnings(as.integer(threads))
    if ((length(threads) != 1) || is.na(threads) || (threads < 1)){
        threads <- 1L
    }

    if (is.null(operation) || !(operation %in% c("save", "load"))){
        flog.info("You must choose between \"load\" or \"save\" as an operation.")
        return(invisible(NULL))
    }

    flog.info(paste("Workspace image is", workspaceimage))

    if (verbose){
        flog.info("Random Access Memory (RAM) status")
        print(RAMbytes_status(RAMbytesavail = RAMbytesavail))
    }

    if (operation == "save"){

        save_workspace_image_safely(workspaceimage = workspaceimage, threads = threads, verify = verify)

    } else if (operation == "load"){

        #n.b. this previously called save.image() here when fastSave was absent, which would have overwritten
        #the image with the current session instead of loading it. Images written by pigz are ordinary gzip
        #streams, so base load() reads them on every OS. Load into the global environment, where they were saved from.
        flog.info("Loading project workspace image. Please be patient...")
        load(workspaceimage, envir = globalenv())

    }

    if (verbose){
        flog.info("Random Access Memory (RAM) status")
        print(RAMbytes_status(RAMbytesavail = RAMbytesavail))
    }

    invisible(workspaceimage)
}


#' save_workspace_image_safely(workspaceimage = NULL, threads = 1, verify = TRUE)
#'
#' Internal. Saves every object in the global environment (like save.image()) to workspaceimage.
#' Multi-threaded gzip compression through the pigz program when available, otherwise single-threaded base R.
#' Writes to a temporary file in the same folder and renames it over workspaceimage only on success.
#' @keywords internal
save_workspace_image_safely <- function(workspaceimage = NULL, threads = 1, verify = TRUE){

    workspaceimage <- normalizePath(workspaceimage, mustWork = FALSE)
    targetdir <- dirname(workspaceimage)
    if (!dir.exists(targetdir)){
        stop(paste("Cannot save the workspace image: folder", targetdir, "does not exist."))
    }
    #Same folder as the target, so the final rename is atomic and never crosses file systems.
    tmpimage <- file.path(targetdir, paste0(".", basename(workspaceimage), ".partial_", Sys.getpid()))
    on.exit(if (file.exists(tmpimage)) unlink(tmpimage), add = TRUE)

    objs <- ls(envir = globalenv(), all.names = TRUE)
    t0 <- Sys.time()
    saved_with <- NULL

    pigz <- Sys.which("pigz")
    if (nzchar(pigz)){
        flog.info(paste("Saving project workspace image with pigz using", threads, "CPUs. Please be patient..."))
        pigz_ok <- tryCatch({
            #pigz reads the serialized workspace from the pipe and writes a standard gzip file.
            con <- pipe(paste(shQuote(pigz), "-p", threads, ">", shQuote(tmpimage)), open = "wb")
            pipe_status <- NULL
            tryCatch(save(list = objs, file = con, envir = globalenv()),
                     finally = { pipe_status <- close(con) })
            #close() on a pipe returns the exit status of pigz (e.g. non-zero if the disk filled up).
            (is.null(pipe_status) || identical(as.integer(pipe_status), 0L)) && file.exists(tmpimage) && (file.size(tmpimage) > 0)
        }, error = function(e){
            flog.warn(paste("Multi-threaded save with pigz failed:", conditionMessage(e)))
            FALSE
        })

        if (pigz_ok && verify){
            flog.info("Verifying integrity of the saved image...")
            test_status <- suppressWarnings(system2(pigz, c("-t", shQuote(tmpimage)), stdout = FALSE, stderr = FALSE))
            pigz_ok <- identical(as.integer(test_status), 0L)
            if (!pigz_ok){
                flog.warn("The pigz-compressed image failed its integrity test.")
            }
        }

        if (pigz_ok){
            saved_with <- paste("pigz with", threads, "CPUs")
        } else {
            unlink(tmpimage)
            flog.warn("Falling back to single-threaded saving with base R. The previous image (if any) is untouched.")
        }
    } else {
        flog.info("The pigz program was not found on your PATH, so the workspace image will be saved single-threaded, which is slow for large projects.")
        flog.info("For multi-threaded saving install pigz (HPC: usually 'module load pigz'; macOS: 'brew install pigz'; Ubuntu/Debian: 'apt install pigz'; conda: 'conda install -c conda-forge pigz').")
    }

    if (is.null(saved_with)){
        flog.info("Saving project workspace image. Please be patient...")
        save(list = objs, file = tmpimage, envir = globalenv(), compress = TRUE)
        saved_with <- "base R (single-threaded)"
    }

    #Replace the previous image only now that the new one is complete.
    if (!file.rename(tmpimage, workspaceimage)){
        stop(paste("The workspace image was written to", tmpimage, "but could not be renamed to", workspaceimage, "- the previous image (if any) is untouched."))
    }

    flog.info(sprintf("Workspace image saved with %s in %s (%.2f GB).", saved_with, format(round(difftime(Sys.time(), t0), 1)), file.size(workspaceimage) / 1e9))

    invisible(workspaceimage)
}


#' RAMbytes_status(RAMbytesavail = NULL)
#' Reports how much RAM memory is being or has maximally been used.
#'
#' @export
RAMbytes_status <- function(RAMbytesavail = NULL){

    if (is.null(RAMbytesavail)){
        RAMbytesavail <- detectHardwareResources()["memory"]
    }

    memstats <- gc(full = TRUE)
    usedMbcolm <- (which(colnames(memstats) == "used") + 1)
    maxusedMbcolm <- (which(colnames(memstats) == "max used") + 1)
    usedRAMbytes <- sum(memstats[ , usedMbcolm]) * (1024 ^ 2)
    maxusedRAMbytes <- sum(memstats[ , maxusedMbcolm]) * (1024 ^ 2)
    RAMbytes <- unname(c(usedRAMbytes, maxusedRAMbytes, RAMbytesavail))
    RAMdf <- data.frame(Bytes = RAMbytes, Gbytes = round((RAMbytes / 1e9), 1), stringsAsFactors = FALSE)
    rownames(RAMdf) <- c("Used", "MaxUsed", "Available")
    RAMdf$ProportionAvail <- round((RAMdf$Bytes / RAMbytesavail), 2)

    return(RAMdf)
}


#' update_ExpObj_metadata(ExpObj = NULL, phenotable = NULL)
#'
#' Updates phenotable of a SummarizedExperiment object
#' @export

update_ExpObj_metadata <- function(ExpObj = NULL, phenotable = NULL){

        #Get appropriate object to work with
        if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
            stop("This function can only take a SummarizedExperiment object as input.")
        }

        #Let's see what we have here
        pheno_original <- colData(ExpObj)
        Samples_have <- rownames(pheno_original)
        if (!all(rownames(phenotable) %in% Samples_have)){
            stop("New phenotable contains samples not present in SummarizedExperiment. Phenotable sample names must match those in SummarizedExperiment. Check your new metadata and try again.")
        }

        #The SummarizedExperiment package does not permit coercing different metadata, so object has to be re-assembled from scratch
        assays_ExpObj <-(assays(ExpObj))
        ftt <- rowData(ExpObj)
        metadata_ExpObj <- metadata(ExpObj)

        ##Create SummarizedExperiment
        newExpObj <- SummarizedExperiment(assays = assays_ExpObj, rowData = ftt, colData = phenotable)
        metadata(newExpObj) <- metadata_ExpObj

        return(newExpObj)
}


#' fixrelpath(JAMSpath = NULL)
#'
#' Fixes path relativity and returns absolute path
#' @export

fixrelpath <- function(JAMSpath = NULL){
    require(R.utils)
    if (!(isAbsolutePath(JAMSpath))){
        fixedpath <- getAbsolutePath(JAMSpath)
    } else {
        fixedpath <- JAMSpath
    }
    #Chomp a "/" from the end of paths
    fixedpath <- gsub("/$", "", fixedpath)

    return(fixedpath)
}


#' name_samples(list.data = NULL)
#'
#' Given a list.data object, returns a vector of sample names present in the object.
#' @export

name_samples <- function(list.data = NULL){
    loadedsamples <- gsub("_abundances", "", (names(list.data)[grep("_abundances", names(list.data))]))

    return(loadedsamples)
}


#' sample_by_category
#'
#' Sample by category
#' @export

sample_by_category <- function(rows, frac = 0.25) {
    #counts <- table(rows$category) %>% .[.>= n]
    counts <- table(rows$category)
    result <- data.frame()
    for (name in names(counts)) {
        result <- rbind(result, sample_frac(rows[rows$category==name,], frac))
    }

    return(result)
}


#' add_shape_to_plot_safely(p = NULL, shapevec = NULL, shapeby = NULL, cdict = NULL, use_letters_as_shapes = FALSE)
#'
#' #Given a vector of classes for adding shapes to a ggplot, attributes shapes safely in the presence or absence of a cdict containing shape info
#' @export

add_shape_to_plot_safely <- function (p = NULL, shapevec = NULL, shapeby = NULL, cdict = NULL, use_letters_as_shapes = FALSE){

    numshapes <- length(unique(shapevec))
    if (!use_letters_as_shapes & (numshapes <= 27)){
        shape_pecking_order <- c(19, 17, 15, 8, 12, 13, 18, 10, 3, 4, 11, 0, 1, 2, 5, 6, 7, 6, 35, 36, 38, 64, 163, 165, 167, 169, 198)
    } else {
        #Use letters
        shape_pecking_order <- c(65:90, 97:122, 35:38, 134:140)
    }

    if (!is.null(cdict)){
        st <- cdict[[shapeby]]
        if ("Shape" %in% colnames(st)){
            groupshapes <- setNames(as.numeric(st$Shape), as.character(st$Name))
            p <- p + scale_shape_manual(values = groupshapes)
        } else {
            p <- p + scale_shape_manual(values = shape_pecking_order[1:numshapes])
        }
    } else {
        p <- p + scale_shape_manual(values = shape_pecking_order[1:numshapes])
    }

    return(p)
}


#' find_clusters_in_matrix(input_matrix = NULL, kmeans_split = NULL)
#'
#' #Splits a matrix with counts in to k clusters by k-means clustering
#' @export

find_clusters_in_matrix <- function(input_matrix = NULL, kmeans_split = NULL){

    kmclusters <- stats::kmeans(input_matrix, centers = kmeans_split)
    clusterdf <- data.frame(Feature = names(fitted(kmclusters, method = "classes")), Cluster_Num = fitted(kmclusters, method = "classes"), stringsAsFactors = FALSE)
    clusterdf <- clusterdf[order(clusterdf$Cluster_Num), ]
    clusterdf$Cluster_Name <- paste("Cluster", as.character(sprintf("%02d", clusterdf$Cluster_Num)), sep = "_")

    return(clusterdf)

}


#' Author: Phillip Burger
#' Date: Sept 3, 2014
#' Purpose: Get the ordinal rank of a number.
#' http://www.phillipburger.net/wordpress/ordinal-number-suffix-function-in-r
#' @export

getOrdinalNumber1 <- function(num) {
    result <- ""
    if (!(num %% 100 %in% c(11, 12, 13))) {
        result <- switch(as.character(num %% 10),
            "1" = {paste0(num, "st")},
            "2" = {paste0(num, "nd")},
            "3" = {paste0(num, "rd")},
            paste0(num, "th"))
    } else {
        result <- paste0(num, "th")
    }
    result
}


#' install_biocon_deps(include_optional = TRUE)
#'
#' Installs Bioconductor JAMS dependencies
#' @export

install_biocon_deps <- function(include_optional = TRUE){
    #Required by JAMS
    bioconductordeps <- c("SummarizedExperiment", "ComplexHeatmap", "Biostrings")
    #Optional: dada2/DECIPHER for JAMS16, ggtree/ggtreeExtra for plot_taxsplit_tree()
    if (include_optional){
        bioconductordeps <- c(bioconductordeps, "dada2", "DECIPHER", "ggtree", "ggtreeExtra")
    }
    if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    #Do not pin a Bioconductor version: BiocManager picks the release that matches the running R.
    BiocManager::install(bioconductordeps, update = FALSE, ask = FALSE)
}


#' make_phenolabels_from_phenotable(phenotable = NULL, class_to_ignore = "N_A")
#'
#' Does exactly what the name says.
#' @export

make_phenolabels_from_phenotable <- function (phenotable = NULL, class_to_ignore = "N_A"){
    Var_label <- colnames(phenotable)
    #Delegate to the exported infer_column_type() rather than keeping a local copy of it.
    #The former local copy referred to a class_to_ignore which was never in scope here.
    Var_type <- sapply(Var_label, function (x) { infer_column_type(phenotable = phenotable, colm = x, class_to_ignore = class_to_ignore) } )
    phenolabels <- data.frame(Var_label = unname(Var_label), Var_type = unname(Var_type), stringsAsFactors = FALSE)

    return(phenolabels)

}

#' generate_filename(title = "ABC_look_at_me", add_date = TRUE, suffix = "pdf")
#'
#' Wrapper for generating file names.
#' @export

generate_filename <- function(title = "ABC_look_at_me", add_date = TRUE, suffix = "pdf"){

    if (add_date){
        dte <- Sys.Date()
        stem <- paste(title, dte, sep = "_")
    } else {
        stem <- title
    }

    fn <- paste(stem, suffix, sep = ".")

    return(fn)
}

#' analyze_fastq_filename(fn = NULL, maintain_Illumina_format = FALSE)
#'
#' Given a fastq filename, will return a list with facts about its anatomy.
#' @export

analyze_fastq_filename <- function(fn = NULL, maintain_Illumina_format = FALSE){

    filefacts <- list()
    fnsplit <- unlist(strsplit(fn, split = "\\."))
    #find out if gzipped and what the filename body is
    if (tail(fnsplit, n = 1) == "gz"){
        filefacts$gunzipped <- TRUE
        fnbody <- paste0(rev(rev(fnsplit)[3:length(fnsplit)]), collapse = ".")
        filefacts$suffix <- paste0(tail(fnsplit, n = 2), collapse = ".")
    } else {
        filefacts$gunzipped <- FALSE
        fnbody <- paste0(rev(rev(fnsplit)[2:length(fnsplit)]), collapse = ".")
        filefacts$suffix <- tail(fnsplit, n = 1)
    }

    #Now delve into fnbody and see what it it about
    #Does it have lane info? Must have L00X in the last 4 elements.
    contains_lane_info <- any(c("L001", "L002", "L003", "L004", "L005", "L006", "L007", "L008") %in% tail(unlist(strsplit(fnbody, split = "_")), n = 4))
    #Does it have paired read info? Must have RX in the last 4 elements.
    contains_paired_read_info <- any(c("R1", "R2") %in% tail(unlist(strsplit(fnbody, split = "_")), n = 4))
    #Is it Illumina output filename style?
    is_Illumina_output <- (tail(unlist(strsplit(fnbody, split = "_")), n = 1) == "001")
    #Get number of body segments
    numfnbodysegments <- length(unlist(strsplit(fnbody, split = "_")))

    #Determine the orgininal sample prefix
    #Case 1, not Illumina style
    if (!is_Illumina_output){
        if (contains_paired_read_info){
            filefacts$OriPrefix <- paste0(unlist(strsplit(fnbody, split = "_"))[1:(numfnbodysegments - 1)], collapse = "_")
            filefacts$Read <- tail(unlist(strsplit(fnbody, split = "_")), n = 1)
            filefacts$Appendage <-  paste(filefacts$Read, filefacts$suffix, sep = ".")
        } else {
            filefacts$OriPrefix <- fnbody
            filefacts$Read <- ""
            filefacts$Appendage <- filefacts$suffix
        }

    } else {
        #OK, is Illumina style prefix
        #Das Read ist immer an zweiter position
        filefacts$Read <- tail(unlist(strsplit(fnbody, split = "_")), n = 2)[1]

        if (contains_lane_info){
            #Case 2, is MiSeq, HiSeq or NextSeq file
            filefacts$OriPrefix <- paste0(unlist(strsplit(fnbody, split = "_"))[1:(numfnbodysegments - 4)], collapse = "_")
            filefacts$Lane <- tail(unlist(strsplit(fnbody, split = "_")), n = 3)[1]
            filefacts$SampNum <- tail(unlist(strsplit(fnbody, split = "_")), n = 4)[1]
            if (maintain_Illumina_format){
                filefacts$Appendage <- paste(paste(filefacts$SampNum, filefacts$Lane, filefacts$Read, "001", sep = "_"), filefacts$Container, sep = ".")
            } else {
                filefacts$Appendage <- paste(filefacts$Read, filefacts$suffix, sep = ".")
            }
        } else {
            #Case 3, is NovaSeq joined file, or contains no lane info
            filefacts$OriPrefix <- paste0(unlist(strsplit(fnbody, split = "_"))[1:(numfnbodysegments - 3)], collapse = "_")
            filefacts$Lane <- ""
            filefacts$SampNum <- tail(unlist(strsplit(fnbody, split = "_")), n = 3)[1]
            if (maintain_Illumina_format){
                filefacts$Appendage <- paste(paste(filefacts$SampNum, filefacts$Read, "001", sep = "_"), filefacts$suffix, sep = ".")
            } else {
                filefacts$Appendage <- paste(filefacts$Read, filefacts$suffix, sep = ".")
            }
        }
    }

    return(filefacts)

}


#' infer_column_type(phenotable = NULL, colm = NULL, class_to_ignore = "N_A")
#'
#' From a phenotable, infers what kind of a column it is, discrete or continuous.
#' @export

infer_column_type <- function(phenotable = NULL, colm = NULL, class_to_ignore = "N_A"){
    if (colm %in% c("Sample", "sample")){
        colmtype <- "Sample"
        return(colmtype)
    }

    #`[[` rather than `[ , colm]`: the latter returns a one-column tibble rather than a vector
    #for a tibble, and is outright wrong for a data.table, where j is evaluated inside the
    #table, so a bare character variable is returned as itself instead of as the column.
    colvals <- phenotable[[colm]]

    if (inherits(colvals, c("Date", "POSIXct", "POSIXlt", "POSIXt", "difftime"))){
        #Date-like columns are numeric under the hood (days or seconds since the epoch),
        #so can_be_made_numeric() happily calls them continuous, and define_kinds_of_variables()
        #then compares them against class_to_ignore. Treat them as categorical instead, which is
        #exactly how ISO date *strings* are already handled by this function. Note that "date" is
        #deliberately NOT returned here: load_metadata_from_file() reserves that Var_type for
        #Excel serial dates, which it pushes through convertToDate().
        colmtype <- "discrete"
    } else {
        if (can_be_made_numeric(colvals, cats_to_ignore = class_to_ignore)){
            colmtype <- "continuous"
        } else {
            colmtype <- "discrete"
        }
    }

    return(colmtype)

}

#' compute_x_font(nx = 1, upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 8, cex = 0.7)
#'
#' Compute the appropriate size of font to be used in JAMS heatmap column (sample) annotations.
#' @export

compute_x_font <- function(nx = 1, upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 8, cex = 0.7){

    m <- ((lower_fs - upper_fs) / (upper_n - lower_n)) * -1
    cnst <- lower_fs - (lower_n * m)
    fontsizex <- round(((m * nx) + cnst), 2)
    fontsizex <- fontsizex * cex

    return(fontsizex)
}

#' hm_fontsize_computer(mat_rownames = NULL, mat_colnames = NULL, upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 8, cex = 0.7)
#'
#' Returns a vector of length 2 with values for the appropriate sizes of fonts to be used in JAMS heatmap columns (samples) and rows (features), respectively.
#' @export

hm_fontsize_computer <- function(mat_rownames = NULL, mat_colnames = NULL, upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 8, cex = 0.7){

    cexy <- cex
    #If not LKT and featname is split, then decrease row font size
    if (!(any(length(grep("^LKT_", mat_rownames) > 0)))){
        if (max(nchar(mat_rownames)) > 60){
            cexy <- cex * 0.7
        }
    }
    fontsizex <- compute_x_font(nx = length(mat_colnames), upper_n = upper_n, upper_fs = upper_fs, lower_n = lower_n, lower_fs = lower_fs, cex = cex)
    fontsizey <- compute_x_font(nx = length(mat_rownames), upper_n = upper_n, upper_fs = upper_fs, lower_n = lower_n, lower_fs = lower_fs, cex = cexy)
    #floor rownames to 0.1
    fontsizey <- max(0.1, fontsizey)
    xy_fs <- c(fontsizex, fontsizey)

    return(xy_fs)
}

#' split_featname(featname = NULL, thresh_featname_split = 40)
#'
#' For a feature name (or any string, for that matter) will add a carriage return to the middle of the string if it exceeds the threshold number of characters. This is used in the heatmap function.
#' @export

split_featname <- function(featname = NULL, thresh_featname_split = 40){
    if (nchar(featname) >= thresh_featname_split){
        len1 <- round((nchar(featname) / 2), 0)
        len2 <- nchar(featname) - len1
        part1 <- paste0(unlist(strsplit(featname, split = ""))[1:len1], collapse = "")
        part2 <- paste0(unlist(strsplit(featname, split = ""))[(1:len2) + len1], collapse = "")
        featname <- paste0(part1, "-\n", part2)
    }

    return(featname)
}


#' post_message(msg = NULL, face = "plain", just = "centre", color = "black", size = 20)
#'
#' Wrapper function for plotting text to a full page on a pdf.
#' @export

post_message <- function(msg = NULL, face = "plain", just = "centre", color = "black", size = 20){
    text_msg <- paste0(msg, collapse = "\n")
    tgrob <- text_grob(text_msg, face = face, color = color, just = just, size = size)
    print(as_ggplot(tgrob))
}


#' print_jams_ASCII(palette = "neon")
#'
#' Prints to terminal a retro ASCII JAMS logo.
#' @export

print_jams_ASCII <- function (palette = NULL){

    palettes <- list(
        #apple2 = c(34, 11, 202, 196, 127, 20),
        apple2rev = rev(c(34, 11, 202, 196, 127, 20)),
        neon = c(27, 93, 129, 201, 207, 99),
        sunset = c(202, 208, 214, 220, 226, 228),
        ocean = c(25, 31, 38, 44, 50, 56),
        fire = c(196, 202, 208, 214, 220, 226),
        grad_dark = c(16, 17, 18, 19, 20, 21),
        candy = c(21, 57, 93, 129, 165, 201),
        reds = 52:57,
        reds2 = (36 * 0:5) + 17,
        skies = (36 * 0:5) + 37,
        green_blue = 22:27,
        greys = 244:249,
        dark2 = c(16, 17, 18, 54, 53, 52),
        mono1 = rep(5, 6),
        mono2 = rep(4, 6),
        mono3 = rep(124, 6),
        mono4 = rep(246, 6)
     )

    if (is.null(palette)){
        #choose random palette
        palette <- sample(names(palettes), 1, replace = TRUE)
        rowvec <- palettes[[palette]]
    } else {
        rowvec <- palettes[[palette]]
    }

    ASCII_lines <- c(sprintf("\033[38;5;%dm         ██╗\033[38;5;%dm     █████╗\033[38;5;%dm     ███╗   ███╗\033[38;5;%dm     ███████╗\033[0m", rowvec[1], rowvec[1], rowvec[1], rowvec[1]), sprintf("\033[38;5;%dm        ██║\033[38;5;%dm    ██╔══██╗\033[38;5;%dm    ████╗ ████║\033[38;5;%dm    ██╔════╝\033[0m", rowvec[2], rowvec[2], rowvec[2], rowvec[2]), sprintf("\033[38;5;%dm       ██║\033[38;5;%dm    ███████║\033[38;5;%dm    ██╔████╔██║\033[38;5;%dm    ███████╗\033[0m", rowvec[3], rowvec[3], rowvec[3], rowvec[3]), sprintf("\033[38;5;%dm  ██   ██║\033[38;5;%dm    ██╔══██║\033[38;5;%dm    ██║╚██╔╝██║\033[38;5;%dm    ╚════██║\033[0m", rowvec[4], rowvec[4], rowvec[4], rowvec[4]),sprintf("\033[38;5;%dm ╚█████╔╝\033[38;5;%dm    ██║  ██║\033[38;5;%dm    ██║ ╚═╝ ██║\033[38;5;%dm    ███████║\033[0m", rowvec[5], rowvec[5], rowvec[5], rowvec[5]), sprintf("\033[38;5;%dm  ╚════╝ \033[38;5;%dm   ╚═╝  ╚═╝\033[38;5;%dm    ╚═╝     ╚═╝\033[38;5;%dm    ╚══════╝\033[0m", rowvec[6], rowvec[6], rowvec[6], rowvec[6]), paste("Microbial genomics analysis software.", paste("Version", packageVersion("JAMS"))))

    cat(paste0(ASCII_lines, collapse = "\n"), "\n")
}

#' extract_NCBI_taxid_from_featname(Taxon = NULL, NCBI_taxonomic_rank = NULL)
#'
#' Securely extracts the NCBI taxid from a JAMS2-style taxonomic feature name
#' @export

extract_NCBI_taxid_from_featname <- function(Taxon = NULL){

    #MetaBAT bin, so determine whether it is one which contains sample name (JAMSbeta style) or not (JAMSalpha style).
    #If LKT or cLKT, this will also work
    split_elements <- unlist(strsplit(Taxon, split = "_"))
    #Get earliest position matching a taxonomic tag
    taxtag_pos <- which(split_elements %in% c("d", "k", "p", "c", "o", "f", "g", "s", "is1"))[1]
    NCBI_taxid <- split_elements[(taxtag_pos + 2)]

    #If not found, then there is no tag, so convert NA to taxid 0 (root)
    if (is.na(NCBI_taxid)){
        NCBI_taxid <- "0"
    }
    #If is not numeric, like "missing", or "Unclassified", return "0"
    if (!can_be_made_numeric(NCBI_taxid)){
         NCBI_taxid <- "0"
    }

    return(NCBI_taxid)
}


#' retrieve_GCA_accession_URL(assemblytargetname = NULL)
#'
#' Given a single valid GCA assembly accession number, this function will return the https URL for downloading the fasta.gz file for that assembly accession.
#' @export

retrieve_GCA_accession_URL <- function(assemblytargetname){

    #First, check you have esearch
    missingdep <- FALSE
    for (dep in c("esearch", "esummary", "xtract")) {
        if (system(str_c("which ", dep), ignore.stdout = TRUE) ==1 ) { 
            flog.warn(str_c("You are missing ", dep))
            missingdep <- TRUE
        }
    }

    if (!missingdep){
        esearchcmd <- paste("esearch -db assembly -query", assemblytargetname, "| esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank >", paste(assemblytargetname, "url", sep = "."))
        system(esearchcmd)
        assemblyURL <- readLines(paste(assemblytargetname, "url", sep = "."), n = 1L, warn = FALSE)
        assemblyURL <- sub("^ftp:", "https:", assemblyURL)
        assemblyURL <- paste(paste(assemblyURL, basename(assemblyURL), sep = "/"), "genomic.fna.gz", sep = "_")

        #if GCA accession is spoof, return NULL
        if (assemblyURL == "_genomic.fna.gz"){
            assemblyURL <- NULL
        }

        return(assemblyURL)

    } else {
        #You don't have what it takes to obtain the accession URL
        return(NULL)
    }
}

#' pretty_time(seconds)
#'
#' Automatically reformats numbers of seconds to a more human-readable format.
#' @export

pretty_time <- function(seconds) {
    if (seconds < 60) {
        sprintf("%.1f s", seconds)
    } else if (seconds < 3600) {
        sprintf("%.1f min", seconds / 60)
    } else {
        sprintf("%.1f h", seconds / 3600)
    }
}


#' tally_metadata(md = NULL, column_to_group_by = NULL, columns_to_tally = NULL)
#'
#' Tallies the number of classes present in a metadata df and returns a data frame with counts.
#' @export

tally_metadata <- function(md = NULL, column_to_group_by = NULL, columns_to_tally = NULL){

    tally_column <- function(md = NULL, column_to_group_by = NULL, col2tally = NULL){ 
        tmp_tally_df <- md %>% dplyr::count(across(all_of(column_to_group_by)), across(all_of(col2tally))) %>% tidyr::pivot_wider(names_from = all_of(col2tally), values_from = n, values_fill = 0) %>% dplyr::rename_with( ~ paste0("Num_Samples_in_", col2tally, ".", .x), -all_of(column_to_group_by)) %>% as.data.frame()

        return(tmp_tally_df)
    }
 
    tallied_column_list <- lapply(columns_to_tally, function (x) { tally_column(md = md, column_to_group_by = column_to_group_by, col2tally = x) })

    #Stitch all together into a single data frame
    md_tally_df <- purrr::reduce(tallied_column_list, dplyr::left_join, by = column_to_group_by)

    return(md_tally_df)
}



#' plot_table_a4_autofit_paged(dataframe = NULL, title = NULL, orientation = NULL, rows_per_page = 35, start_fontsize = 10, min_fontsize = 5, margins_portrait_mm = c(top=6, right=6, bottom=6, left=6), margins_landscape_mm = c(top=6, right=8, bottom=6, left=8), gap_mm = 2, cell_padding_mm = c(v = 0.6, h = 1.2), header_padding_mm = c(v = 0.8, h = 1.2), lineheight = 0.95, wrap = FALSE, wrap_chars = NULL, wrap_min = 12, wrap_max = 60)
#' EXPERIMENTAL novel plotting function for printing data frames into a pdf.
#' Usage:
#' pdf("my_landscape.pdf", height = 210/25.4, width = 297/25.4)  # A4 landscape
#' plot_table_a4_autofit_paged(dataframe = opt$featuredata[1:80, ], title="Featuretable", rows_per_page = 35, wrap = TRUE, wrap_chars = 20, start_fontsize = 15, min_fontsize = 2,  orientation = "landscape")
#' dev.off()
#'
#' @export

plot_table_a4_autofit_paged <- function(dataframe = NULL, title = NULL, orientation = "landscape", rows_per_page = 35, start_fontsize = 10, min_fontsize = 5, margins_portrait_mm = c(top=6, right=6, bottom=6, left=6), margins_landscape_mm = c(top=6, right=8, bottom=6, left=8), gap_mm = 2, cell_padding_mm = c(v = 0.6, h = 1.2), header_padding_mm = c(v = 0.8, h = 1.2), lineheight = 0.95, wrap = FALSE, wrap_chars = NULL, wrap_min = 12, wrap_max = 60) {

    if (!is.data.frame(dataframe)) {
        stop("`dataframe` must be a data.frame")
    }
    orientation <- match.arg(orientation)

    # define margins based on orientation
    margins_mm <- if (orientation == "portrait") margins_portrait_mm else margins_landscape_mm

    if (length(margins_mm) == 1) {
        margins_mm <- c(top=margins_mm, right=margins_mm, bottom=margins_mm, left=margins_mm)
    } else if (is.null(names(margins_mm))) {
        stop("margins_*_mm must be named (top/right/bottom/left) or a single number")
    }

    # optional wrapping (applied once to full df so it stays consistent across pages)
    df <- dataframe
  if (wrap) {
    if (is.null(wrap_chars)) {
      nc <- max(1, ncol(df))
      base <- round(90 / nc)
      bump <- round((start_fontsize - 8) * 2)
      wrap_chars <- max(wrap_min, min(wrap_max, base + bump))
    }
    wrap_one <- function(x) {
      x <- as.character(x)
      vapply(x, function(s) paste(strwrap(s, width = wrap_chars), collapse = "\n"), character(1))
    }
    df[] <- lapply(df, wrap_one)
  }

  # device size in mm (device must be open)
  dev_in <- dev.size("in")
  page_w_mm <- dev_in[1] * 25.4
  page_h_mm <- dev_in[2] * 25.4

  avail_w_mm <- page_w_mm - (margins_mm["left"] + margins_mm["right"])
  avail_h_mm <- page_h_mm - (margins_mm["top"] + margins_mm["bottom"])
  if (avail_w_mm <= 0 || avail_h_mm <= 0) stop("Margins too large for device size")

  # build title grob once; show it on every page by default
  title_g <- NULL
  title_h_mm <- 0
  if (!is.null(title) && nzchar(title)) {
    title_g <- textGrob(title, gp = gpar(fontsize = max(10, start_fontsize + 2)))
    title_h_mm <- convertHeight(grobHeight(title_g), "mm", valueOnly = TRUE)
  }

  make_theme <- function(fs) {
    ttheme_default(
      base_size = fs,
      core = list(
        padding = unit(c(cell_padding_mm["v"], cell_padding_mm["h"]), "mm"),
        fg_params = list(lineheight = lineheight)
      ),
      colhead = list(
        padding = unit(c(header_padding_mm["v"], header_padding_mm["h"]), "mm"),
        fg_params = list(lineheight = lineheight)
      )
    )
  }

  measure_tbl_mm <- function(tbl) {
    w <- convertWidth(sum(tbl$widths), "mm", valueOnly = TRUE)
    h <- convertHeight(sum(tbl$heights), "mm", valueOnly = TRUE)
    c(width = as.numeric(w), height = as.numeric(h))
  }

  # split into pages (35 rows per page if needed)
  n <- nrow(df)
  if (n == 0) {
    warning("Empty data frame: nothing to plot.")
    return(invisible(NULL))
  }

  starts <- seq(1, n, by = rows_per_page)
  ends <- pmin(starts + rows_per_page - 1, n)

  # For each page, fit fontsize (can vary per page, but we’ll keep it consistent:
  # choose the *smallest* fontsize required across all pages so pages match visually.)
  fitted_fontsizes <- integer(length(starts))

  # First pass: compute needed fontsize per chunk
  for (i in seq_along(starts)) {
    chunk <- df[starts[i]:ends[i], , drop = FALSE]
    fs <- start_fontsize
    repeat {
      tbl <- tableGrob(chunk, rows = NULL, theme = make_theme(fs))
      dims <- measure_tbl_mm(tbl)

      avail_table_h_mm <- avail_h_mm - title_h_mm - ifelse(is.null(title_g), 0, gap_mm)
      if (avail_table_h_mm <= 0) avail_table_h_mm <- avail_h_mm - title_h_mm

      fits <- (dims["width"] <= avail_w_mm) && (dims["height"] <= avail_table_h_mm)

      if (fits || fs <= min_fontsize) break
      fs <- fs - 1
    }
    fitted_fontsizes[i] <- fs
  }

  # Use a single fontsize across pages for consistency (smallest that worked)
  final_fs <- min(fitted_fontsizes)

  # Second pass: draw pages
  for (i in seq_along(starts)) {
    chunk <- df[starts[i]:ends[i], , drop = FALSE]
    tbl <- tableGrob(chunk, rows = NULL, theme = make_theme(final_fs))

    grid.newpage()

    # margins in npc
    left_npc   <- margins_mm["left"]   / page_w_mm
    right_npc  <- margins_mm["right"]  / page_w_mm
    top_npc    <- margins_mm["top"]    / page_h_mm
    bottom_npc <- margins_mm["bottom"] / page_h_mm

    avail_w_npc <- 1 - left_npc - right_npc

    # title
    title_h_npc <- 0
    gap_npc <- 0
    if (!is.null(title_g)) {
      title_h_npc <- convertHeight(grobHeight(title_g), "npc", valueOnly = TRUE)
      gap_npc <- gap_mm / page_h_mm

      pushViewport(viewport(
        x = left_npc + avail_w_npc / 2,
        y = 1 - top_npc - title_h_npc / 2,
        width = avail_w_npc,
        height = title_h_npc,
        just = c("center", "center")
      ))
      grid.draw(title_g)
      popViewport()
    }

    # table area below title
    table_top <- 1 - top_npc - title_h_npc - gap_npc
    table_bottom <- bottom_npc
    table_h <- table_top - table_bottom
    table_center_y <- table_bottom + table_h / 2

    pushViewport(viewport(
      x = left_npc + avail_w_npc / 2,
      y = table_center_y,
      width = avail_w_npc,
      height = table_h,
      just = c("center", "center"),
      clip = "on"
    ))
    grid.draw(tbl)
    popViewport()
  }

  invisible(list(
    pages = length(starts),
    rows_per_page = rows_per_page,
    final_fontsize = final_fs,
    wrap_chars = if (wrap) wrap_chars else NA
  ))
}

#' resolve_dichotomization_threshold(x = NULL, threshold = NULL, func_for_threshold = "median")
#'
#' Internal helper shared by dichotomize_variable_in_md() and plot_Ordination(). Returns the numeric cut-point used to split a continuous vector into two groups, plus a short human-readable description of how it was obtained (for plot titles and log messages). If threshold is given, it is validated and used as is; otherwise func_for_threshold (a function, or the name of one, e.g. "median", "mean") is applied to the non-missing values of x. The function is called with the missing values already removed and WITHOUT an na.rm argument, so user functions that do not accept na.rm, such as function(v) quantile(v, 0.75), work.
#' @keywords internal
#' @noRd

resolve_dichotomization_threshold <- function(x = NULL, threshold = NULL, func_for_threshold = "median"){

    x_ok <- x[!is.na(x)]
    if (length(x_ok) == 0){
        stop("Unable to dichotomize: there are no usable (non-missing, numeric) values in the variable.")
    }

    fmt_num <- function(v){ format(signif(v, 6), scientific = FALSE, trim = TRUE) }

    if (!is.null(threshold)){
        if (!(is.numeric(threshold) && (length(threshold) == 1) && is.finite(threshold))){
            stop("The dichotomization threshold must be a single finite number.")
        }
        thr <- unname(threshold)
        thrdesc <- paste0("threshold ", fmt_num(thr))
    } else {
        if (is.function(func_for_threshold)){
            thrfun <- func_for_threshold
            funcname <- "custom function"
        } else if (is.character(func_for_threshold) && (length(func_for_threshold) == 1)){
            thrfun <- tryCatch(match.fun(func_for_threshold), error = function(e){
                stop(paste0("Unable to find a function called \"", func_for_threshold, "\" for computing the dichotomization threshold."), call. = FALSE)
            })
            funcname <- func_for_threshold
        } else {
            stop("func_for_threshold must be a function, or the name of a function (e.g. \"median\"), returning a single number.")
        }

        thr <- tryCatch(thrfun(x_ok), error = function(e){
            stop(paste0("Applying ", funcname, " to compute the dichotomization threshold failed: ", conditionMessage(e)), call. = FALSE)
        })
        if (!(is.numeric(thr) && (length(thr) == 1) && is.finite(thr))){
            stop(paste0("The function used to compute the dichotomization threshold (", funcname, ") must return a single finite number, but it returned something else (length ", length(thr), ", class ", class(thr)[1], ")."))
        }
        thr <- unname(thr)
        thrdesc <- paste0(funcname, " (", fmt_num(thr), ")")
    }

    return(list(threshold = thr, description = thrdesc, formatted_threshold = fmt_num(thr)))
}


#' dichotomize_variable_in_md(md = NULL, column_to_dichotomize = NULL, threshold = NULL, func_for_threshold = "median", class_to_ignore = "N_A", high_label = "High", low_label = "Low", missing_label = NULL)
#'
#' Translates a continuous variable in metadata into a dichotomous discrete variable, either by a set threshold, or by being above or below a particular value computed by func_for_threshold, e.g. the median. Values greater than or equal to the threshold are labelled high_label ("High"), values strictly below it low_label ("Low"). Ties with the threshold therefore always go to the High group, which matters for the median of variables with many repeated values.
#'
#' In keeping with the JAMS convention that all metadata cells are character, numbers stored as strings (e.g. "52.5") are handled transparently, and factors are converted through their labels rather than their internal level codes. Cells that are NA or that belong to class_to_ignore are treated as missing, do not contribute to the threshold, and are returned as missing_label. Any other cell that cannot be read as a number is also returned as missing_label, with a warning naming the offending values.
#'
#' @param md data.frame (or S4 DataFrame, e.g. colData of a SummarizedExperiment) containing the variable.
#' @param column_to_dichotomize String giving the name of the column to dichotomize.
#' @param threshold Optional single number to use as the cut-point. If given, func_for_threshold is ignored. Default is NULL.
#' @param func_for_threshold Function, or the name of a function, applied to the non-missing values to obtain the cut-point when threshold is NULL. It is called without an na.rm argument (missing values are removed beforehand) and must return a single number. Default is "median".
#' @param class_to_ignore String or vector of strings to be treated as missing data. Default is "N_A".
#' @param high_label String used for values at or above the threshold. Default is "High".
#' @param low_label String used for values below the threshold. Default is "Low".
#' @param missing_label Value returned for missing or unreadable cells. If NULL (the default), the first element of class_to_ignore is used (i.e. "N_A"), or NA if class_to_ignore is NULL.
#' @return Unnamed character vector, the same length and order as nrow(md).
#' @export

dichotomize_variable_in_md <- function(md = NULL, column_to_dichotomize = NULL, threshold = NULL, func_for_threshold = "median", class_to_ignore = "N_A", high_label = "High", low_label = "Low", missing_label = NULL){

    if (!(is.data.frame(md) || inherits(md, c("DataFrame", "DFrame")))){
        stop("md must be a data.frame (or an S4 DataFrame such as the colData of a SummarizedExperiment).")
    }
    if (!(is.character(column_to_dichotomize) && (length(column_to_dichotomize) == 1))){
        stop("column_to_dichotomize must be a single column name.")
    }
    if (!(column_to_dichotomize %in% colnames(md))){
        stop(paste0("Column \"", column_to_dichotomize, "\" was not found in md."))
    }
    if (!(is.character(high_label) && is.character(low_label) && (length(high_label) == 1) && (length(low_label) == 1) && (high_label != low_label))){
        stop("high_label and low_label must be two different single strings.")
    }
    if (is.null(missing_label)){
        missing_label <- if (is.null(class_to_ignore)) NA_character_ else as.character(class_to_ignore[1])
    }

    #`[[` works for data.frame, tibble, data.table and S4 DataFrame alike, and always returns a plain vector.
    rawvals <- md[[column_to_dichotomize]]

    if (is.numeric(rawvals) || is.logical(rawvals) || inherits(rawvals, c("Date", "POSIXt", "difftime"))){
        #Already numeric-like; use directly, avoiding a lossy round trip through character.
        xnum <- as.numeric(rawvals)
        #Honour class_to_ignore here too (e.g. a numeric sentinel such as -999)
        if (!is.null(class_to_ignore)){
            xnum[as.character(rawvals) %in% as.character(class_to_ignore)] <- NA
        }
    } else {
        #Character (the JAMS norm) or factor. as.character() first, so a factor yields its labels, not its level codes.
        rawchar <- as.character(rawvals)
        ignored <- is.na(rawchar) | (rawchar %in% as.character(class_to_ignore))
        xnum <- suppressWarnings(as.numeric(rawchar))
        xnum[ignored] <- NA
        unreadable <- (!ignored) & is.na(xnum)
        if (any(unreadable)){
            badvals <- unique(rawchar[unreadable])
            flog.warn(paste0("dichotomize_variable_in_md: ", sum(unreadable), " value(s) of \"", column_to_dichotomize, "\" could not be read as numbers and will be set to ", missing_label, ": ", paste0(utils::head(badvals, 10), collapse = ", "), if (length(badvals) > 10) ", ..." else "", ". If these represent missing data, consider adding them to class_to_ignore."))
        }
    }

    thrinfo <- resolve_dichotomization_threshold(x = xnum, threshold = threshold, func_for_threshold = func_for_threshold)

    dichot_vector <- ifelse(xnum >= thrinfo$threshold, high_label, low_label)
    dichot_vector[is.na(xnum)] <- missing_label

    return(unname(as.character(dichot_vector)))
}


#' Calculate safe number of threads to avoid Out-Of-Memory errors when using contextualize_taxonomy
#'
#' @param requested_threads The maximum number of CPUs requested by the user.
#' @param totmembytes The total available system memory in bytes.
#' @param list.data The list.data object that will be exported to workers.
#' @param BinsDF The BinsDF object that will be exported to workers.
#' @return An integer representing the safe number of threads to use.
#' @export

calculate_safe_threads <- function(requested_threads= NULL, totmembytes = NULL, list.data = NULL, BinsDF = NULL) {
    
    #If totmembytes is missing or invalid, fallback to requested threads
    if (is.null(totmembytes) || is.na(totmembytes) || totmembytes <= 0) {
        return(requested_threads)
    }

    size_list_data <- as.numeric(utils::object.size(list.data))
    size_BinsDF <- as.numeric(utils::object.size(BinsDF))

    #overhead for each worker
    r_overhead_bytes <- 400 * 1024^2 

    #Calculate memory per worker, with 25% margin
    mem_per_worker <- (size_list_data + size_BinsDF + r_overhead_bytes) * 1.25
    available_mem_for_workers <- totmembytes - mem_per_worker

    if (available_mem_for_workers <= 0) {
        safe_threads <- 1
    } else {
        safe_threads <- floor(available_mem_for_workers / mem_per_worker)
    }

    final_threads <- min(requested_threads, safe_threads)
    final_threads <- max(1, final_threads)

    return(final_threads)
}


#' resolve_tax_annotation_colours(feature_table = NULL, row_order = NULL, want_phylum = FALSE, want_gram = FALSE)
#'
#' Internal helper shared by JAMS heatmap/correlation plotting functions to resolve
#' per-feature Phylum and Gram annotation colours consistently. Phylum colours are
#' sourced from the package PhyCols dictionary (data(PhyCols)), with collision-aware
#' fallback colours assigned to any phylum not in the dictionary. Gram is read directly
#' from the feature table as a per-feature property (a phylum may contain both Gram
#' positive and Gram negative species), and mapped to a fixed, total colour scheme.
#'
#' Returns a list with: phylum (character vector, ordered to row_order), gram (character
#' vector, ordered to row_order), phylum_cols (named colour vector) and gram_cols (named
#' colour vector). Either annotation can be requested independently via want_phylum /
#' want_gram; unrequested components are returned as NULL. If a requested component cannot
#' be resolved (missing column), it degrades gracefully to NULL with a warning rather than
#' erroring, so callers can simply skip that annotation.
#'
#' Do not attempt to use this outside of the JAMS plotting functions.
#' @export

resolve_tax_annotation_colours <- function(feature_table = NULL, row_order = NULL, want_phylum = FALSE, want_gram = FALSE){

    out <- list(phylum = NULL, gram = NULL, phylum_cols = NULL, gram_cols = NULL)

    if (!any(c(want_phylum, want_gram))){
        return(out)
    }

    tt <- as.data.frame(feature_table)

    #Order to the matrix rows if supplied; otherwise leave as-is.
    if (!is.null(row_order)){
        row_order <- row_order[row_order %in% rownames(tt)]
        tt <- tt[row_order, , drop = FALSE]
    }

    #--- Phylum colours ---
    if (want_phylum){
        if ("Phylum" %in% colnames(tt)){
            data(PhyCols)
            phyla_vec <- as.character(tt$Phylum)
            phyla_vec[is.na(phyla_vec)] <- "p__Missing"
            phyla_present <- unique(phyla_vec)

            phycols <- PhyCols[names(PhyCols) %in% phyla_present]
            missing_phyla <- setdiff(phyla_present, names(phycols))
            if (length(missing_phyla) > 0){
                #Fallback colours that steer clear of the fixed PhyCols hues already in use.
                used_cols <- toupper(unname(phycols))
                fallback_pool <- toupper(grDevices::rainbow(length(missing_phyla) + length(used_cols) + 1))
                fallback_pool <- fallback_pool[!(fallback_pool %in% used_cols)]
                fallback_cols <- fallback_pool[seq_along(missing_phyla)]
                if (any(is.na(fallback_cols))){
                    fallback_cols[is.na(fallback_cols)] <- "#BCC2C2"
                }
                names(fallback_cols) <- missing_phyla
                phycols <- c(phycols, fallback_cols)
            }
            #Conventional fixed colours for the catch-all categories, if present.
            if ("p__Unclassified" %in% names(phycols)) phycols["p__Unclassified"] <- "#000000"
            if ("p__Missing" %in% names(phycols))      phycols["p__Missing"]      <- "#757575"

            out$phylum <- phyla_vec
            out$phylum_cols <- phycols
        } else {
            flog.warn("resolve_tax_annotation_colours: Phylum annotation requested but no Phylum column is present in the feature table (expected when agglomerating to Phylum or above). Skipping Phylum annotation.")
        }
    }

    #--- Gram colours ---
    if (want_gram){
        if ("Gram" %in% colnames(tt)){
            #Preserve all four canonical Gram categories from JAMStaxtable:
            #positive, negative, N_A (does not apply, e.g. eukaryotes/viruses) and
            #undetermined (unknown). Anything unexpected (including NA) folds into
            #"undetermined", the honest "we don't know" bucket.
            canonical_gram <- c("positive", "negative", "N_A", "undetermined")
            gram_vec <- as.character(tt$Gram)
            gram_vec[is.na(gram_vec)] <- "undetermined"
            gram_vec[!(gram_vec %in% canonical_gram)] <- "undetermined"

            gramcols_all <- c(positive = "#7D00C4", negative = "#FC0345", "N_A" = "#000000", undetermined = "#BCC2C2")
            #Keep only categories actually present, in canonical order.
            gramcols <- gramcols_all[canonical_gram[canonical_gram %in% unique(gram_vec)]]

            out$gram <- gram_vec
            out$gram_cols <- gramcols
        } else {
            flog.warn("resolve_tax_annotation_colours: Gram annotation requested but no Gram column is present in the feature table. Skipping Gram annotation.")
        }
    }

    return(out)
}


#' fix_Gram_in_SEobj(SEobj = NULL)
#' Use this function to safely correct (or add) Gram to the featuretable of a JAMS2 compatible SummarizedExperimetn object based off of what information is available on the JAMStaxtable of the version of JAMS you have loaded (please always used the latest).
#'
#' @param SEobj A JAMS2-compatible SummarizedExperiment object.
#' @export

fix_Gram_in_SEobj <- function(SEobj = NULL){
    tt <- as.data.frame(rowData(SEobj))
    if (all(c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT") %in% colnames(tt))){
        flog.info("Fixing Gram column in SEobj rowData")
        data(JAMStaxtable)
        tt_cols <- colnames(tt)
        if (!("Gram" %in% tt_cols)){
            flog.warn("A GRAM column was not found on the original featuretable of the SEobj. Adding one now.")
            tt_cols <- c("Gram", tt_cols)
        }
        #Delete spurious Gram information
        tt$Gram <- NULL
        #Add Taxid
        tt$Taxid <- sapply(tt$LKT, function(x) { extract_NCBI_taxid_from_featname(Taxon = x) } )
        #Fill in dark matter Taxids with "0"
        tt$Taxid[which(is.na(tt$Taxid))] <- "0"
        #Add Gram information
        Taxid2gram <- JAMStaxtable[which(JAMStaxtable$Taxid %in% tt$Taxid), c("Taxid", "Gram")]
        Taxid2gram <- Taxid2gram[!(duplicated(Taxid2gram$Taxid)), ]
        tt <- left_join(as.data.frame(tt), Taxid2gram, by = "Taxid")
        tt[which(is.na(tt[ , "Gram"])), "Gram"] <- "undetermined"
        tt <- tt[ , tt_cols]
        rownames(tt) <- tt[ , "LKT"]
        rowData(SEobj) <- tt
    } else {
        flog.warn(paste("Your SEobj of input must be a taxonomic object. Featuretable with all taxonomy levels was not found."))
    }

    return(SEobj)
}

#' basename_safe(path)
#' basename() that returns NA rather than erroring on NULL/NA/empty input.
#' @export
basename_safe <- function(path){
    if (is.null(path) || length(path) == 0 || all(is.na(path)) || all(path == "")) return(NA_character_)
    basename(path)
}

#' export_df(df = NULL, output_basename = NULL, export_path = NULL, export_format = "xlsx", write_rownames = FALSE, add_date = TRUE)
#'
#' Quick wrapper for exporting a data frame or a (named) list of data frames to the system
#' in JAMS-preferred formatting, without having to remember the write.xlsx / write.table argument
#' differences each time. For "xlsx", a list is written as one workbook with one tab per element.
#' For "tsv", a list is written serially to one tsv file per element.
#' @export

export_df <- function(df = NULL, output_basename = NULL, export_path = NULL, export_format = "xlsx", write_rownames = FALSE, add_date = TRUE){

    #output_basename is required
    if (is.null(output_basename)){
        flog.warn("You must supply an output_basename. This is a required argument. Nothing was exported.")
        return(invisible(NULL))
    }

    #Vet format
    if (!(export_format %in% c("xlsx", "tsv"))){
        flog.warn("export_format must be either \"xlsx\" or \"tsv\". Nothing was exported.")
        return(invisible(NULL))
    }

    #Check there is something to export
    if (is.null(df)){
        flog.warn("You must supply a data frame or a list of data frames to df. Nothing was exported.")
        return(invisible(NULL))
    }

    #Standardise input into a named list of data frames
    if (is.data.frame(df)){
        df_list <- list(df)
        names(df_list) <- output_basename
    } else if (is.list(df)){
        #Ensure every element is a data frame
        if (!all(sapply(df, is.data.frame))){
            flog.warn("When df is a list, every element must be a data frame. Nothing was exported.")
            return(invisible(NULL))
        }
        df_list <- df
        #Make sure elements are named, as these become xlsx tabs or tsv filename parts
        if (is.null(names(df_list)) || any(names(df_list) == "")){
            flog.info("Some list elements are unnamed. Naming them sequentially.")
            blanks <- which(names(df_list) == "" | is.na(names(df_list)))
            if (is.null(names(df_list))){
                names(df_list) <- paste0("Sheet_", seq_along(df_list))
            } else {
                names(df_list)[blanks] <- paste0("Sheet_", blanks)
            }
        }
    } else {
        flog.warn("df must be either a data frame or a list of data frames. Nothing was exported.")
        return(invisible(NULL))
    }

    #Resolve and vet export path
    if (is.null(export_path)){
        export_path <- getwd()
    } else {
        export_path <- fixrelpath(export_path)
        if (!dir.exists(export_path)){
            flog.warn(paste("The export_path supplied does not exist:", export_path, "- defaulting to current working directory instead."))
            export_path <- getwd()
        } else if (file.access(export_path, mode = 2) != 0){
            flog.warn(paste("HEADS UP: the export_path supplied is not writeable:", export_path, "- defaulting to current working directory instead."))
            export_path <- getwd()
        } else {
            flog.warn(paste("HEADS UP: output is being written to", export_path, "and NOT to the current working directory."))
        }
    }

    if (export_format == "xlsx"){
        #Cap tab names to 31 characters because of Microsoft Excel.
        if (any(nchar(names(df_list)) > 31)){
            flog.info("Capping tab names to 31 characters because of Microsoft Excel.")
            newnames <- sapply(names(df_list), function(x){ stringr::str_trunc(x, 29) })
            newnames <- make.unique(newnames)
            names(df_list) <- newnames
        }
        fn <- generate_filename(title = output_basename, add_date = add_date, suffix = "xlsx")
        fullfn <- file.path(export_path, fn)
        flog.info(paste("Saving spreadsheet as", fullfn))
        write.xlsx(df_list, file = fullfn, asTable = TRUE, rowNames = write_rownames, colNames = TRUE, colWidths = "auto", borders = "all")
    } else {
        #tsv: one file per data frame
        for (nm in names(df_list)){
            #Fold the element name into the basename so serial files are distinguishable
            curr_title <- if (length(df_list) > 1) paste(output_basename, nm, sep = "_") else output_basename
            fn <- generate_filename(title = curr_title, add_date = add_date, suffix = "tsv")
            fullfn <- file.path(export_path, fn)
            flog.info(paste("Saving tsv as", fullfn))
            write.table(df_list[[nm]], file = fullfn, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = write_rownames, col.names = TRUE)
        }
    }

    return(invisible(NULL))
}

#' jams_rowSds(x, ...)
#'
#' Internal. Row-wise standard deviation for dense or sparse (Matrix) matrices and data frames.
#' Same arithmetic as genefilter::rowSds(), which used to mask matrixStats::rowSds() when JAMS
#' attached genefilter; inlined so genefilter is no longer needed.
#' @keywords internal
jams_rowSds <- function(x, ...){
    n <- rowSums(!is.na(x))
    n[n <= 1] <- NA
    rowvars <- rowSums((x - rowMeans(x, ...))^2, ...) / (n - 1)

    return(sqrt(rowvars))
}


#' merge_option_lists(x = NULL, y = NULL)
#'
#' Returns list x with any named elements of list y that x lacks added to it (x wins on conflicts).
#' Used to fill command line options with defaults. Same behaviour as RCurl::merge.list(), which JAMS
#' used to get by attaching RCurl; RCurl is no longer a JAMS dependency.
#' @export
merge_option_lists <- function(x = NULL, y = NULL){
    if (length(x) == 0){
        return(y)
    }
    if (length(y) == 0){
        return(x)
    }
    i <- is.na(match(names(y), names(x)))
    if (any(i)){
        x[names(y)[which(i)]] <- y[which(i)]
    }

    return(x)
}
