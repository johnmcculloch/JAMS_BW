#' generate_JAMS16_methods(opt = NULL, rawfastqsdf = NULL, trunclen_used = NULL,
#'                         samples_requested = NULL, samples_used = NULL,
#'                         subsampling_triggered = FALSE, outfile = NULL)
#'
#' Synthesises a prose Methods paragraph describing the choices applied in a
#' JAMS16 run, omitting any sentence whose inputs were not supplied.
#' @export

generate_JAMS16_methods <- function(opt = NULL, rawfastqsdf = NULL, trunclen_used = NULL,
                                    samples_requested = NULL, samples_used = NULL,
                                    subsampling_triggered = FALSE, outfile = NULL){

    # ---- helpers -------------------------------------------------------------
    pkg_ver <- function(p){
        tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
    }
    # join sentences into a space-separated paragraph, dropping NULL/empty
    para <- function(x){ paste0(stats::na.omit(unlist(x)), collapse = " ") }

    S <- list()   # sentences, added conditionally

    # ---- 1. Cohort and input -------------------------------------------------
    if (!is.null(samples_used)){
        n_used <- length(samples_used)
        base_cohort <- paste0("Amplicon sequence data from ", n_used, " sample",
                              if (n_used != 1) "s" else "", " were processed.")
        if (!is.null(samples_requested) && length(samples_requested) > n_used){
            n_drop <- length(samples_requested) - n_used
            base_cohort <- paste0("Amplicon sequence data from ", n_used,
                " of ", length(samples_requested), " samples described in the metadata were processed (",
                n_drop, " sample", if (n_drop != 1) "s were" else " was",
                " excluded for missing read files or for failing read-count filtering).")
        }
        S$cohort <- base_cohort
    }

    # ---- 2. Read depth -------------------------------------------------------
    if (!is.null(rawfastqsdf) && "Number_Read_Pairs_Input" %in% colnames(rawfastqsdf)){
        depths <- suppressWarnings(as.numeric(unique(
            rawfastqsdf[, c("Prefix", "Number_Read_Pairs_Input")])$Number_Read_Pairs_Input))
        depths <- depths[!is.na(depths)]
        if (length(depths) > 0){
            S$depth <- paste0("The median number of input read pairs per sample was ",
                format(round(median(depths)), big.mark = ","), " (range ",
                format(min(depths), big.mark = ","), " to ",
                format(max(depths), big.mark = ","), ").")
        }
    }

    # ---- 3. Subsampling (only if a cap was set AND it fired) -----------------
    if (!is.null(opt$max_read_count_threshold) && isTRUE(subsampling_triggered)){
        S$subsample <- paste0("Samples exceeding ",
            format(opt$max_read_count_threshold, big.mark = ","),
            " read pairs were randomly subsampled to that depth (seqtk) prior to denoising.")
    }

    # ---- 4. Primer removal (only if primers supplied) ------------------------
    primers_present <- !is.null(opt$forward_primer) && !is.null(opt$reverse_primer)
    if (primers_present){
        S$primers <- paste0("Primers were removed with the DADA2 removePrimers function (forward primer ",
            opt$forward_primer, ", reverse primer ", opt$reverse_primer, ").")
    }

    # ---- 5. Quality filtering and truncation ---------------------------------
    if (!is.null(trunclen_used) && length(trunclen_used) == 2){
        trunc_clause <- paste0("Reads were quality-filtered and truncated with DADA2 filterAndTrim ",
            "(truncLen = ", trunclen_used[1], " for R1 and ", trunclen_used[2], " for R2, ",
            "maxEE = 2,2, truncQ = 2, maxN = 0, with PhiX removal)")
        if (primers_present){
            trunc_clause <- paste0(trunc_clause,
                "; truncation lengths were reduced where necessary so as not to exceed read length after primer removal")
        }
        S$filter <- paste0(trunc_clause, ".")
    }
    if (!is.null(opt$min_read_count_threshold)){
        S$minreads <- paste0("Samples with fewer than ",
            format(opt$min_read_count_threshold, big.mark = ","),
            " read pairs after filtering were excluded from downstream analysis.")
    }

    # ---- 6. Core DADA2 pipeline (always runs) --------------------------------
    S$dada2 <- paste0("Error rates were learned per run and amplicon sequence variants (ASVs) were ",
        "inferred with DADA2; forward and reverse reads were merged, a sequence table constructed, ",
        "and chimeras removed with removeBimeraDenovo (consensus method).")

    # ---- 7. Taxonomy (branch on method) --------------------------------------
    db_db   <- if (!is.null(opt$trainingset)) basename(opt$trainingset) else NA
    db_train <- if (!is.null(opt$silva_trainset)) basename(opt$silva_trainset) else NA
    db_sp    <- if (!is.null(opt$silva_species_assignment)) basename(opt$silva_species_assignment) else NA

    if (!is.null(opt$taxamethod)){
        S$taxonomy <- switch(opt$taxamethod,
            "dada2" = {
                parts <- "Taxonomy was assigned with the DADA2 naive Bayesian classifier (assignTaxonomy) to genus level"
                if (!is.na(db_train)) parts <- paste0(parts, " against the ", db_train, " reference")
                parts <- paste0(parts, ", followed by exact-match species assignment (addSpecies")
                if (!is.na(db_sp)) parts <- paste0(parts, " using ", db_sp)
                paste0(parts, ").")
            },
            "dada2pr2" = paste0("Taxonomy was assigned with the DADA2 naive Bayesian classifier (assignTaxonomy) ",
                "using the PR2 reference database", if (!is.na(db_db)) paste0(" (", db_db, ")") else "",
                ", following the PR2 rank scheme (Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species)."),
            "IDtaxa" = paste0("Taxonomy was assigned with the IDTAXA algorithm (DECIPHER) at a confidence ",
                "threshold of 60", if (!is.na(db_db)) paste0(" against the ", db_db, " training set") else "", "."),
            "IDtaxapr2" = paste0("Taxonomy was assigned with the IDTAXA algorithm (DECIPHER) at a confidence ",
                "threshold of 60 using the PR2 training set",
                if (!is.na(db_db)) paste0(" (", db_db, ")") else "",
                ", following the PR2 rank scheme."),
            NULL)
    }

    # ---- 8. Software versions ------------------------------------------------
    vers <- c(dada2 = pkg_ver("dada2"))
    if (!is.null(opt$taxamethod) && opt$taxamethod %in% c("IDtaxa", "IDtaxapr2")){
        vers <- c(vers, DECIPHER = pkg_ver("DECIPHER"))
    }
    vers <- c(vers, JAMS = pkg_ver("JAMS"))
    vers <- vers[!is.na(vers)]
    if (length(vers) > 0){
        verstr <- paste0(names(vers), " v", vers, collapse = ", ")
        S$versions <- paste0("Analyses were performed with ", verstr,
            " under ", R.version.string, ".")
    }

    # ---- assemble ------------------------------------------------------------
    methods_text <- para(S)

    if (!is.null(outfile)){
        writeLines(methods_text, con = outfile)
        flog.info(paste("Wrote methods narrative to", outfile))
    }

    return(invisible(list(text = methods_text, sentences = S, versions = vers)))
}


#' tabulate_JAMS16_methods(opt = NULL, rawfastqsdf = NULL, trunclen_used = NULL,
#'                         trunclen_source = "user_specified",
#'                         samples_requested = NULL, samples_used = NULL,
#'                         subsampling_triggered = FALSE, outfile = NULL)
#'
#' Builds a tidy parameter-value data frame of every methods-relevant choice in a
#' JAMS16 run, for machine-readable provenance alongside the prose narrative.
#' Rows whose underlying value is absent are omitted, so the table never carries
#' empty or misleading cells. Optionally writes it to TSV.
#' @export

tabulate_JAMS16_methods <- function(opt = NULL, rawfastqsdf = NULL, trunclen_used = NULL,
                                    trunclen_source = "user_specified",
                                    samples_requested = NULL, samples_used = NULL,
                                    subsampling_triggered = FALSE, outfile = NULL){

    pkg_ver <- function(p){ tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_) }

    rows <- list()
    #Small accumulator: only lands a row when the value is non-NULL / non-NA / non-empty.
    add <- function(parameter, value, category = NA_character_){
        if (is.null(value)) return(invisible(NULL))
        value <- as.character(value)
        if (length(value) == 0 || all(is.na(value)) || all(value == "")) return(invisible(NULL))
        rows[[length(rows) + 1]] <<- data.frame(
            Category  = category,
            Parameter = parameter,
            Value     = paste0(value, collapse = ", "),
            stringsAsFactors = FALSE)
    }

    # ---- Project / cohort ----------------------------------------------------
    add("Project_name", opt$projectname, "Project")
    if (!is.null(samples_requested)) add("Samples_in_metadata", length(samples_requested), "Cohort")
    if (!is.null(samples_used))      add("Samples_analysed",    length(samples_used),      "Cohort")
    if (!is.null(samples_requested) && !is.null(samples_used)){
        add("Samples_excluded", (length(samples_requested) - length(samples_used)), "Cohort")
    }

    # ---- Read depth ----------------------------------------------------------
    if (!is.null(rawfastqsdf) && "Number_Read_Pairs_Input" %in% colnames(rawfastqsdf)){
        d <- suppressWarnings(as.numeric(unique(
            rawfastqsdf[, c("Prefix", "Number_Read_Pairs_Input")])$Number_Read_Pairs_Input))
        d <- d[!is.na(d)]
        if (length(d) > 0){
            add("Median_input_read_pairs", round(median(d)), "Read_depth")
            add("Min_input_read_pairs",    min(d),           "Read_depth")
            add("Max_input_read_pairs",    max(d),           "Read_depth")
        }
    }

    # ---- Read-count thresholds ----------------------------------------------
    add("Min_read_count_threshold", opt$min_read_count_threshold, "Filtering")
    add("Max_read_count_threshold", opt$max_read_count_threshold, "Filtering")
    add("Subsampling_applied", if (!is.null(opt$max_read_count_threshold) && isTRUE(subsampling_triggered)) "yes" else "no", "Filtering")

    # ---- Primers -------------------------------------------------------------
    add("Forward_primer", opt$forward_primer, "Primers")
    add("Reverse_primer", opt$reverse_primer, "Primers")
    if (!is.null(opt$forward_primer) && !is.null(opt$reverse_primer)){
        #primer_source lets a reader see whether these came from -V4 or were typed in.
        add("Primer_source", if (!is.null(opt$V4)) "V4_preset" else "user_supplied", "Primers")
    }

    # ---- Truncation and filtering --------------------------------------------
    if (!is.null(trunclen_used) && length(trunclen_used) == 2){
        add("TruncLen_R1", trunclen_used[1], "Filtering")
        add("TruncLen_R2", trunclen_used[2], "Filtering")
        add("TruncLen_source", trunclen_source, "Filtering")
    }
    #filterAndTrim constants, hard-wired in the script; recording them makes the row self-contained.
    add("MaxEE_R1_R2", "2,2", "Filtering")
    add("TruncQ",      "2",   "Filtering")
    add("MaxN",        "0",   "Filtering")
    add("PhiX_removed","yes", "Filtering")

    # ---- Taxonomy ------------------------------------------------------------
    add("Taxonomy_method", opt$taxamethod, "Taxonomy")
    if (!is.null(opt$taxamethod)){
        if (opt$taxamethod == "dada2"){
            add("Taxonomy_trainset",  basename_safe(opt$silva_trainset),           "Taxonomy")
            add("Taxonomy_species_db", basename_safe(opt$silva_species_assignment), "Taxonomy")
            add("Taxonomy_trainset_path",  opt$silva_trainset,            "Taxonomy")
            add("Taxonomy_species_db_path", opt$silva_species_assignment, "Taxonomy")
        } else {
            add("Taxonomy_trainset",      basename_safe(opt$trainingset), "Taxonomy")
            add("Taxonomy_trainset_path", opt$trainingset,                "Taxonomy")
            if (opt$taxamethod %in% c("IDtaxa", "IDtaxapr2")){
                add("IDTAXA_confidence_threshold", "60", "Taxonomy")
            }
        }
    }

    # ---- Merged amplicon length (free, informative) --------------------------
    if (!is.null(rawfastqsdf) && "Median_Amplicon_Length_bp" %in% colnames(rawfastqsdf)){
        aml <- suppressWarnings(as.numeric(rawfastqsdf$Median_Amplicon_Length_bp))
        aml <- aml[!is.na(aml)]
        if (length(aml) > 0) add("Median_merged_amplicon_length_bp", round(median(aml)), "Amplicon")
    }

    # ---- Software versions ---------------------------------------------------
    add("Version_dada2", pkg_ver("dada2"), "Software")
    if (!is.null(opt$taxamethod) && opt$taxamethod %in% c("IDtaxa", "IDtaxapr2")){
        add("Version_DECIPHER", pkg_ver("DECIPHER"), "Software")
    }
    add("Version_JAMS", pkg_ver("JAMS"), "Software")
    add("Version_R",    paste(R.version$major, R.version$minor, sep = "."), "Software")
    add("Run_timestamp", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "Software")

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
