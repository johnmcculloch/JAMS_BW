#' load_metadata_from_file(xlsxFile = NULL, phenotable_tsv = NULL, phenolabels_tsv = NULL)
#'
#' Loads an Excel(TM) file or tsv files and returns a list containing 1) the metadata itself and 2) description of label types for the metadata, provided the Excel file contains these data on sheets 1 and 2, respsctively.
#'
#' Column names are read exactly as written and then conformed to JAMS style by enforce_JAMS_style_metadata() (q.v.), which also trims whitespace, fills missing cells with class_to_ignore, and drops empty rows, rows without a sample name and duplicated samples. Every column rename is logged, one line each (so it appears in the _JAMS.log file of JAMSbeta and JAMS16), and, when opt is passed, the table of original -> new names is written to <outdir>/<projectname>_metadata_column_name_changes.tsv and kept in opt$metadata_column_name_changes. phenolabels entries are renamed to match, whether they were written with the original names or with the dotted names produced by earlier JAMS versions.
#' @export

load_metadata_from_file <- function(opt = NULL, xlsxFile = NULL, phenotable_tsv = NULL, phenolabels_tsv = NULL, class_to_ignore = "N_A"){

    ctable <- NULL

    if (!is.null(opt)){
        xlsxFile <- opt$excel_metadata
        phenotable_tsv <- opt$phenotable
        phenolabels_tsv <- opt$phenolabels
        class_to_ignore <- opt$class_to_ignore
    }

    #opt does not necessarily carry a class_to_ignore, and the line above would then set it to
    #NULL, which propagates into every downstream sentinel comparison. Fall back to the default.
    if (is.null(class_to_ignore) || length(class_to_ignore) < 1){
        class_to_ignore <- "N_A"
    }

    if (!is.null(xlsxFile)) {

        xlMD <-  list()
        for (sheet in 1:3){
            cursheet <- tryCatch(read.xlsx(xlsxFile, sheet = sheet, startRow = 1, colNames = TRUE, rowNames = FALSE, detectDates = TRUE, skipEmptyRows = FALSE, skipEmptyCols = TRUE), error = function(e) {
                flog.info(paste0("Metadata spreadsheet sheet ", sheet, " could not be read (", conditionMessage(e), "); assuming it is absent or empty."))
                NULL
            })
            if (is.null(cursheet) || (ncol(cursheet) == 0)){
                next
            }

            #read.xlsx turns spaces in the header into dots. Recover the header exactly as typed, so that all name
            #fixing happens in one place, enforce_JAMS_style_metadata() below, and is logged there. The raw header is
            #only used if it is consistent with what read.xlsx returned (same columns, differing only by spaces vs dots);
            #otherwise read.xlsx's names are kept.
            rawhdr <- tryCatch(read.xlsx(xlsxFile, sheet = sheet, rows = 1, colNames = FALSE, rowNames = FALSE, skipEmptyCols = TRUE), error = function(e) { NULL })
            if (!is.null(rawhdr) && (ncol(rawhdr) == ncol(cursheet))){
                rawnames <- trimws(as.character(unlist(rawhdr[1, ], use.names = FALSE)))
                if (!anyNA(rawnames) && identical(gsub(" ", ".", rawnames, fixed = TRUE), gsub(" ", ".", trimws(colnames(cursheet)), fixed = TRUE))){
                    colnames(cursheet) <- rawnames
                } else {
                    flog.info(paste0("Could not recover the exact header of metadata sheet ", sheet, "; using the column names as returned by read.xlsx (spaces turned into dots)."))
                }
            }

            #All character, with blanks set to the missing-data sentinel
            cursheet[] <- lapply(cursheet, function(v){
                v <- as.character(v)
                v[is.na(v)] <- class_to_ignore[1]
                v
            })
            xlMD[[length(xlMD) + 1]] <- cursheet
        }
        if (length(xlMD) == 0){
            stop(paste("No readable sheets were found in metadata spreadsheet", xlsxFile))
        }
        #Find out who is who
        find_metadata_type <- function(xlMD = NULL, sheet = NULL){
            if (all(c("Class_label", "Class_colour") %in% colnames(xlMD[[sheet]]))){
                mdtype <- "ctable"
            } else if (all(c("Var_label", "Var_type") %in% colnames(xlMD[[sheet]]))){
                mdtype <- "phenolabels"
            } else {
                mdtype <- "phenotable"
            }

            return(mdtype)
        }

        metadata_types <- sapply(1:length(xlMD), function (x) { find_metadata_type(xlMD, sheet = x) })

        #Get a phenotable
        if (unname(table(metadata_types)["phenotable"]) != 1){
            stop("Unable to determine which sheet is the phenotable. Review metadata and try again.")
        } else {
            phenotable <- xlMD[[which(metadata_types == "phenotable")]]
        }

        #See if colour table exists
        #Get a phenotable
        if (!is.na(unname(table(metadata_types)["ctable"]))) {
            flog.info("Found a colour table for mapping classes onto colours (and shapes)")
            ctable <- xlMD[[which(metadata_types == "ctable")[1]]]
            ctable <- trim_whitespace_from_df(ctable)
            fix_hex_cols <- function(colour){
                if (length(grep("#", colour, fixed = TRUE)) == 0){
                    #Equivalent of gplots::col2hex(), using grDevices only.
                    colour <- grDevices::rgb(t(grDevices::col2rgb(colour)), maxColorValue = 255)
                }
                return(colour)
            }
            ctable$Class_colour <- sapply(ctable$Class_colour, function (x) { fix_hex_cols(x) } )
        }

        #See if phenolabels exists, else make one
        if (!is.na(unname(table(metadata_types)["phenolabels"]))) {
            phenolabels <- xlMD[[which(metadata_types == "phenolabels")[1]]]
            phenolabels <- trim_whitespace_from_df(phenolabels)
        } else {
            Var_label <- colnames(phenotable)
            Var_type <- sapply(Var_label, function (x) { infer_column_type(phenotable = phenotable, colm = x, class_to_ignore = class_to_ignore) } )
            phenolabels <- data.frame(Var_label = unname(Var_label), Var_type = unname(Var_type), stringsAsFactors = FALSE)
        }

    } else {
        #check.names = FALSE keeps the column names exactly as written (the default runs make.names(), turning
        #"IL-6 (pg/mL)" into "IL.6..pg.mL."); name fixing is done, and logged, by enforce_JAMS_style_metadata() below.
        #colClasses = "character" keeps values as written (e.g. "007" stays "007"); comment.char = "" stops a "#"
        #in a column name or cell from silently truncating the rest of the line; quote = "\"" (rather than the default,
        #which also includes the apostrophe) stops a "Crohn's" in a cell from swallowing tabs and newlines.
        strip_bom <- function(df){ colnames(df) <- sub("^\ufeff", "", colnames(df)); df }
        phenotable <- read.table(file = phenotable_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", comment.char = "", quote = "\"")
        phenotable <- strip_bom(phenotable)
        phenotable[] <- lapply(phenotable, as.character)
        if(!is.null(phenolabels_tsv)){
            phenolabels <- read.table(file = phenolabels_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", comment.char = "", quote = "\"")
            phenolabels <- strip_bom(phenolabels)
            phenolabels[] <- lapply(phenolabels, as.character)
        } else {
            flog.info("No phenolabels file found. Attempting to guess which one is the Sample column.")
            #Not ideal, but guess everything as being discrete
            phenolabels <- data.frame(Var_label = colnames(phenotable), Var_type = rep("discrete", ncol(phenotable)), stringsAsFactors = FALSE)
            sampcol <- which(colnames(phenotable) == "Sample")
            if (length(sampcol) == 0){
                #Exact match ignoring case first (so "SampleID" cannot beat "sample"), then any column containing "sample"
                sampcol <- which(tolower(trimws(colnames(phenotable))) == "sample")
                if (length(sampcol) == 0){
                    sampcol <- grep("Sample", colnames(phenotable), ignore.case = TRUE)
                }
            }
            if (!is.na(sampcol[1])){
                phenolabels$Var_type[sampcol[1]] <- "Sample"
            }
        }
    }

    #### Enforce JAMS style on the metadata ####
    #One place, one set of rules, for every way metadata enters JAMS: column names made safe (and every rename logged,
    #one line each, to the _JAMS.log file), whitespace trimmed, every cell character with blanks/NA set to
    #class_to_ignore, empty rows and rows without a sample name dropped, duplicated samples dropped (as before).
    #phenolabels is renamed in step but NOT reconciled, as JAMSbeta uses it as the list of columns to keep.
    sampcolname <- phenolabels$Var_label[which(phenolabels$Var_type == "Sample")]
    sampcolname <- if (length(sampcolname) >= 1) sampcolname[1] else "Sample"
    flog.info(paste0("Enforcing JAMS style on the metadata (sample column \"", sampcolname, "\")."))
    enforced <- enforce_JAMS_style_metadata(pheno = phenotable, phenolabels = phenolabels, sample_column = sampcolname, class_to_ignore = class_to_ignore, fix_column_names = TRUE, reconcile_phenolabels = FALSE, on_duplicate_samples = "drop")
    phenotable <- enforced$phenotable
    phenolabels <- enforced$phenolabels
    namechanges <- attr(phenotable, "JAMS_column_name_changes")
    dupes <- attr(phenotable, "JAMS_duplicated_samples_dropped")
    final_sampcol <- attr(phenotable, "JAMS_sample_column")
    if (length(dupes) > 0){
        flog.warn("There were duplicated sample prefixes in the metadata (see above); those samples have been excluded. If you still want them, fix the metadata and relaunch.")
    }

    #JAMSbeta and JAMS16 address the sample column literally as "Sample", so whatever it was called, call it that
    if (!identical(final_sampcol, "Sample")){
        if ("Sample" %in% colnames(phenotable)){
            flog.warn(paste0("The sample names are in column \"", final_sampcol, "\", but there is also a different column called \"Sample\". JAMS expects the sample names to be in a column called \"Sample\"; please rename one of them in the metadata."))
        } else {
            flog.warn(paste0("enforce_JAMS_style_metadata: metadata column renamed: \"", final_sampcol, "\" -> \"Sample\" (JAMS requires the sample names to be in a column called \"Sample\")."))
            colnames(phenotable)[colnames(phenotable) == final_sampcol] <- "Sample"
            phenolabels$Var_label[phenolabels$Var_label %in% final_sampcol] <- "Sample"
            if (final_sampcol %in% namechanges$New){
                namechanges$New[namechanges$New == final_sampcol] <- "Sample"
            } else {
                namechanges <- rbind(namechanges, data.frame(Original = final_sampcol, New = "Sample", stringsAsFactors = FALSE))
            }
            attr(phenotable, "JAMS_column_name_changes") <- namechanges
        }
    }

    #Keep a permanent, machine-readable record of any column renames next to the log file
    if (!is.null(opt) && !is.null(opt$outdir) && (nrow(namechanges) > 0)){
        prefix <- if (!is.null(opt$projectname)) opt$projectname else "JAMS"
        namechanges_fn <- file.path(opt$outdir, paste0(prefix, "_metadata_column_name_changes.tsv"))
        wrote <- tryCatch({
            utils::write.table(namechanges, file = namechanges_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
            TRUE
        }, error = function(e){ FALSE })
        if (wrote){
            flog.info(paste0("The full table of original -> new metadata column names has been written to ", namechanges_fn))
        } else {
            flog.warn(paste0("Could not write the table of metadata column name changes to ", namechanges_fn, "; the renames are listed above in this log."))
        }
    }

    if ("date" %in% phenolabels$Var_type){
        #Adjust metadata to date. Only cells which actually hold data are converted: the
        #class_to_ignore sentinel must survive untouched, as convertToDate() would otherwise
        #turn it into a real NA and undo the normalisation performed immediately above.
        colsToDate <- phenolabels$Var_label[which(phenolabels$Var_type == "date")]
        missing_date_cols <- colsToDate[!(colsToDate %in% colnames(phenotable))]
        if (length(missing_date_cols) > 0){
            flog.warn(paste0("phenolabels declares Var_type \"date\" for column(s) not present in the metadata, which will be skipped: ", paste0(missing_date_cols, collapse = ", ")))
        }
        colsToDate <- colsToDate[colsToDate %in% colnames(phenotable)]
        for (colm in seq_along(colsToDate)){
            curr_col <- phenotable[ , colsToDate[colm]]
            usable <- which(!(curr_col %in% class_to_ignore))
            if (length(usable) > 0){
                converted <- as.character(convertToDate(curr_col[usable]))
                #convertToDate() coerces its input with as.numeric(), so it yields NA for any
                #cell which is already an ISO date string ("2019-08-09") rather than an Excel
                #serial number. Say so rather than silently emptying the column.
                if (any(is.na(converted))){
                    flog.warn(paste0("Metadata column ", colsToDate[colm], " is declared with Var_type \"date\", but ", sum(is.na(converted)), " of its ", length(usable), " value(s) could not be read as an Excel serial date and have been set to ", class_to_ignore[1], ". If this column already holds ISO dates such as 2019-08-09, declare it as discrete instead."))
                    converted[is.na(converted)] <- class_to_ignore[1]
                }
                curr_col[usable] <- converted
            }
            phenotable[ , colsToDate[colm]] <- curr_col
        }
    }

    if (!is.null(opt)){
        opt$phenotable <- phenotable
        opt$phenolabels <- phenolabels
        opt$dupes <- dupes
        opt$metadata_column_name_changes <- namechanges
        opt$ctable <- ctable

        return(opt)
    } else {
        metadata <- list()
        metadata[[1]] <- phenotable
        metadata[[2]] <- phenolabels
        metadata[[3]] <- ctable
        names(metadata) <- c("phenotable", "phenolabels", "ctable")[unname(sapply(c("phenotable", "phenolabels", "ctable"), function(x) { !is.null(get(x))} ))]

        return(metadata)
    }
}
