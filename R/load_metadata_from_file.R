#' load_metadata_from_file(xlsxFile = NULL, phenotable_tsv = NULL, phenolabels_tsv = NULL)
#'
#' Loads an Excel(TM) file or tsv files and returns a list containing 1) the metadata itself and 2) description of label types for the metadata, provided the Excel file contains these data on sheets 1 and 2, respsctively.
#' @export

load_metadata_from_file <- function(opt = NULL, xlsxFile = NULL, phenotable_tsv = NULL, phenolabels_tsv = NULL, class_to_ignore = "N_A"){

    ctable <- NULL

    if (!is.null(opt)){
        xlsxFile <- opt$excel_metadata
        phenotable_tsv <- opt$phenotable
        phenolabels_tsv <- opt$phenolabels
        class_to_ignore <- opt$class_to_ignore
    }

    if (!is.null(xlsxFile)) {

        xlMD <-  list()
        for (sheet in 1:3){
            xlMD[[sheet]] <- tryCatch(read.xlsx(xlsxFile, sheet = sheet, startRow = 1, colNames = TRUE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE, skipEmptyCols = TRUE), error = function(e) { message(paste("Sheet", sheet, "is empty.")) } )
##Not working at the moment so will manually fix dates after beta
            #Coerce dates, numeric, all to character. JAMS fixes class on the fly but need this for Summarized Experiment object to work.
            #if(nrow(xlMD[[sheet]]) > 0){
            #    for (colm in 1:ncol(xlMD[[sheet]])){
            #        xlMD[[sheet]][ , colm] <- as.character(xlMD[[sheet]][ , colm])
          #          xlMD[[sheet]][is.na(xlMD[[sheet]][ , colm]), colm] <- "N_A"
        #        }
        #    }
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
                    colour <- col2hex(colour)
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
        phenotable <- read.table(file = phenotable_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        phenotable[] <- lapply(phenotable, as.character)
        if(!is.null(phenolabels_tsv)){
            phenolabels <- read.table(file = phenolabels_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            phenolabels[] <- lapply(phenolabels, as.character)
        } else {
            flog.info("No phenolabels file found. Attempting to guess which one is the Sample column.")
            #Not ideal, but guess everything as being discrete
            phenolabels <- data.frame(Var_label = colnames(phenotable), Var_type = rep("discrete", ncol(phenotable)), stringsAsFactors = FALSE)
            sampcol <- which(colnames(phenotable) == "Sample")
            if (length(sampcol) == 0){
                #Try grepping
                sampcol <- grep("Sample", colnames(phenotable), ignore.case = TRUE)[1]
            }
            phenolabels$Var_type[sampcol[1]] <- "Sample"
        }
    }

    #If there are any empty cells, fill them with JAMS "N_A" for missing data
    phenotable[is.na(phenotable)] <- class_to_ignore

    #Trim whitespace (common in Excel spreadsheets)
    phenotable <- trim_whitespace_from_df(df = phenotable)
    phenolabels <- trim_whitespace_from_df(df = phenolabels)

    #Coerce all of the metadata initially to character.
    phenotable[] <- lapply(phenotable, as.character)
    if ("date" %in% phenolabels$Var_type){
        #Adjust metadata to date
        colsToDate <- phenolabels$Var_label[which(phenolabels$Var_type == "date")]
        for (colm in 1:length(colsToDate)){
            phenotable[ , colsToDate[colm]] <- convertToDate(phenotable[ , colsToDate[colm]])
            phenotable[ , colsToDate[colm]] <- as.character(phenotable[ , colsToDate[colm]])
        }
    }

    #Check for and eliminate duplicatesq
    sampcolname <- phenolabels$Var_label[which(phenolabels$Var_type == "Sample")]
    dupes <- phenotable[ , sampcolname][duplicated(phenotable[ , sampcolname])]
    if (length(dupes) > 0){
        phenotable <- phenotable[which(!(phenotable[ , sampcolname] %in% dupes)), ]

        flog.warn(paste("There are duplicated sample prefixes in the metadata. Samples", paste0(dupes, collapse = ", "), "have been excluded from metadata. Will continue with building experiments. If you still want these samples, fix metadata and relaunch JAMSbeta."))
    }

    if (!is.null(opt)){
        opt$phenotable <- phenotable
        opt$phenolabels <- phenolabels
        opt$dupes <- opt$dupes
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
