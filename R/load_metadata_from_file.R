#' load_metadata_from_file(xlsxFile = NULL, phenotable_tsv = NULL, phenolabels_tsv = NULL)
#'
#' Loads an Excel(TM) file or tsv files and returns a list containing 1) the metadata itself and 2) description of label types for the metadata, provided the Excel file contains these data on sheets 1 and 2, respsctively.
#' @export

load_metadata_from_file <- function(opt = NULL, xlsxFile = NULL, phenotable_tsv = NULL, phenolabels_tsv = NULL, class_to_ignore = "N_A"){

    if (!is.null(opt)){
        xlsxFile <- opt$excel_metadata
        phenotable_tsv <- opt$phenotable
        phenolabels_tsv <- opt$phenolabels
        class_to_ignore <- opt$class_to_ignore
    }

    if (!is.null(xlsxFile)) {
        #Load metadata table from an Excel-style spreadsheet.
        phenotable <- read.xlsx(xlsxFile, sheet = 1, startRow = 1, colNames = TRUE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE, skipEmptyCols = TRUE)
        phenolabels <- read.xlsx(xlsxFile, sheet = 2, startRow = 1, colNames = TRUE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE, skipEmptyCols = TRUE)
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
        return(opt)
    } else {
        metadata <- list()
        metadata[[1]] <- phenotable
        metadata[[2]] <- phenolabels
        return(metadata)
    }

}
