#' load_metadata_from_xl
#'
#' Loads an Excel(TM) file and returns a vector containing 1) the metadata itself and 2) description of label types for the metadata, provided the Excel file contains these data on sheets 1 and 2, respsctively.
#' @export

load_metadata_from_xl<-function(xlsxFile=NULL){
    metadata<-NULL
    metadata<-vector("list",length=2)
    #Load metadata table from an Excel spreadsheet (Ugh, I know!!). Sigh.
    phenotable <- read.xlsx(xlsxFile, sheet = 1, startRow = 1, colNames = TRUE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE, skipEmptyCols = TRUE, na.strings = "NA")
    phenolabels <- read.xlsx(xlsxFile, sheet = 2, startRow = 1, colNames = TRUE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE, skipEmptyCols = TRUE, na.strings = "NA")

    #Trim whitespace (common in Excel spreadsheets)
    phenotable <- trim_whitespace_from_df(df = phenotable)
    phenolabels <- trim_whitespace_from_df(df = phenolabels)

    #Adjust phenotable for class of data
    phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="Sample")[1,1]))]<-as.character(phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="Sample")[1,1]))])

    for (v in 1:length(subset(phenolabels, Var_type=="discrete")$Var_label)){
        phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="discrete")$Var_label[v]))]<-as.character(phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="discrete")$Var_label[v]))])
    }

    for (v in 1:length(subset(phenolabels, Var_type=="ID")$Var_label)){
        phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="ID")$Var_label[v]))]<-as.character(phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="ID")$Var_label[v]))])
    }

    for (v in 1:length(subset(phenolabels, Var_type=="temporal")$Var_label)){
        phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="temporal")$Var_label[v]))]<-as.character(phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="temporal")$Var_label[v]))])
    }

    for (v in 1:length(subset(phenolabels, Var_type=="continuous")$Var_label)){
        phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="continuous")$Var_label[v]))]<-as.numeric(phenotable[ ,which(colnames(phenotable)==(subset(phenolabels, Var_type=="continuous")$Var_label[v]))])
    }

    metadata[[1]] <- phenotable
    metadata[[2]] <- phenolabels

    return(metadata)
}
