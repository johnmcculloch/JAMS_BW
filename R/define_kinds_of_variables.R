#' define_kinds_of_variables(phenolabels=NULL, phenotable=NULL, maxclass=10, maxsubclass=4)
#'
#' From the metadata returns a list of variable names with the desired number of classes belonging to "sample", "discrete", "continuous", "binary" and "subsettable" types.
#' @export

#returns a list of variable names
define_kinds_of_variables <- function(metadataXL = NULL, phenolabels = NULL, phenotable = NULL, maxclass = 10, maxsubclass = 4){

    if (!(is.null(metadataXL))){
        metadata <- load_metadata_from_xl(xlsxFile = metadataXL)
        phenotable <- metadata[[1]]
        phenolabels <- metadata[[2]]
    }

    validcols <-as.character(phenolabels$Var_label)
    print(paste("Phenolabels contains", length(validcols), "categories."))

    #Stop if you are asking for more than you have.
    if(length(validcols) > ncol(phenotable)){
        stop("Phenolabels has more categories than the metadata in phenotable. Review and try again.")
    }

    ptsampcol <- as.character(phenolabels$Var_label[which(phenolabels$Var_type == "Sample")])
    if (length(ptsampcol) != 1){
        stop("You must define exactly one column as being the Sample column in your metadata. Please see documenation.")
    } else {
        print(paste("Samples are in the", ptsampcol, "column in the metadata."))
    }

    print(paste("Metadata classes specified in phenolabels are:", paste0(validcols, collapse = ", ")))
    print("Adjusting metadata to contain only these columns.")

    #Stop if the classes you want in phenolabels do not exist in phenotable.
    if (all(validcols %in% colnames(phenotable))){
        phenotable <- phenotable[ , validcols]
    } else {
        stop("One or more classes specified in phenolabels is missing as (a) column(s) in phenotable. Please check files.")
    }

    vartypes <- c("sample", "discrete", "continuous", "binary", "subsettable")
    variable_list <- vector("list", length = length(vartypes))
    names(variable_list) <- vartypes

    #Define kinds of variables
    #Find the sample variable
    variable_list$sample <- as.character(phenolabels$Var_label[which((phenolabels$Var_type == "Sample"))])

    #Find and validate subsettable variables
    variables_subs <- as.character(phenolabels$Var_label[which((phenolabels$Var_type == "subsettable"))])
    if (length(variables_subs) > 0){
        numcl <- NULL
        for(v in 1:length(variables_subs)){
            numcl[v] <- length(unique(phenotable[ , (variables_subs[v])]))
        }
        print(paste("Found subsettable variables", paste0(variables_subs, collapse=", "), "with", paste0(numcl, collapse=", "),  "classes, respectively."))
        variables_subs <- variables_subs[which((numcl > 1) & (numcl < (maxsubclass + 1)))]
        print(paste("Keeping subsettable variables", paste0(variables_subs, collapse=", ")))
    } else {
        print("Found no subsettable variables.")
    }
    variable_list$subsettable <- as.character(variables_subs)

    #Find and validate discrete variables
    variables_disc <- as.character(phenolabels$Var_label[which((phenolabels$Var_type %in% c("discrete", "subsettable")))])
    if (length(variables_disc) > 0){
        numcl <- NULL
        for(v in 1:length(variables_disc)){
            numcl[v] <- length(unique(phenotable[ , (variables_disc[v])]))
        }
        print(paste("Found discrete variables", paste0(variables_disc, collapse=", "), "with", paste0(numcl, collapse=", "),  "classes, respectively."))
        variables_disc <- variables_disc[which((numcl > 1) & (numcl < (maxclass + 1)))]
        print(paste("Keeping discrete variables", variables_disc))
    } else {
        print("Found no discrete variables.")
    }
    variable_list$discrete <- as.character(variables_disc)

    #Find and validate continuous variables
    variables_cont <- phenolabels$Var_label[which((phenolabels$Var_type == "continuous"))]
    if (length(variables_cont) > 0){
        print(paste("Found continuous variables", paste0(variables_cont, collapse=", ")))
    } else {
        print("Found no continuous variables.")
    }
    variable_list$continuous <- as.character(variables_cont)

    #Find binary variables
    variables_bin <- as.character(phenolabels$Var_label[which((phenolabels$Var_type != "Sample") & (phenolabels$Var_type != "continuous"))])
    if (length(variables_bin) > 0){
        numcl <- NULL
        for(v in 1:length(variables_bin)){
            numcl[v] <- length(unique(phenotable[ , (variables_bin[v])]))
        }
        variables_bin <- variables_bin[which(numcl == 2)]
        if (length(variables_bin) > 0){
            print(paste("Found binary variables", variables_bin))
        }
    }
    variable_list$binary <- as.character(variables_bin)

    return(variable_list)
}
