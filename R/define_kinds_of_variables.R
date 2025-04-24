#' define_kinds_of_variables(metadataXL = NULL, phenolabels = NULL, phenotable = NULL, maxclass = 10, maxsubclass = 4)
#'
#' From the metadata returns a list of variable names with the desired number of classes belonging to "sample", "discrete", "continuous", "binary" and "subsettable" types.
#' @export

#returns a list of variable names
define_kinds_of_variables <- function(metadataXL = NULL, phenolabels = NULL, phenotable = NULL, maxclass = 15, maxsubclass = 4, class_to_ignore = "N_A", verbose = TRUE){

    if (!(is.null(metadataXL))){
        metadata <- load_metadata_from_file(xlsxFile = metadataXL)
        phenotable <- metadata[[1]]
        phenolabels <- metadata[[2]]
    }

    if (is.null(phenolabels)){
        if (verbose){
            flog.info("Imputing types of columns on the metadata.")
        }
        Var_label <- colnames(phenotable)
        Var_type <- sapply(Var_label, function (x) { infer_column_type(phenotable = phenotable, colm = x, class_to_ignore = class_to_ignore) } )
        phenolabels <- data.frame(Var_label = unname(Var_label), Var_type = unname(Var_type), stringsAsFactors = FALSE)
    }

    validcols <- as.character(phenolabels$Var_label)
    if (verbose){
        flog.info(paste("Phenolabels contains", length(validcols), "categories."))
    }

    #Stop if you are asking for more than you have.
    #if (length(validcols) > ncol(phenotable)){
    #    stop("Phenolabels has more categories than the metadata in phenotable. Review and try again.")
    #}

    ptsampcol <- as.character(phenolabels$Var_label[which(phenolabels$Var_type == "Sample")])
    if (length(ptsampcol) != 1){
        stop("You must define exactly one column as being the Sample column in your metadata. Please see documenation.")
    } else {
        if (verbose){
            flog.info(paste("Samples are in the", ptsampcol, "column in the metadata."))
        }
    }

    if (verbose){
        flog.info(paste("Metadata classes specified in phenolabels are:", paste0(validcols, collapse = ", ")))
        flog.info("Adjusting metadata to contain only these columns.")
    }

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

    if (verbose) {
        flog.info(paste("Ignoring cells containing", class_to_ignore))
    }

    #Find and validate subsettable variables
    variables_subs <- as.character(phenolabels$Var_label[which((phenolabels$Var_type == "subsettable"))])
    if (length(variables_subs) > 0){
        numcl <- NULL
        for(v in 1:length(variables_subs)){
            numcl[v] <- length(unique(phenotable[ , (variables_subs[v])][which(phenotable[ , (variables_subs[v])] != class_to_ignore)]))
        }
        if (verbose){
            flog.info(paste("Found subsettable variables", paste0(variables_subs, collapse=", "), "with", paste0(numcl, collapse=", "), "classes, respectively."))
        }
            variables_subs <- variables_subs[which((numcl > 1) & (numcl < (maxsubclass + 1)))]
        if (verbose){
            flog.info(paste("Keeping subsettable variables", paste0(variables_subs, collapse=", ")))
        }
    } else {
        if (verbose){
            flog.info("Found no subsettable variables.")
        }
    }
    variable_list$subsettable <- as.character(variables_subs)

    #Find and validate discrete variables
    variables_disc <- as.character(phenolabels$Var_label[which((phenolabels$Var_type %in% c("discrete", "subsettable", "date")))])
    if (length(variables_disc) > 0){
        numcl <- NULL
        for(v in 1:length(variables_disc)){
            numcl[v] <- length(unique(phenotable[ , (variables_disc[v])][which(phenotable[ , (variables_disc[v])] != class_to_ignore)]))
        }
        if (verbose){
            flog.info(paste("Found discrete variables", paste0(variables_disc, collapse = ", "), "with", paste0(numcl, collapse = ", "),  "classes, respectively."))
        }
        variables_disc <- variables_disc[which((numcl > 0) & (numcl < (maxclass + 1)))]
        if (verbose){
            flog.info(paste("Keeping discrete variables", variables_disc))
        }
    } else {
        if (verbose){
            flog.info("Found no discrete variables.")
        }
    }
    variable_list$discrete <- as.character(variables_disc)

    #Find and validate continuous variables
    variables_cont <- phenolabels$Var_label[which((phenolabels$Var_type == "continuous"))]
    if (length(variables_cont) > 0){
        #Check that valid (non-ignore values can be made numeric)
        is_numeric <- NULL
        for (v in 1:length(variables_cont)){
            is_numeric[v] <- length(which(is.na(phenotable[ , (variables_cont[v])][which(phenotable[ , (variables_cont[v])] != class_to_ignore)]) == FALSE))
            if (is_numeric[v] == 0){
                if (verbose){
                    flog.info(paste("No values which are numeric can be found in", variables_cont[v]))
                }
            }
        }
        variables_cont <- variables_cont[which(is_numeric > 0)]
        if (length(variables_cont) > 0){
            if (verbose){
                flog.info(paste("Found continuous variables", paste0(variables_cont, collapse = ", ")))
            }
        }
    } else {
        if (verbose){
            flog.info("Found no continuous variables.")
        }
    }
    variable_list$continuous <- as.character(variables_cont)

    #Find binary variables
    variables_bin <- as.character(phenolabels$Var_label[which((phenolabels$Var_type != "Sample") & (phenolabels$Var_type != "continuous"))])
    if (length(variables_bin) > 0){
        numcl <- NULL
        for(v in 1:length(variables_bin)){
            numcl[v] <- length(unique(phenotable[ , (variables_bin[v])][which(phenotable[ , (variables_bin[v])] != class_to_ignore)]))
        }
        variables_bin <- variables_bin[which(numcl == 2)]
        if (length(variables_bin) > 0){
            if (verbose){
                flog.info(paste("Found binary variables", variables_bin))
            }
        }
    }
    variable_list$binary <- as.character(variables_bin)

    return(variable_list)
}
