#' make_colour_dictionary(variable_list = NULL, pheno = NULL, phenolabels = NULL, class_to_ignore = "N_A", colour_of_class_to_ignore = "#bcc2c2", colour_table = NULL, shuffle = FALSE, print_palette = FALSE, object_to_return = "cdict")
#'
#' Returns a dictionary of colours for each class within each variable_list or phenotable
#' @export

make_colour_dictionary <- function(variable_list = NULL, pheno = NULL, phenolabels = NULL, class_to_ignore = "N_A", colour_of_class_to_ignore = "#bcc2c2", colour_table = NULL, shuffle = FALSE, print_palette = FALSE, object_to_return = "cdict"){

    atoms <- c("0200FC", "ff7f02", "468B00", "CD2626", "267ACD", "640032", "FC00FA", "000000", "478559", "161748", "f95D9B", "39A0CA")

    JAMSpalette <- paste0("#", atoms)
    for (bset in c("Set1", "Paired", "Set2", "Set3")){
        JAMSpalette <- c(JAMSpalette, colorRampPalette(brewer.pal(8, bset))(length(atoms)))
    }

    if (shuffle == TRUE){
        JAMSpalette <- JAMSpalette[shuffle(JAMSpalette)]
    }

    if (print_palette){
        pdf("JAMSpalette_values.pdf", paper = "a4r")
        scales::show_col(JAMSpalette)
        dev.off()
    }

    if (is.null(variable_list)){
        if (is.null(phenolabels)){
            flog.info("Imputing types of columns on the metadata.")
            Var_label <- colnames(pheno)
            Var_type <- sapply(Var_label, function (x) { infer_column_type(phenotable = pheno, colm = x, class_to_ignore = class_to_ignore) } )
            phenolabels <- data.frame(Var_label = unname(Var_label), Var_type = unname(Var_type), stringsAsFactors = FALSE)
        }

        #Define kinds of variables
        variable_list <- define_kinds_of_variables(phenolabels = phenolabels, phenotable = pheno, maxclass = 15, maxsubclass = 6, class_to_ignore = class_to_ignore)
    }

    #See how many different things there are
    discretes <- unique(unname(unlist(variable_list[which(!(names(variable_list) %in% c("sample", "continuous")))])))
    uniqueclasses <- unique(unlist(pheno[ , discretes]))

    #Attribute colours to a non-redundant list of classes
    if (length(uniqueclasses) > length(JAMSpalette)){
        cat("There are far too many different classes in the metadata. Humans cannot possibly distinguish this number of colours.\n")

        return(NULL)

    } else {
        #Start by attributing colours to things that are named the same
        discrete2colour <- data.frame(Name = uniqueclasses, Colour = JAMSpalette[1:length(uniqueclasses)], stringsAsFactors = FALSE)
        rownames(discrete2colour) <- uniqueclasses

        #Make class_to_ignore gray
        if (!is.null(class_to_ignore)){
            if (class_to_ignore %in% discrete2colour$Name){
                discrete2colour$Colour[(discrete2colour$Name %in% class_to_ignore)] <- colour_of_class_to_ignore
            }
        }

        if (!is.null(colour_table)){
            rownames(colour_table) <- colour_table$Class_label
            colour_table <- colour_table[rownames(colour_table) %in% discrete2colour$Name, ]
            if (nrow(colour_table) > 0){
                discrete2colour$Colour[match(colour_table$Class_label, discrete2colour$Name)]<-colour_table$Class_colour
            }
        }

        cdict <- NULL
        cdict <- list()
        for(d in 1:length(discretes)){
            classnames <- sort(unique(pheno[ , discretes[d]]))
            #if (!is.null(class_to_ignore)){
            #    classnames <- classnames[!(classnames %in% class_to_ignore)]
            #}

            ctab <- subset(discrete2colour, Name %in% classnames)
            ctab$Hex <- ctab$Colour
            cdict[[d]] <- ctab

            names(cdict)[d] <- discretes[d]
        }

        if (object_to_return == "ctable"){

            ctable_new <- plyr::rbind.fill(cdict)
            ctable_new <- ctable_new[!duplicated(ctable_new$Name), ]
            rownames(ctable_new) <- ctable_new$Name

            return(ctable_new)

        } else {
            return(cdict)
        }
    }
}
