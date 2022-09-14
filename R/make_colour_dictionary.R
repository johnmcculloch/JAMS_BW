#' make_colour_dictionary(variable_list = NULL, pheno = NULL, phenolabels = NULL, class_to_ignore = "N_A", colour_of_class_to_ignore = "#bcc2c2", colour_table = NULL, shuffle = FALSE, showcoloursonly = FALSE)
#'
#' Returns a dictionary of colours for each class within each variable_list or phenotable
#' @export

make_colour_dictionary <- function(variable_list = NULL, pheno = NULL, phenolabels = NULL, class_to_ignore = "N_A", colour_of_class_to_ignore = "#bcc2c2", colour_table = NULL, shuffle = FALSE, showcoloursonly = FALSE){

    atoms <- c("0200FC", "ff7f02", "468B00", "CD2626", "267ACD", "640032", "FC00FA", "000000", "478559", "161748", "f95d9b", "39a0ca")

    JAMSpalette <- paste0("#", atoms)

    offsets <- c(2,4,6,8)
    for(offsetamplitude in offsets){
        offset <- as.hexmode(offsetamplitude)
        atomsoff <- as.character(as.hexmode(atoms) + offset)
        JAMSpalette_ext <- paste0("#", atomsoff)
        JAMSpalette <- c(JAMSpalette, JAMSpalette_ext)
    }

    if (shuffle == TRUE){
        JAMSpalette <- JAMSpalette[shuffle(JAMSpalette)]
    }

    if (showcoloursonly == TRUE){
        #make example pies if requested
        sliceValues <- as.data.frame(rep(10, length(atoms)))
        colnames(sliceValues) <- "Slice"
        colourpies <- list()
        for (pies in 1:(length(offsets) + 1)){
            currsliceValues <- sliceValues
            pallist <- split(1:length(JAMSpalette), ceiling(seq_along(1:length(JAMSpalette)) / length(atoms)))
            currpalette <- JAMSpalette[pallist[[pies]]]
            currsliceValues$Colour <- factor(currpalette, levels=as.character(currpalette))
            currcols <- as.character(currsliceValues$Colour)
            names(currcols) <- as.character(currsliceValues$Colour)
            colourpies[[pies]] <- ggplot(currsliceValues, aes(x = "", y = Slice, fill = Colour)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = currcols)  + coord_polar("y") + labs(title = paste0("JAMSpalette", pies))
        }
        return(colourpies)
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
            discrete2colour$Colour[(discrete2colour$Name %in% class_to_ignore)] <- colour_of_class_to_ignore
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

        return(cdict)
    }

}
