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

    if (is.null(variable_list)){
        if (is.null(phenolabels)){
            flog.info("Imputing types of columns on the metadata.")
            Var_label <- colnames(pheno)
            Var_type <- sapply(Var_label, function (x) { infer_column_type(phenotable = pheno, colm = x, class_to_ignore = class_to_ignore) } )
            phenolabels <- data.frame(Var_label = unname(Var_label), Var_type = unname(Var_type), stringsAsFactors = FALSE)
        }

        #Define kinds of variables
        variable_list <- define_kinds_of_variables(phenolabels = phenolabels, phenotable = pheno, maxclass = 30, maxsubclass = 100, class_to_ignore = class_to_ignore, verbose = FALSE)
    }

    #See how many different things there are
    discretes <- unique(unname(unlist(variable_list[which(!(names(variable_list) %in% c("sample", "continuous")))])))
    uniqueclasses <- unique(unlist(pheno[ , discretes]))

    #Attribute colours to a non-redundant list of classes
    #Start by attributing colours to things that are named the same

    if (length(uniqueclasses) > length(JAMSpalette)){
        flog.warn(paste0("The number of unique classes on metadata (", length(uniqueclasses), ") is greater than the number of colors in the standard JAMS color palette. JAMS will recycle colors for this colour dictionary, but you may want to ponder over the magnitude of your metadata..."))

        JAMSpalette <- base::rep(JAMSpalette, times = ceiling(length(uniqueclasses) / length(JAMSpalette)))
    }
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
        ctab <- subset(discrete2colour, Name %in% classnames)
        ctab$Hex <- ctab$Colour
        cdict[[d]] <- ctab

        names(cdict)[d] <- discretes[d]
    }

    ctable_new <- plyr::rbind.fill(cdict)
    ctable_new <- ctable_new[!duplicated(ctable_new$Name), ]
    rownames(ctable_new) <- ctable_new$Name


    if (FALSE){
        pdf("JAMSpalette_values.pdf", paper = "a4r")

        show_ncol <- 4 %||% ceiling(sqrt(nrow(ctable_new)))
        show_nrow <- ceiling(nrow(ctable_new) / show_ncol)

        show_colours <- c(ctable_new$Colour, rep(NA, show_nrow * show_ncol - nrow(ctable_new)))
        show_colours_mat <- matrix(show_colours, ncol = show_ncol, byrow = TRUE)
        show_labels <- c(ctable_new$Name, rep(NA, show_nrow * show_ncol - nrow(ctable_new)))
        show_labels <- unname(sapply(show_labels, function(x){ split_featname(featname = x, thresh_featname_split = 20)}))
        show_labels_mat <- matrix(show_labels, ncol = show_ncol, byrow = TRUE)
        hcl <- farver::decode_colour(show_colours, "rgb", "hcl")
        label_col <- ifelse(hcl[, "l"] > 50, "black", "white")

        show_colours_df <- data.frame(fill = show_colours, colour = label_col, label = show_labels)

        gl <- mapply(function(f,l,c) grobTree(rectGrob(gp=gpar(fill=f, col="white",lwd=2)), textGrob(l, gp=gpar(col=c))), f = show_colours_df$fill, l = show_colours_df$label, c = show_colours_df$colour, SIMPLIFY = FALSE)

        print(grid.arrange(grobs=gl, nrow = NULL, ncol = NULL))


        dev.off()
    }

    if (object_to_return == "ctable"){
        return(ctable_new)
    } else {
        return(cdict)
    }
}
