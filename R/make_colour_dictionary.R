#' make_colour_dictionary(variable_list=NULL, pheno=NULL)
#'
#' Returns a dictionary of colours for each class within each variable_list
#' @export

make_colour_dictionary <- function(variable_list = NULL, pheno = NULL, class_to_ignore = NULL, colour_of_class_to_ignore = "#bcc2c2", shuffle = FALSE){

    #Build a palette
    #atoms <- c("07437A", "96CABB", "FBEBA8", "FAB370", "ED7E65", "41303E", "486746", "203C51", "AC1629", "4E4346", "2D734E", "D6D77C", "2AD1BF", "FBFAF4", "1E1A11", "C48556", "38433F", "06FBFA", "E15631")

    #atoms <- c("ce4bd6", "23b6cc", "308ceb", "f7bf59", "491544", "cc9e00", "ff5575", "005086",  "277647",  "7a83ff", "818181", "42d98f",  "676a78", "9f5c19")

    atoms <- c("0075f2", "096b72", "ff7300", "00f2f2",  "fc03f0", "e0d12b", "fc0303", "0df505", "000000")

    JAMSpalette <- paste0("#", atoms)

    offsets <- c(48, -71)

    for(offsetamplitude in offsets){
        offset <- as.hexmode(offsetamplitude)
        atomsoff <- as.character(as.hexmode(atoms) + offset)
        JAMSpalette_ext <- paste0("#", atomsoff)
        JAMSpalette <- c(JAMSpalette, JAMSpalette_ext)
    }

    if (shuffle == TRUE){
        JAMSpalette <- JAMSpalette[shuffle(JAMSpalette)]
    }

    #See how many different things there are
    discretes <- unique(unname(unlist(variable_list[which(!(names(variable_list) %in% c("sample", "continuous")))])))
    uniqueclasses <- unique(unlist(pheno[ , discretes]))

    #Attribute colours to a non-redundant list of classes
    if (length(uniqueclasses) > length(JAMSpalette)){
        cat("There are far too many different classes in the metadata. Humans cannot possibly distinguigh this number of colours.\n")

        return(NULL)

    } else {
        #Start by attributing colours to things that are named the same
        discrete2colour <- data.frame(Name = uniqueclasses, Colour = JAMSpalette[1:length(uniqueclasses)], stringsAsFactors = FALSE)

        #Make class_to_ignore black
        if (!is.null(class_to_ignore)){
            discrete2colour$Colour[(discrete2colour$Name %in% class_to_ignore)] <- colour_of_class_to_ignore
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
