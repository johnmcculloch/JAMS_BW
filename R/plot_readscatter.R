#' plot_readscatter(opt = NULL)
#'
#' Makes a scatterplot of percentage of reads assembled into contigs by sequencing depth.
#' @export

plot_readscatter <- function(opt = NULL){

    readdata <- opt$readdata
    phenotable <- opt$phenotable
    phenolabels <- opt$phenolabels
    readdata <- readdata[opt$phenotable$Sample, ]
    varlist <- define_kinds_of_variables(phenolabels = opt$phenolabels, phenotable = opt$phenotable, maxclass = opt$maxnumclasses, maxsubclass = 10000, class_to_ignore = opt$class_to_ignore, verbose = FALSE)

    #Make a copy of readdata, excluding samples with input as contigs
    if ("Proj_type" %in% colnames(readdata)){
        readdf <- subset(readdata, Proj_type == "Assemble_from_reads")
    } else {
        readdf <- readdata
    }

    gvec <- NULL
    gvec <- vector("list", length = length(varlist$discrete))
    g = 1

    for (discrete in varlist$discrete){

        readdf$Usable_Gb <- round(readdf$NonHost_bases / 1000000000, 2)
        readdf$Type <- phenotable[match(readdf$Sample, phenotable[ , "Sample"]), which(colnames(phenotable) == discrete)]

        assembly.plot <- ggplot(readdf, aes(x = Usable_Gb, y = PctAss)) + geom_point(aes(color = factor(Type))) + geom_smooth(method = "loess", se = FALSE, span = 2)
        assembly.plot <- assembly.plot + labs(color = discrete, title = paste(opt$project, "Contig Assembly Statistics"))
        assembly.plot <- assembly.plot + xlab("Gigabases of assemblable data") + ylab("Percentage Bases assembled into Contigs")

        gvec[[g]] <- assembly.plot
        g <- g + 1
    }

    #Redefine graphics list as ones only containing plots
    gvec <- gvec[sapply(gvec, function(x){ !(is.null(x)) })]

    return(gvec)
}
