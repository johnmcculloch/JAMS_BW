#' plot_readdata(readdata=NULL)
#'
#' Plots available read data of sequencing statistics.
#' @export

plot_readdata <- function(readdata = NULL){

    #Make a copy of readdata, excluding samples with input as contigs
    if("Proj_type" %in% colnames(readdata)){
        dfr <- subset(readdata, Proj_type == "Assemble_from_reads")
    } else {
        dfr <- readdata
    }

    #Get appropriate size of font for x-axis
    relsize <- round(((-1/250) * (nrow(readdata))) + 1, 1)

    gvec <- NULL
    gvec <- vector("list", length = 4)
    g <- 1

    Bar.reads.QC <- ggplot(dfr, aes(Sample, PctQCpass)) + geom_bar(stat = "identity", fill = "blue") + labs(title = paste(opt$project, "Reads passing QC (%)", sep = " - "))
    Bar.reads.QC <- Bar.reads.QC + labs(x = "Sample", y = expression("Percentage bases passing QC"))
    Bar.reads.QC <- Bar.reads.QC + theme(axis.text.x = element_text(angle = 90, size = rel(relsize), colour = "black"))
    gvec[[g]] <- Bar.reads.QC
    g <- g + 1

    if(!(is.null(readdata$NonHost_bases))){
        Bar.reads.PctHost <- ggplot(readdata, aes(Sample, PctHost)) + geom_bar(stat = "identity", fill = "blue") + labs(title = paste(opt$project, "Bases belonging to Host Species (%)"))
        Bar.reads.PctHost <- Bar.reads.PctHost+labs(x = "Sample", y = expression("Percentage bases belonging to Host Species"))
        Bar.reads.PctHost <- Bar.reads.PctHost+theme(axis.text.x = element_text(angle = 90, size = rel(relsize), colour = "black"))
        gvec[[g]] <- Bar.reads.PctHost
        g <- g + 1

        dfr$GbNAHS <- round((dfr$NonHost_bases / 1000000000), 2)
        Bar.reads.NAHS <- ggplot(dfr, aes(Sample, GbNAHS)) + geom_bar(stat = "identity", fill = "blue") + labs(title = paste(opt$project, "Gigabytes of usable (non-host) data sequenced"))
        Bar.reads.NAHS <- Bar.reads.NAHS+labs(x = "Sample", y = expression("Non-Host Bases Sequenced (Gbp)"))
        Bar.reads.NAHS <- Bar.reads.NAHS+theme(axis.text.x = element_text(angle = 90, size = rel(relsize), colour = "black"))
        gvec[[g]] <- Bar.reads.NAHS
        g <- g + 1

    } else {

        dfr$GbTrim <- round((dfr$Trim_bases/1000000000),2)
        Bar.reads.trim <- ggplot(dfr, aes(Sample, GbTrim)) + geom_bar(stat = "identity", fill = "blue") + labs(title = paste(opt$project, "Gigabytes of QC-pass data sequenced"))
        Bar.reads.trim <- Bar.reads.trim + labs(x = "Sample", y = expression("QC-pass Bases Sequenced (Gbp)"))
        Bar.reads.trim <- Bar.reads.trim+theme(axis.text.x = element_text(angle = 90, size = rel(relsize), colour = "black"))
        gvec[[g]] <- Bar.reads.trim
        g <- g + 1
    }

    Bar.reads.PctAss <- ggplot(dfr, aes(Sample, PctAss)) + geom_bar(stat = "identity", fill = "blue") + labs(title = paste(opt$project, "Bases Assembled into Contigs (%)"))
    Bar.reads.PctAss <- Bar.reads.PctAss + labs(x = "Sample", y = expression("Percentage Bases Assembled into Contigs"))
    Bar.reads.PctAss <- Bar.reads.PctAss+theme(axis.text.x = element_text(angle = 90, size = rel(relsize), colour = "black"))
    gvec[[g]] <- Bar.reads.PctAss

    #Redefine graphics list as ones only containing plots
    gvec <- gvec[sapply(gvec, function(x){!(is.null(x))})]

    return(gvec)
}
