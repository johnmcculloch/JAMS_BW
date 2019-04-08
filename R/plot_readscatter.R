#' plot_readscatter(readdata=NULL, phenotable=NULL, discrete=NULL)
#'
#' Makes a scatterplot of percentage of reads assembled into contigs by sequencing depth.
#' @export

plot_readscatter<-function(readdata=NULL, phenotable=NULL, discrete=NULL){
    
    #Make a copy of readdata, excluding samples with input as contigs
    if("Proj_type" %in% colnames(readdata)){
        readdf<-subset(readdata, Proj_type == "Assemble_from_reads")
    } else {
        readdf<-readdata
    }

    readdf$Usable_Gb<-round(readdf$NonHost_bases/1000000000, 2)
    readdf$Type<-phenotable[match(readdf$Sample, phenotable[, which(colnames(phenotable)==(subset(phenolabels, Var_type=="Sample"))[]$Var_label)]), which(colnames(phenotable)==discrete)]
    assembly.plot<-ggplot(readdf, aes(x=Usable_Gb, y=PctAss)) + geom_point(aes(color=factor(Type))) + geom_smooth(method=loess, se=FALSE)
    assembly.plot<-assembly.plot+labs(color = discrete, title = "Contig Assembly Statistics")
    assembly.plot<-assembly.plot+xlab("Gigabases of assemblable data") + ylab("Percentage Bases assembled into Contigs")

    return(assembly.plot)
}
