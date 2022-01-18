#' pareto_abundance(LKTdose=NULL, samplename=NULL, taxlevel=NULL, assemblystats=NULL, showcompletenesstext=TRUE, showcompletenessline=TRUE, showcontigline=FALSE, showparetocurve=TRUE, showpercentagetext=TRUE)
#'
#' Plots a pareto chart of the relative abundance of a taxon.
#' @export

pareto_abundance <- function(LKTdose = NULL, samplename = NULL, list.data = NULL, taxlevel = NULL, assemblystats = NULL, showcompletenesstext = TRUE, showcompletenessline = TRUE, showcontigline = FALSE, showparetocurve = TRUE, showpercentagetext = TRUE){

    if (is.null(LKTdose)){
        if (is.null(list.data)){
            stop("You must provide an LKTdose or list.data object as input.")
        }
        if (is.null(samplename)){
            stop("Given a list.data object, you must provide the sample name for which you want to plot.")
        }
        LKTdose <- list.data[[paste(samplename, "LKTdose", sep="_")]]
        assemblystats <- list.data[[paste(samplename, "assemblystats", sep="_")]]
    }

    taxlvlspresent <- colnames(LKTdose)[colnames(LKTdose) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT")]

    #Make a colour dictionary by phylum
    tt<-LKTdose[, taxlvlspresent]
    tt[] <- lapply(tt, as.character)
    tt$Colour<-rainbow(length(tt$Phylum), start=.5, end=.1, alpha=1)[rank(tt$Phylum)]

    #Addpctfromctgs if applicable
    if("PctFromCtg" %in% colnames(LKTdose)){
        LKTdose$BasesFromCtg<-round(((LKTdose$PctFromCtg/100)*LKTdose$NumBases), 0)
    } else {
        LKTdose$BasesFromCtg<-0
    }

    #aggregate by desired taxlevel
    dfdose<-dplyr::group_by(LKTdose, get(taxlevel)) %>% dplyr::summarise(NumBases=sum(as.numeric(NumBases)), BasesFromCtg=sum(as.numeric(BasesFromCtg))) %>% dplyr::arrange(-NumBases)
    colnames(dfdose)<-c("Taxon", "NumBases", "BasesFromCtg")
    dfdose$PctFromCtg<-round((dfdose$BasesFromCtg/dfdose$NumBases)*100, 2)
    totbases<-sum(dfdose$NumBases)
    dfdose$Percentage<-round(((dfdose$NumBases/totbases)*100), 2)
    taxa<-unique(dfdose$Taxon)
    dfdose <- dfdose[order(dfdose$Percentage, decreasing=TRUE), ]
    dfdose$Pctcum <- cumsum(dfdose$Percentage)

    #Filter out points with 0.00%
    dfdose<-subset(dfdose, Percentage != 0.00)

    #Add assembly statistics, if applicable
    if(!(is.null(assemblystats))){
        assemblystatsinterest<-subset(assemblystats, TaxLevel==taxlevel)
        png<-assemblystatsinterest[, c("Taxon", "ProbNumGenomes")]
        n16s<-assemblystatsinterest[, c("Taxon", "Num16S")]
        dfdose<-left_join(dfdose, png)
        dfdose$ProbNumGenomes[is.na(dfdose$ProbNumGenomes)] <- 0.00
        dfdose<-left_join(dfdose, n16s)
        dfdose$Num16S[is.na(dfdose$Num16S)] <- 0
        dfdose$Gencomp<-paste0(dfdose$ProbNumGenomes, "(", dfdose$Num16S, ")")
    }

    if (nrow(dfdose)>30) {
        df.top <- dfdose[1:30, ]
        df.rest <- dfdose[31:(length(dfdose$Taxon)), ]
        if(!(is.null(assemblystats))){
            df.tail <- data.frame("Remainder", sum(df.rest$NumBases), sum(df.rest$BasesFromCtg),  round((sum(df.rest$BasesFromCtg) / sum(df.rest$NumBases))*100, 2),  sum(df.rest$Percentage), max(df.rest$Pctcum), sum(df.rest$ProbNumGenomes), sum(df.rest$Num16S), paste0(sum(df.rest$ProbNumGenomes), "(", sum(df.rest$Num16S), ")"))
        } else {
            df.tail <- data.frame("Remainder", sum(df.rest$NumBases), sum(df.rest$Percentage), max(df.rest$Pctcum))
        }
        colnames(df.tail) <- colnames(df.rest)
        dfdose <- rbind(df.top, df.tail)
    }

    #Add phylum information
    dfdose$PhylumInfo<-tt[ , "Phylum"][match(dfdose$Taxon, tt[ , taxlevel])]
    dfdose$PhylumInfo[which(dfdose$Taxon == "Remainder")]<-"Remainder"
    dfdose$Colour<-tt[ , "Colour"][match(dfdose$Taxon, tt[ , taxlevel])]
    dfdose$Colour[which(dfdose$Taxon == "Remainder")]<-"#000000"
    PhyColours<-dfdose$Colour
    names(PhyColours)<-dfdose$PhylumInfo
    #Set Genome completeness of remainder to 0
    dfdose[which(dfdose$PhylumInfo == "Remainder"), "Gencomp"]<-"0(0)"
    dfdose[which(dfdose$PhylumInfo == "Remainder"), "ProbNumGenomes"]<-as.numeric(0)

    pareto.title<-paste("Taxonomic relative abundance", taxlevel, sep = " - ")
    if(!(is.null(samplename))){
        pareto.title<-paste(pareto.title, "in", samplename)
    }

    if(showpercentagetext == TRUE){
        lbls<-paste(dfdose$Taxon, "=", dfdose$Percentage, "%", sep ="")
    } else {
        lbls<-dfdose$Taxon
    }
    #lbls<-paste0(dfdose$Taxon, " [", dfdose$PctFromCtg, "%ctgs]")
    dfdose$lbls <- lbls
    dfdose$lbls <- factor(dfdose$lbls, levels=dfdose$lbls)

    #Draw relabund bars
    if(taxlevel %in% c("Domain", "Kingdom", "Phylum")){
        plot.df <- ggplot(dfdose, aes(x=lbls)) + geom_bar(aes(y=Percentage), stat="identity", fill=dfdose$Colour)
    } else {
        plot.df <- ggplot(dfdose, aes(x=lbls)) + geom_bar(aes(y=Percentage, fill=PhylumInfo), stat="identity") +  scale_fill_manual(values = PhyColours)
    }

    #Add percent relabund up into a pareto curve
    if(showparetocurve == TRUE){
        plot.df <- plot.df + geom_point(aes(y=Pctcum)) + geom_path(aes(y=Pctcum, group=1))
    }

    #Add genome size and num of 16Ss on top of each bar if applicable.
    if(!(is.null(assemblystats))){
        if(showcompletenesstext == TRUE){
            plot.df <- plot.df + geom_text(data = dfdose, aes(x = lbls, y = Percentage, label = Gencomp), vjust = 0, angle = 90, nudge_y = 6.3, size = rel(1.8))
        }
        if(showcompletenessline == TRUE){
            plot.df <- plot.df + geom_path(aes(y = ProbNumGenomes*100, group=1, colour = "Expected Genome Completeness"))
            plot.df <- plot.df + labs(colour="Confidence stats")
        }
       #plot.df <- plot.df + scale_y_continuous(sec.axis = sec_axis(~., name = "Genome Completeness [%]"))
    }

    if(showcontigline == TRUE){
        plot.df <- plot.df + geom_path(aes(y = PctFromCtg, group=1, colour = "Percent ID from Contigs"))
        plot.df <- plot.df + labs(colour="Confidence stats")
    }

    #Add main and axis titles
    plot.df<-plot.df+ggtitle(pareto.title)
    plot.df<-plot.df+labs(x=taxlevel, y=expression("Percentage"))
    plot.df<-plot.df+theme(axis.text.x=element_text(angle=90, size = rel(0.9), colour="black"))

    return(plot.df)
}
