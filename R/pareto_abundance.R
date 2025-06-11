#' pareto_abundance(LKTdose = NULL, max_num_taxa_to_show = 30, orderby = "Percentage", readstats = NULL, read_denominator = "Assembled", samplename = NULL, taxlevel = NULL, showrelabundtext = TRUE, showcompletenesstext = TRUE, showcompletenessline = TRUE, showcontaminationline = FALSE, showparetocurve = TRUE, rotate_plot = FALSE)
#'
#' Plots a pareto chart of the relative abundance of a taxon.
#' @export

pareto_abundance <- function(LKTdose = NULL, max_num_taxa_to_show = 30, orderby = "Percentage", readstats = NULL, read_denominator = "Assembled", samplename = NULL, taxlevel = NULL, showrelabundtext = TRUE, showcompletenesstext = TRUE, showcompletenessline = TRUE, showcontaminationline = FALSE, showparetocurve = TRUE, rotate_plot = FALSE){

    #Determine sequencing effort in basepairs
    if (!is.null(readstats)){
        totbases <- as.numeric(readstats[read_denominator, "Base_counts"])
    } else {
        totbases <- sum(LKTdose$NumBases)
    }

    taxlvlspresent <- colnames(LKTdose)[colnames(LKTdose) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT", "Contig_LKT", "MB2bin", "ConsolidatedGenomeBin")]

    #Make a colour dictionary by phylum
    tt <- LKTdose[, taxlvlspresent]
    tt[] <- lapply(tt, as.character)
    tt$HexCol <- rainbow(length(tt$Phylum), start=.5, end=.1, alpha=1)[rank(tt$Phylum)]

    #Make a colour dictionary by Quality
    Qual_colour_dict <- data.frame(Quality = c("HQ", "MHQ", "MQ", "LQ", "Contaminated", "N_A"), HexCol = c("#18bf06", "#037c91", "#f2b407", "#f2fa05", "#f71105", "#bcc2c2"))
    rownames(Qual_colour_dict) <- Qual_colour_dict$Quality

    LKTdose <- LKTdose[ , c(taxlvlspresent, "NumBases", "Completeness", "Contamination", "Quality", "RefScore")]
    #Ensure dark matter is marked as such
    LKTdose[which(LKTdose$LKT == "LKT__Unclassified"), c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1")] <- c("d__Unclassified", "k__Unclassified", "p__Unclassified", "c__Unclassified", "o__Unclassified", "f__Unclassified", "g__Unclassified", "s__Unclassified", "is1__Unclassified")[c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1") %in% colnames(LKTdose)]

    #aggregate by desired taxlevel
    if (!(taxlevel %in% c("LKT", "Contig_LKT", "MB2bin", "ConsolidatedGenomeBin"))){
        dfdose <- LKTdose %>% group_by_at(taxlevel) %>% dplyr::summarise(NumBases = sum(as.numeric(NumBases)))
        dfdose[ , c("Completeness", "Contamination")] <- 0
        dfdose[ , c("Quality", "RefScore")] <- "N_A"
        #default to ordering by relabund
        orderby <- "Percentage"
        colourby <- "Phylum"
        showcompletenesstext <- FALSE
        showcompletenessline <- FALSE
        showcontaminationline <- FALSE
    } else {
        dfdose <- LKTdose
        colourby <- "Quality"
    }

    dfdose <- as.data.frame(dfdose)
    dfdose$Percentage <- round(((dfdose$NumBases / totbases) * 100), 2)
    dfdose$Taxon <- dfdose[ , taxlevel]

    #Filter out points with 0.00%
    dfdose <- subset(dfdose, Percentage != 0.00)
    #Order by requested criterion
    dfdose <- dfdose[order(dfdose[ , orderby], decreasing = TRUE), ]
    dfdose$Pctcum <- cumsum(dfdose$Percentage)
    dfdose <- dfdose[ , c("Taxon", "NumBases", "Percentage", "Pctcum", "Completeness", "Contamination", "Quality", "RefScore")]

    if (nrow(dfdose) > max_num_taxa_to_show) {
        df.top <- dfdose[1:max_num_taxa_to_show, ]
        df.rest <- dfdose[(max_num_taxa_to_show + 1):(length(dfdose$Taxon)), ]
        df.tail <- data.frame(Taxon = "Remainder", NumBases = sum(df.rest$NumBases), Percentage = round(sum(df.rest$Percentage), 2), Pctcum = round(max(df.rest$Pctcum), 2), Completeness = 0, Contamination = 0, Quality = "N_A", RefScore = "N_A")
        colnames(df.tail) <- colnames(df.rest)
        rownames(df.tail) <- "Remainder"
        dfdose <- rbind(df.top, df.tail)
    }

    #Add bar colour information
    if (!(taxlevel %in% c("LKT", "Contig_LKT", "MB2bin", "ConsolidatedGenomeBin"))){
        if (!(taxlevel %in% c("Domain", "Kingdom", "Phylum"))){
            dfdose$ColourInfo <- tt[ , "Phylum"][match(dfdose$Taxon, tt[ , taxlevel])]
            dfdose$ColourInfo[which(dfdose$Taxon == "Remainder")] <- "Remainder"
            dfdose$Colour <- tt[ , "HexCol"][match(dfdose$ColourInfo, tt[ , "Phylum"])]
            dfdose$Colour[which(dfdose$Taxon == "Remainder")] <- "#000000"
            cd <- dfdose[ , c("ColourInfo", "Colour")]
            cd <- cd[!duplicated(cd$ColourInfo), ]
            BarColours <- cd$Colour
            names(BarColours) <- cd$ColourInfo
            BarName <- "Phylum"
        } else {
            dfdose$ColourInfo <- NULL
            dfdose$Colour <- NULL
        }
    } else {
        dfdose$ColourInfo <- dfdose$Quality
        dfdose$Colour <- Qual_colour_dict[ , "HexCol"][match(dfdose$Quality, Qual_colour_dict[ , "Quality"])]
        dfdose$Colour[which(dfdose$Taxon == "Remainder")] <- "#000000"
        BarColours <- Qual_colour_dict$HexCol
        names(BarColours) <- Qual_colour_dict$Quality
        BarName <- "Quality"
    }

    #Fill in NAs
    dfdose$Completeness[which(is.na(dfdose$Completeness))] <- 0
    dfdose$Contamination[which(is.na(dfdose$Contamination))] <- 0
    dfdose$RefScore[which(is.na(dfdose$RefScore))] <- "N_A"

    pareto.title <- paste("Taxonomic relative abundance", taxlevel, sep = " - ")
    if (!(is.null(samplename))){
        pareto.title <- paste0(pareto.title, "\n", samplename)
    }

    if (showcompletenesstext == TRUE){
        lbls <- paste0(dfdose$Taxon, "; ", dfdose$Completeness, "% GenComp")
    } else {
        lbls <- dfdose$Taxon
    }

    dfdose$lbls <- lbls
    dfdose$lbls <- factor(dfdose$lbls, levels = dfdose$lbls)

    #Draw relabund bars
    if (!(taxlevel %in% c("Domain", "Kingdom", "Phylum"))){
        plot.df <- ggplot(dfdose, aes(x = lbls)) + geom_bar(aes(y = Percentage, fill = ColourInfo), stat = "identity") +  scale_fill_manual(values = BarColours, breaks = names(BarColours), name = BarName)
    } else {
        plot.df <- ggplot(dfdose, aes(x = lbls)) + geom_bar(aes(y = Percentage), fill = "blue", stat = "identity")
    }

    fontsize_x <- compute_x_font(nx = nrow(dfdose), upper_n = 40, upper_fs = 5, lower_n = 5, lower_fs = 12, cex = 0.72)

    #Add percent relabund up into a pareto curve
    if (showparetocurve == TRUE){
        plot.df <- plot.df + geom_point(aes(y = Pctcum)) + geom_path(aes(y = Pctcum, group = 1))
    }

    if (showrelabundtext){
        plot.df <- plot.df + geom_text(data = dfdose, aes(x = lbls, y = Percentage, label = Percentage), vjust = 0, angle = 90, nudge_y = 6.2, size = rel(1.75))
    }

    if (showcompletenessline){
        plot.df <- plot.df + geom_path(aes(y = Completeness, group = 1, colour = "Genome Completeness"))
        plot.df <- plot.df + labs(colour = "Confidence stats")
    }

    if (showcontaminationline){
        plot.df <- plot.df + geom_path(aes(y = Contamination, group = 1, colour = "Contamination"))
        plot.df <- plot.df + labs(colour = "Confidence stats")
    }

    #Add main and axis titles
    plot.df <- plot.df + ggtitle(pareto.title)
    plot.df <- plot.df + labs(x = taxlevel, y = expression("Percentage"))
    plot.df <- plot.df + theme(plot.title = element_text(size = 12, hjust = 0.5, colour = "black") , axis.text.x = element_text(angle = 90, size = fontsize_x, colour = "black"))

    if (rotate_plot){
        #plot.df <- plot.df + coord_flip()
    }

    return(plot.df)
}
