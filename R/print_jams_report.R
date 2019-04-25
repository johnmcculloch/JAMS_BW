#' print_jams_report(opt=opt)
#'
#' JAMSalpha function
#' @export

print_jams_report <- function(opt=opt, outputdir=NULL, elements=c("readplots", "taxonomy", "SixteenSid", "resfinder", "plasmidfinder", "vfdb", "func"), plotucobias=FALSE){

    #See what is available
    elementsinopt <- elements[elements %in% names(opt)]
    sampletype <- as.character(opt$projinfo[which(opt$projinfo$Run_info == "Run_type"), ]$Run_value)
    samplename <- as.character(opt$projinfo[which(opt$projinfo$Run_info == "Sample_name"), ]$Run_value)

    #Authorship message
    authors <- as.character(as.person(packageDescription("JAMS")$Author))
    authorshipmessage <- c(packageDescription("JAMS")$Title, paste("JAMS version", packageVersion("JAMS")), authors, paste("Contact:", "john.mcculloch@nih.gov"), "National Cancer Institute", "National Institutes of Health", "Bethesda, MD, USA")

    if (missing(outputdir)){
        outputdir <- opt$outdir
    }

    pdffn <- file.path(outputdir, paste(paste(opt$prefix, "JAMSalpha", "report", sep="_"), "pdf", sep="."))
    pdf(pdffn, paper="a4r")

    #Header with run info
    plot.new()
    grid.table(opt$projinfo, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 12))

    plot.new()
    grid.table(authorshipmessage, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 15))

    #Print read and assembly plots, if available
    if ("readplots" %in% elementsinopt){
        print(opt$readplots)
    }

    #Print taxonomic analysis, if applicable
    if (("taxonomy" %in% elements)){
        taxlvlspresent <- colnames(opt$LKTdose)[colnames(opt$LKTdose) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT")]
        for (taxl in taxlvlspresent){
            print(pareto_abundance(LKTdose=opt$LKTdose, samplename=samplename, taxlevel=taxl, assemblystats=opt$assemblystats, showcompletenesstext=TRUE, showcompletenessline=FALSE, showcontigline=FALSE, showparetocurve=TRUE, showpercentagetext=TRUE))
        }
        #Print an extra pareto chart for LKT with genome completeness
        print(pareto_abundance(LKTdose=opt$LKTdose, assemblystats=opt$assemblystats, samplename=samplename, taxlevel="LKT", showcompletenesstext=TRUE, showcompletenessline=TRUE, showcontigline=FALSE, showparetocurve=TRUE, showpercentagetext=TRUE))
        if ("ucoplots" %in% names(opt)){
            print(opt$ucoplots)
        }
    }

    #If transcriptome, report proportion of CDS to rRNA
    #if(opt$analysis %in% c("metatranscriptome", "isolaternaseq")){
        RNArelabund <- subset(opt$featuredose, Analysis == "FeatType")
        RNArelabund <- RNArelabund[, c("Accession", "Description", "NumBases")]
        totbases <- sum(RNArelabund$NumBases)
        RNArelabund$PPM <- round((RNArelabund$NumBases/totbases)*1000000, 0)
        dfRNA <- RNArelabund[,c("Accession", "PPM")]
        colnames(dfRNA) <- c("FeatType", "PPM")
        dfRNA$Pct <- round(dfRNA$PPM/10000, 1)
        dfRNA$lbls <- paste0(dfRNA$FeatType, " = ", dfRNA$Pct, "%")
        dfRNA$lbls <- as.factor(dfRNA$lbls)
        dfRNA <- subset(dfRNA, Pct > 0)
        Featcolours <- c("green", "red", "blue", "yellow", "brown", "black")[1:nrow(dfRNA)]
        rnaP <- ggplot(dfRNA, aes(x = "", y = Pct, fill = lbls)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = Featcolours) + coord_polar("y") + labs(title = paste(opt$prefix, "% Feature type found in Contigs", sep=" - "))
        rnaP <- rnaP + labs(fill = "Feature Type")
        rnaP <- rnaP + theme( axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks <- element_blank(), plot.title=element_text(size=14, face="bold"))
        print(rnaP)
    #}

    #Print consolidated 16S table, if available
    if (("SixteenSid" %in% colnames(opt$featuredata)) && ("SixteenSid" %in% elements)){
       taxa_16S_cons <- subset(opt$featuredata, SixteenSid != "none")
       if (nrow(taxa_16S_cons) > 0){
           taxa_16S_cons <- taxa_16S_cons[, c("Feature", "LengthDNA", "Contig", "SixteenSid")]
           taxa_16S_cons <- left_join(taxa_16S_cons, opt$contigsdata)
           taxa_16S_cons <- taxa_16S_cons[, c("Feature", "LengthDNA", "Contig", "Length", "SixteenSid", "LKT")]
           colnames(taxa_16S_cons) <- c("Gene", "Length_16S", "Contig", "Contig_Length", "LKT_dada2", "LKT_kraken")
           ti <- "Consolidated 16S rRNA table"
           print_table(tb=taxa_16S_cons, tabletitle=ti, fontsize=5, numrows=20)
       } else {
           print("No 16S rRNA features found in contigs")
       }
    }

    availblastanalyses <- elementsinopt[!(elementsinopt %in% c("readplots", "SixteenSid"))]

    batp <- NULL
    for (ba in availblastanalyses){
        batp <- opt[[ba]]
        if(ba == "resfinder"){
            contig2length <- opt$contigsdata[, c("Contig","Length")]
            batp <- batp[ , c("Accession","Class","Pident","Contig","LKT")]
            batp <- left_join(batp, contig2length, by="Contig")
            colnames(batp) <- c("DB_Hit","Function","Percent_ID","Contig","LKT", "Length")
            batp <- batp[,c("DB_Hit","Function","Percent_ID","Contig","Length","LKT")]
        } else {
            batp <- batp[ , c("Accession","Product","Pident","LKT")]
            colnames(batp) <- c("DB_Hit","Function","Percent_ID","LKT")
        }

        #Shorten hit because it is too verbose from resfinder and vfdb
        batp$DB_Hit <- gsub(":.*", "", batp$DB_Hit)
        ti <- switch(ba, "resfinder" = "Known antibiotic resistance genes", "vfdb" = "Virulence factors present in VFDB \n www.ncbi.nlm.nih.gov/pubmed/15608208", "plasmidfinder" = "Plasmid Replicon-associated genes present in PlasmidFinder \n https://www.ncbi.nlm.nih.gov/pubmed/24777092")
        print_table(tb=batp, tabletitle=ti, fontsize=6, numrows=20)
    }

    if("func" %in% elements){
        funcanalyses <- unique(opt$featuredose$Analysis)[!(unique(opt$featuredose$Analysis) %in% c(availblastanalyses, "napdos", "serofinderO", "serofinderH"))]
        for(fa in funcanalyses){
            plot_wordcloud_sample(opt=opt, analysis=fa, removeunclassifieds=TRUE)
        }
    }

    plot.new()
    grid.table(authorshipmessage, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 15))
    dev.off()
}
