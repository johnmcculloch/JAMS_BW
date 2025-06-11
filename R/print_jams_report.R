#' print_jams_report(opt = opt)
#'
#' JAMSalpha function
#' @export

print_jams_report <- function(opt = opt, outputdir = NULL, elements = c("readplots", "taxonomy", "resfinder", "plasmidfinder", "func"), report_style = NULL){

    if (is.null(report_style)){
        report_style <- opt$analysis
    }

    QualPieColours <- c("#18bf06", "#037c91", "#f2b407", "#f2fa05", "#f71105", "#bcc2c2", "#bcc2c2")
    names(QualPieColours) <- c("HQ", "MHQ", "MQ", "LQ", "Contaminated", "N_A", "Unused")

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

    pdffn <- file.path(outputdir, paste(paste(opt$prefix, "JAMSalpha", "report", sep="_"), "pdf", sep = "."))
    pdf(pdffn, paper = "a4r")

    #Header with run info
    plot.new()
    info_to_print <- c("JAMS_version", "Run_type", "JAMS_Kdb_Version", "Sample_name", "Process", "Host_species", "Run_start", "Run_end", "Duration_hours", "CPUs_used", "Contigs_supplied", "GenBank_assembly_accession_number", "NCBI_SRA_accession_number", "Reads_source", "Read_chemistry", "Library_strategy", "Contig_assembler","Contig_assembler_version", "Input_number_bp_cap")
    info_to_print <- info_to_print[info_to_print %in% rownames(opt$projinfo)]

    grid.table(opt$projinfo[info_to_print, c("Run_info", "Run_value")], rows = NULL, cols = NULL, theme = ttheme_default(base_size = 12))

    plot.new()
    grid.table(authorshipmessage, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 15))

    #Print read and assembly plots, if available
    if ("readplots" %in% elementsinopt){
        print(opt$readplots)
    }

    #Print taxonomic analysis, if applicable
    if (("taxonomy" %in% elements)){
        taxonomic_spaces <- c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin")[c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin") %in% names(opt$abundances$taxonomic)]

        #Print an equivalent report for each available taxonomic space
        if (report_style %in% c("metagenome", "metatranscriptome")){
            for (curr_taxonomic_space in taxonomic_spaces){
                curr_qual_df <- opt$taxonomic_Quality_split_list[[curr_taxonomic_space]]
                #Account for unused contigs in space
                if (sum(curr_qual_df$Pct) < 100){
                    curr_qual_df["Unused", ] <- c("Unused", round((100 - sum(sum(curr_qual_df$Pct))), 2), 100)
                }
                curr_qual_df$Pct <- as.numeric(curr_qual_df$Pct)
                curr_qual_df$CumPct <- as.numeric(curr_qual_df$CumPct)
                curr_qual_df <- curr_qual_df[rev(curr_qual_df$Quality), ]
                curr_qual_df$Pos <- (curr_qual_df$Pct / 2) + lead(curr_qual_df$CumPct, 1)
                curr_qual_df$Pos[which(is.na(curr_qual_df$Pos))] <- curr_qual_df$Pct[which(is.na(curr_qual_df$Pos))] / 2
                curr_qual_df$Percentage <- paste0(curr_qual_df$Pct, "%")
                curr_qual_df$Quality <- factor(curr_qual_df$Quality, levels = curr_qual_df$Quality)

                QCplot <- ggplot(curr_qual_df, aes(x = "", y = Pct, fill = Quality)) + geom_bar(width = 1, stat = "identity")
                QCplot <- QCplot + scale_fill_manual(values = QualPieColours, breaks = names(QualPieColours), name = "Genome Quality")
                QCplot <- QCplot + geom_bar(stat = "identity", width = 1, colour = "white")
                QCplot <- QCplot + coord_polar("y", start = 0)
                QCplot <- QCplot + geom_label_repel(aes(y = Pos, label = Percentage), size = 4.5, nudge_x = 1, show.legend = FALSE) 
                QCplot <- QCplot + guides(fill = guide_legend(title = "Group"))
                QCplot <- QCplot + theme_void()
                QCplot <- QCplot + labs(title = paste("Quality of MAGs in", curr_taxonomic_space, "space"))
                QCplot <- QCplot + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

                print(QCplot)

                for (ordby in c("Percentage", "Completeness")){
                    print(pareto_abundance(LKTdose = opt$abundances$taxonomic[[curr_taxonomic_space]], max_num_taxa_to_show = 40, readstats = opt$readstats, read_denominator = "Assembled", orderby = ordby, samplename = opt$prefix, taxlevel = curr_taxonomic_space, showcompletenesstext = TRUE, showcompletenessline = TRUE, showcontaminationline = FALSE, showparetocurve = TRUE, showrelabundtext = TRUE))
                }

            }

            #Plot aggregated higher taxonomic levels
            higher_taxlvls_present <- rev(c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))[rev(c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) %in% colnames(opt$abundances$taxonomic$Contig_LKT)]
            for (taxlvls in higher_taxlvls_present){
                print(pareto_abundance(LKTdose = opt$abundances$taxonomic$Contig_LKT, max_num_taxa_to_show = 35, readstats = opt$readstats, read_denominator = "Assembled", orderby = "Percentage", samplename = opt$prefix, taxlevel = taxlvls, showcompletenesstext = FALSE, showcompletenessline = FALSE, showcontaminationline = FALSE, showparetocurve = TRUE, showrelabundtext = TRUE))
            }
        }
    }

    availblastanalyses <- elementsinopt[!(elementsinopt %in% c("readplots", "SixteenSid"))]

    batp <- NULL
    taxonomic_spaces_avail <- c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin")[c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin") %in% names(opt$abundances$functional)]
    for (ba in availblastanalyses){
        batp <- opt[[ba]]
        if (ba == "resfinder"){
             batp <- batp[,c("Accession", "Product", "Pident", "GeneNameNR", "Feature", taxonomic_spaces_avail)]
             colnames(batp) <- c("DB_Hit", "Function", "Percent_ID", "Gene_symbol", "Gene", taxonomic_spaces_avail)
        } else {
            batp <- batp[ , c("Accession", "Product", "Pident", taxonomic_spaces_avail)]
            colnames(batp) <- c("DB_Hit", "Function", "Percent_ID", taxonomic_spaces_avail)
        }

        #Shorten hit because it is too verbose from resfinder and vfdb
        #batp$DB_Hit <- gsub(":.*", "", batp$DB_Hit)
        ti <- switch(ba, "resfinder" = "Known antibiotic resistance genes", "vfdb" = "Virulence factors present in VFDB \n www.ncbi.nlm.nih.gov/pubmed/15608208", "plasmidfinder" = "Plasmid Replicon-associated genes present in PlasmidFinder \n https://www.ncbi.nlm.nih.gov/pubmed/24777092")
        print_table(tb = batp, tabletitle = ti, fontsize = 4, numrows = 25)
    }

    #if ("func" %in% elements){
    #    funcanalyses <- unique(opt$featuredose$Analysis)[!(unique(opt$featuredose$Analysis) %in% c(availblastanalyses, "napdos", "serofinderO", "serofinderH"))]
    #    for(fa in funcanalyses){
    #        plot_wordcloud_sample(opt=opt, analysis=fa, removeunclassifieds=TRUE)
    #    }
    #}
    plot.new()
    grid.table(authorshipmessage, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 15))
    dev.off()
}
