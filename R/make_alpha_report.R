#' make_alpha_report(project=NULL, expvec=NULL, usefulexp=NULL, variable_list=NULL, measures=c("Observed", "InvSimpson"), cdict=NULL)
#'
#' Generates standard comparative plots for alpha diversity measures.
#' @export

make_alpha_report <- function(project = NULL, expvec = NULL, usefulexp = NULL, appendtofilename = NULL, variable_list = NULL, measures = c("Observed", "InvSimpson", "GeneCounts"), cdict = NULL, stratify_by_kingdoms = TRUE, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, addtit = NULL, max_pairwise_cats = 4, ignoreunclassified = TRUE, class_to_ignore = "N_A", ...){

    if (is.null(usefulexp)){
        usefulexp <- names(expvec)[!(names(expvec) %in% c("FeatType", "vfdb", "SFLD", "Coils", "Gene3D", "Phobius", "ProSitePatterns", "SMART", "resfinder", "ProDom"))]
    } else {
        usefulexp <- names(expvec)[(names(expvec) %in% usefulexp)]
    }
    expvec2 <- expvec[usefulexp]

    #if (is.null(applyfilters)){
    #    applyfilters <- "none"
    #}

    #Set variables
    pt <-as.data.frame(colData(expvec2[[1]]))
    readscols <- (colnames(pt)[(colnames(pt) %in% c("GbNAHS", "PctAss"))])

    variables_disc <- variable_list$discrete
    variables_bin <- variable_list$binary
    variables_cont <- variable_list$continuous
    variables_subs <- variable_list$subsettable
    variables_all <- unique(c(variables_disc, variables_subs, variables_bin, variables_cont))
    variables_discont <- variables_all[!(variables_all %in% variables_cont)]
    variables_cont <- append(variables_cont, readscols, after = length(variables_all))

    pt <- as.data.frame(colData(expvec2[[1]]))

    if  (!is.null(class_to_ignore)){
        flog.info(paste("Metadata value(s)", paste0(class_to_ignore, collapse = ", "), "will be omitted." ))
    }

    if (is.null(project)){
        project <- "MyProject"
    }

    Pn <- 0
    plotname <- "Alpha_Diversity"
    comparisons <- variables_all

    basepdffn <- paste("JAMS", project, plotname, sep="_")
    if (!(is.null(appendtofilename))){
        basepdffn <- paste(basepdffn, appendtofilename, sep = "_")
    }
    pdffn <- paste(basepdffn, "pdf", sep=".")
    #Check if fn exists to avoid overwriting.
    ffn=1
    while(file.exists(pdffn)){
        pdffn <- paste(paste("JAMS", project, plotname, ffn, sep = "_"), "pdf", sep = ".")
        ffn = ffn + 1
        if (ffn > 100){
            stop("I think you have enough plots already. May I suggest you try looking through them.")
        }
    }

    pdf(pdffn, paper = "a4r")
    authors = as.character(as.person(packageDescription("JAMS")$Author))
    intromessage = c(packageDescription("JAMS")$Title, paste("JAMS version", packageVersion("JAMS")), authors, paste("Contact:", "john.mcculloch@nih.gov"),"National Cancer Institute", "National Institutes of Health", "Bethesda, MD, USA")
    plot.new()
    grid.table(intromessage, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 15))
    plot.new()
    grid.table(c("JAMS", plotname, project), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

    for(cmp in comparisons){

        for(a in 1:length(expvec2)){

            plot.new()
            grid.table(c(names(expvec2)[a], plotname, paste("between", cmp), "No subsetting"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

            print(plot_alpha_diversity(ExpObj = expvec2[[a]], measures = measures, stratify_by_kingdoms = stratify_by_kingdoms, glomby = glomby, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, subsetby = NULL, compareby = cmp, colourby = NULL, shapeby = NULL, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, cdict = cdict, addtit = addtit, signiflabel = "p.format", max_pairwise_cats = max_pairwise_cats, ignoreunclassified = ignoreunclassified, class_to_ignore = class_to_ignore, ...))

            if (length(variables_subs) > 0){
                #If there is any subsettable data, then do that.
                validsubs <- NULL
                validsubs <- variables_subs[!(variables_subs %in% cmp)]
                for (vs in validsubs){
                    plot.new()
                    grid.table(c(names(expvec2)[a], plotname, paste("between", cmp), paste("Subset by", vs)), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                    print(plot_alpha_diversity(ExpObj = expvec2[[a]], measures = measures, stratify_by_kingdoms = stratify_by_kingdoms, glomby = glomby, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, subsetby = vs, compareby = cmp, colourby = NULL, shapeby = NULL, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, cdict = cdict, addtit = addtit, signiflabel = "p.format", max_pairwise_cats = max_pairwise_cats, ignoreunclassified = ignoreunclassified, class_to_ignore = class_to_ignore, ...))
                }
            }
        }
    }

    dev.off()

    flog.info(paste(plotname, "report complete"))
}
