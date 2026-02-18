#' Plots relative abundance heatmaps annotated by the metadata using as input a SummarizedExperiment object

#' @param ExpObj JAMS-style SummarizedExperiment object

#' @param glomby String giving the taxonomic level at which to agglomerate counts. This argument should only be used with taxonomic SummarizedExperiment objects. When NULL (the default), there is no agglomeration

#' @param hmtype Type of heatmap to plot. Options are "exploratory" or "comparative". When "exploratory" is passed, features are ranked by variance across all samples. When "comparative" is passed, the metadata variable name for grouping samples must be passed using the compareby argument (see below).

#' @param compareby String specifying the metadata variable name for grouping samples when the hmtype argument is set to "comparative". This will calculate p-values for each feature using the Mann-Whitney-Wilcoxon U-test when there are exactly two classes within the variable, and the log2 foldchange between the two groups will be calculated. When there are three or more classes within the variable, the p-value will be calculated using ANOVA. If there is only a single class within the variable, hmtype will default to "exploratory" and features will be ranked by variance across samples.

#' @param wilcox_paired_by String specifying metadata variable name for pairing samples if compareby has exactly two classes. If passed and not NULL, the Mann-Whitney-Wilcoxon U-test will be paired with samples paired by duplicate values in the wilcox_paired_by column. In this case, a test will be made to ascertain whether there are exactly two samples, one in each compareby class, bearing an identical value for wilcox_paired_by. If there is an odd number of samples, or somehow they cannot be unambiguously paired, the Mann-Whitney-Wilcoxon test will be run in the paired=FALSE mode, with a warning message.

#' @param colcategories Vector with variables to include on the header sample annotation. A key legend for each variable will be included. If set to NULL, all variables contained in the SummarizedExperiment object metadata containing between 2 and 10 classes will be included. Variables containing continuous data will be plot with a gradient scale.

#' @param samplesToKeep Vector with sample names to keep. If NULL, all samples within the SummarizedExperiment object are kept. Default is NULL.

#' @param featuresToKeep Vector with feature names to keep. If NULL, all features within the SummarizedExperiment object are kept. Default is NULL. Please note that when agglomerating features with the glomby argument (see above), feature names passed to featuresToKeep must be post-agglomeration feature names. For example, if glomby="Family", featuresToKeep must be family names, such as "f__Enterobacteriaceae", etc.

#' @param subsetby String specifying the metadata variable name for subsetting samples. If passed, multiple plots will be drawn, one plot for samples within each different class contained within the variable.  If NULL, data is not subset. Default is NULL.

#' @param invertbinaryorder Requires a logical value. If set to TRUE, when using compareby when there are exactly two classes within the variable, the log2 foldchange signs of the two groups will be inverted. Default is FALSE.

#' @param hmasPA Requires a logical value. If set to TRUE, cells within a heatmap will be plot as present (red) or absent (black) instead of using a continuous colour scale. Default is FALSE.

#' @param threshPA Numeric value setting the threshold for absence in Presence/Absence heatmaps. Default is 0.

#' @param ntop Numeric value setting the cap for the maxiumum number of features to plot. If NULL, all features surviving filtration settings will be plot. This argument is usually used when plotting "exploratory" heatmaps. For "comparative" heatmaps, it is often more practical to cap the number of features show on a heatmap with showonlypbelow (see below). Default is NULL.

#' @param splitcolsby String specifying the metadata variable for splitting the heatmap into column groups. Samples are thus grouped by each existing class within the variable. Default is NULL.

#' @param ordercolsby String specifying the metadata variable for ordering the samples sequentially. If variable is continuous, order will be from smaller to greater values; if variable is discrete, order will be alphabetical. If NULL, samples will be clustered naively and a dendrogram generated. Default is NULL.

#' @param column_split_group_order When used with splitcolsby, column_split_group_order will take a vector setting the explicit order of classes by which to order the split column groups. For clustering samples within each split column group see cluster_column_slices.

#' @param cluster_column_slices Requires a logical value. If set to TRUE and columns are split into column groups (see splitcolsby), clustering within the slice is performed. Default is TRUE. To turn off clustering, thereby set to FALSE.

#' @param cluster_samples_per_heatmap Requires a logical value. If set to TRUE, samples will be clustered within each heatmap using only the features shown on that single heatmap. If set to FALSE, pre-clusterization of all samples will be performed without drawing a heatmap so that the sample order can then be maintained on all heatmaps until all features surviving filtration have been plot. The number of features per heatmap can be set by max_rows_in_heatmap. Please note that if set to FALSE, plotting of heatmaps may be substantially slower as clusterization of the master matrix will be computed first. Default is TRUE.

#' @param cluster_rows Requires a logical value. If set to TRUE, features will be clustered together on the heatmap based on sample structure. See cluster_features_per_heatmap for more options. Default is TRUE.

#' @param cluster_features_per_heatmap Requires a logical value. If set to TRUE, features will be clustered within each heatmap using only the samples shown on that single heatmap. If set to FALSE, pre-clusterization of all features surviving filtration will be performed without drawing a heatmap as to define feature order on all heatmaps until all features surviving filtration have been plot. The number of features per heatmap can be set by max_rows_in_heatmap. Please note that if set to FALSE, plotting of heatmaps may be substantially slower as clusterization of the master matrix will be computed first. Default is TRUE.

#' @param textby String specifying the metadata variable containing classes for annotating samples with text on the top of each column. Default is NULL. Please note that label_samples must be set to TRUE for metadata text to be added.

#' @param label_samples Requires a logical value. If set to TRUE (the default), the sample names will be printed at the bottom of each sample column.

#' @param max_rows_in_heatmap Numeric value setting the maximum number of rows to be plot on a single heatmap. Default is 50. If there are more features surviving filtration than this number, they will be plot in subsequent heatmaps until all features have been plot. See ntop for capping the number of features. See maxnumheatmaps for capping the number of heatmaps to be plot in a comparison.

#' @param maxnumheatmaps Numeric value setting the maximum number of heatmaps to be plot within a comparison. If there is a larger number of features surviving filtration than whatever value is set by max_rows_in_heatmap, multiple heatmaps will be generated untill all features have been plot. When maxnumheatmaps is set to NULL (the default) there is no cap on the number of heatmaps to be generated within a comparison. See also ntop and max_rows_in_heatmap.

#' @param featcutoff Requires a numeric vector of length 2 for specifying how to filter out features by relative abundance. The first value of the vector specifies the minimum relative abundance in Parts per Million (PPM) and the second value is the percentage of samples which must have at least that relative abundance. Thus, passing c(250, 10) to featcutoff would filter out any feature which does not have at least 250 PPM (= 0.025 percent) of relative abundance in at least 10 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the featcutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param GenomeCompletenessCutoff Requires a numeric vector of length 2 for specifying how to filter out features by genome completeness. This is, of course, only applicble for taxonomic shotgun SummarizedExperiment objects. When passed to non-taxonomic shotgun SummarizedExperiment objects, GenomeCompletenessCutoff will be ignored. The first value of the vector specifies the minimum genome completeness in percentage  and the second value is the percentage of samples which must have at least that genome completeness. Thus, passing c(50, 5) to GenomeCompletenessCutoff would filter out any taxonomic feature which does not have at least 50 percent of genome completeness in at least 5 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the GenomeCompletenessCutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param applyfilters Optional string specifying filtration setting "combos", used as a shorthand for setting the featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff arguments in JAMS plotting functions. If NULL, none of these arguments are set if not specified. Permissible values for applyfilters are "light", "moderate" or "stringent". The actual values vary whether the SummarizedExperiment object is taxonomical (LKT) or not. For a taxonomical SummarizedExperiment object, using "light" will set featcutoff=c(50, 5), GenomeCompletenessCutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(250, 15), GenomeCompletenessCutoff=c(10, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(2000, 15), GenomeCompletenessCutoff=c(30, 10), minl2fc=2, minabscorrcoeff=0.8. For non-taxonomical (i.e. functional) SummarizedExperiment objects, using "light" will set featcutoff=c(0, 0), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(50, 15), minl2fc=2.5, minabscorrcoeff=0.8. When using applyfilters, one can still set one or more of featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff, which will then take the user set value in lieu of those set by the applyfilters shorthand. Default is light.

#' @param discard_SDoverMean_below Numeric value setting the minimum standard deviation over mean (sd/mean) value cutoff for a feature to be kept. Features with an sd/mean value smaller than this threshold will be discarded. When NULL, this filtration is not applied. Default is NULL.

#' @param showpval Requires a logical value. If set to TRUE, a text label with the p-values (raw and adjusted) for each feature will be shown on the right of each row of the heatmap. Default is TRUE.

#' @param showroundedpval Requires a logical value. If set to TRUE, the text label with the p-values (raw and adjusted) for each feature will be rounded to three decimal places. Default is TRUE.

#' @param showonlypbelow Numeric value setting the maximum p-value for a feature to be shown on a heatmap. This only applies to "comparative" heatmaps (see hmtype). Features with p-values above this threshold will not be shown. When NULL, all features surviving other filtration criteria will be shown. Default is NULL.

#' @param adj_pval_for_threshold Requires a logical value. If set to TRUE, then if setting a numerical value for p-value cutoff with showonlypbelow (q.v.), only **adjusted** p-values below or equal to the threshold will be shown. Default is FALSE, i.e. use a raw p-value for cutoff.

#' @param minl2fc Numeric value setting the minimum absolute log2 fold change value to report within a heatmap. This only applies to "comparative" heatmaps (see hmtype) when the variable being compared (see compareby) has exactly two classes. When NULL, this filtration is not applied and features with minimum log2 foldchanges of any value down to 0 will be shown. Default is NULL.

#' @param maxl2fc Numeric value setting the maximum absolute log2 fold change value to report within a heatmap. This only applies to "comparative" heatmaps (see hmtype) when the variable being compared (see compareby) has exactly two classes. When NULL, this filtration is not applied and features with maximum log2 foldchanges of any value up to Inf will be shown. Default is NULL.

#' @param fun_for_l2fc String specifying the mathematical function for aggregating samples within each group when calculating the log2 fold change for each feature. This only applies to "comparative" heatmaps (see hmtype) when the variable being compared (see compareby) has exactly two classes. Permissible entries are "sum", "mean", "median" or "geom_mean". Default is geom_mean. For "geom_mean", the geometric mean within each group is calculated with exp(mean(log((x1, x2, x3, ...) + 1))).

#' @param showl2fc Requires a logical value. If set to TRUE, a text label with the log2 fold change between groups of each feature when applicable will be shown on the right of each row of the heatmap. Default is TRUE. See compareby and fun_for_l2fc.

#' @param showGram Requires a logical value. If set to TRUE, if the SummarizedExperiment object is taxonomical (see ExpObj), annotations with the Phylum and predicted Gram cell wall category of each feature will be plot to the left of each row of the heatmap. Default is FALSE.

#' @param show_GenomeCompleteness Requires a logical value. When TRUE (the default), if the SummarizedExperiment object is taxonomical (see ExpObj), annotations with the Phylum and predicted Gram cell-wall category of each feature will be plot to the left of each row of the heatmap. Default is FALSE.

#' @param addtit Optional string with text to append to heatmap main title. Default is NULL.

#' @param assay_for_matrix String specifying the SummarizedExperiment assay to be used for the heatmap. Permissible values are "BaseCounts" or "GeneCounts". "BaseCounts" (the default) will use the basepair counts for each feature (either taxonomical or functional). These values will be normalized into relative abundance in PPM unless specified by the normalization argument (see normalization and PPM_normalize_to_bases_sequenced). When using "GeneCounts" (only available in non-taxonomical SummarizedExperiment objects) the *number of genes* annotated as each feature will be used. The heatmap will be plot with a scale of 0 to the maximum number of genes for a single feature on the heatmap. For instance, using "GeneCounts" for, let's say, an ECNumber SummarizedExperiment will plot the number of genes bearing each Enzyme Commission Number annotation within each sample. Default is "BaseCounts".

#' @param normalization String specifying if the BaseCounts for the assay should be normalized or not. Permissible values are "relabund" and "compositions". When using "relabund", the relative abundance of each feature will be calculated in Parts per Million (PPM) by dividing the number of bases covering each feature by the sum of each sample column **previous to any filtration**. See also PPM_normalize_to_bases_sequenced for details.  When using "compositions", the counts matrix will be transformed using the clr function of the compositions package. Please install this package independently of JAMS as it is not a JAMS dependency.

#' @param PPM_normalize_to_bases_sequenced Requires a logical value. Non-filtered JAMS feature counts tables (the BaseCounts assay within SummarizedExperiment objects) always includes unclassified taxonomical features (for taxonomical SummarizedExperiment objects) or unknown/unattributed functional features (for non-taxonomical SummarizedExperiment objects), so the relative abundance for each feature (see normalization) will be calculated in Parts per Million (PPM) by dividing the number of bases covering each feature by the sum of each sample column **previous to any filtration**. Relative abundances are thus representative of the entirety of the genomic content for taxonomical objects, whereas for non-taxonomical objects, strictly speaking, it is the abundance of each feature relative to only the coding regions present in the metagenome, even if these are annotationally unatributed. In other words, intergenic regions are not taken into account. In order to relative-abundance-normalize a **non-taxonomical** SummarizedExperiment object with the total genomic sequencing content, including non-coding regions, set PPM_normalize_to_bases_sequenced = TRUE. Default is FALSE.

#' @param scaled Requires a logical value. If set to TRUE the z-scores for each row (each feature) will be plot on the heatmap rather than their relative abundances.

#' @param numthreads Numeric value setting the number of threads to use for any multi-threaded process within this function. The default is 1.

#' @param statsonlog Requires a logical value. If set to TRUE, sample statistics will be calculated on a log2 transformed relative abundance table. Default is FALSE, meaning statistics (p-values or variance) will be calculated on relative abundance values in PPM. See also hmtype.

#' @param ignoreunclassified Requires a logical value. If set to TRUE, for taxonomical SummarizedExperiment objects, the feature "LKT__Unclassified" will be omitted from being shown. In the case of non-taxonomical SummarizedExperiment objects, the completely unannotated features will be omitted. For example, for an ECNumber SummarizedExperiment object, genes **without** an Enzyme Commission Number annotation (feature "EC_none") will not be shown. Statistics are, however, computed taking the completely unclassifed feature into account, so p-values will not change.

#' @param returnstats Requires a logical value. If set to TRUE, this function will return a named list of the statistical computations obtained for **all** features within each subset (see subsetby) with the relevant statistic (MannWhitneyWilcoxon p-value, ANOVA p-value, variance, SD, rank, etc). Default is FALSE. Note that when FALSE, this function does not return anything, rather it plots to the device using the ComplexHeatmaps package.

#' @param class_to_ignore String or vector specifying any classes which should lead to samples being excluded from the comparison within the variable passed to compareby. Default is N_A. This means that within any metadata variable passed to compareby containing the "N_A" string within that specific variable, the sample will be dropped from that comparison.

#' @param no_underscores Requires a logical value. If set to TRUE, removes underscores from taxonomical feature names on a heatmap. For example, "LKT__s__Staphylococcus_aureus" would be plot as "LKT s Staphylococcus aureus". Default is FALSE.

#' @export

plot_relabund_heatmap <- function(ExpObj = NULL, glomby = NULL, hmtype = "exploratory", samplesToKeep = NULL, featuresToKeep = NULL, subsetby = NULL, compareby = NULL, wilcox_paired_by = NULL, invertbinaryorder = FALSE, hmasPA = FALSE, threshPA = 0, ntop = NULL, splitcolsby = NULL, cluster_column_slices = TRUE, column_split_group_order = NULL, ordercolsby = NULL, cluster_samples_per_heatmap = TRUE, cluster_features_per_heatmap = TRUE, colcategories = NULL, textby = NULL, label_samples = TRUE, cluster_rows = TRUE, row_order = NULL, max_rows_in_heatmap = 50, applyfilters = "light", featcutoff = NULL, GenomeCompletenessCutoff = NULL, discard_SDoverMean_below = NULL, maxl2fc = NULL, minl2fc = NULL, fun_for_l2fc = "geom_mean", adj_pval_for_threshold = FALSE, showonlypbelow = NULL, showpval = TRUE, showroundedpval = TRUE, showl2fc = TRUE, showGram = FALSE, show_GenomeCompleteness = TRUE, use_checkM2_style_GenomeCompleteness = TRUE, addtit = NULL, assay_for_matrix = "BaseCounts", normalization = "relabund", asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, scaled = FALSE, cdict = NULL, maxnumheatmaps = NULL, numthreads = 1, statsonlog = FALSE, ignoreunclassified = TRUE, returnstats = FALSE, class_to_ignore = "N_A", no_underscores = FALSE, ...){

    #Account for JAMS2 spaces
    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")

    #Hardwire the minimum number of samples for a side by side secondary heatmap. This may be in a later version added as a function argument.
    threshold_for_double_plot <- 7

    secondaryheatmap <- "none"

    #Test for silly stuff
    if ((hmtype %in% c("comparative", "PA")) && (is.null(compareby))){
        stop("If rows are to be selected by highest significant difference between classes in a discrete category, you must determine the category using the argument *compareby*")
    }

    #Fix metadata classes and remove classes to ignore, if present
    if (hmtype %in% c("comparative", "PA")){
        variables_to_fix <- c(compareby, subsetby, ordercolsby)
    } else {
        variables_to_fix <- c(subsetby, ordercolsby)
    }

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, glomby = glomby, variables_to_fix = variables_to_fix, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    #Set applyfilters to null if featuresToKeep is set
    if (!is.null(featuresToKeep)){
        applyfilters <- NULL
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, maxl2fc = maxl2fc, minl2fc = minl2fc)

    if (assay_for_matrix == "GeneCounts"){
        flog.warn("Counts matrix used for heatmap will represent the *number of genes* for each feature, rather than its relative abundance. For using relative abundance (default), set assay_for_matrix = \"BaseCounts\"")
        asPPM <- FALSE
        global_FTK <- rownames(obj)
        hmtypemsg <- "Gene count Heatmap"
    } else {
        hmtypemsg <- "Relative Abundance Heatmap"

        #Obtain "global" list of features if applyfilters is not null, or any of the filtering arguments is not c(0,0) as to maintain same features within subsets, if any.
        if (any(c(all(presetlist$featcutoff != c(0,0)), all(presetlist$GenomeCompletenessCutoff != c(0,0)), !is.null(applyfilters)))){
            currobj <- filter_experiment(SEobj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = NULL, featuresToKeep = NULL, normalization = "relabund", PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, discard_SDoverMean_below = discard_SDoverMean_below)
            global_FTK <- rownames(currobj)
        } else {
            global_FTK <- rownames(obj)
        }
    }

    if (!(is.null(subsetby))){
        subset_list <- multiple_subsetting_sample_selector(SEobj = obj, phenotable = NULL, subsetby = subsetby, compareby = compareby, cats_to_ignore = class_to_ignore)
        #Ensure there are only valid subsets which will plot on a heatmap.
        subset_df <- subset_list$Subsets_stats
        subset_df <- subset_df[which(subset_df$Subset_Tier_Level != 0), , drop = FALSE]
        if (any(subset_df$Num_samples_in_subset < 2)){
            #Some subsets have less than 2 samples. Eliminate and report.
            LowSampSubsets <- subset_df[which(subset_df$Num_samples_in_subset < 2), "Subset_Tier_Class_Name"]
            flog.warn(paste("Subsets", paste0(LowSampSubsets, collapse = ", "), "contain less than 2 samples, and will thus be omitted."))
            subset_df <- subset_df[which(subset_df$Num_samples_in_subset >= 2), , drop = FALSE]
        }
        #Test for valid subsets
        if (nrow(subset_df) < 1){
            flog.warn("There are no surviving subset points. Defaulting to no subsetting.")
            subset_points <- "none"
            subsetby <- NULL
        } else {
            subset_points <- subset_df$Subset_Tier_Class_Name
        }
    } else {
        subset_points <- "none"
    }

    if (is.null(colcategories)){
        if (!is.null(cdict)){
            colcategories <- names(cdict)
        } else {
            #Include anything with between 2 and 10 classes in it
            colcategories <- colnames(colData(obj))[which((sapply(1:ncol(colData(obj)), function (x) { length(levels(as.factor(colData(obj)[, x]))) })) < 10 & (sapply(1:ncol(colData(obj)), function (x) { length(levels(as.factor(colData(obj)[, x]))) })) > 1)]
        }
    }

    #ensure that what the user is asking for is within the metadata
    colcategories <- colcategories[colcategories %in% colnames(colData(obj))]

    #Initialize Stats Vector list
    svec <- list()
    s <- 1
    n <- 1

    #Stop if there is no data to plot
    if ((dim(obj)[1] * dim(obj)[2]) < 4){
        flog.warn("There are fewer than 4 cells in the SummarizedExperiment object matrix, impossible to draw a heatmap. Check your input and try again.")

        return(NULL)
    }

    #subset by metadata column
    for (sp in 1:length(subset_points)){
        if (!(is.null(subsetby))){
            samplesToKeep_sp <- subset_list[[subset_points[sp]]]
            flog.info(paste("Plotting within", subset_points[sp]))
            subsetname <- subset_points[sp]
        } else {
            samplesToKeep_sp <- rownames(colData(obj))
            subsetname <- "no_sub"
        }

        #See if there are enough samples and features to go ahead
        proceed <- TRUE

        if (length(samplesToKeep_sp) < 2){
            skip_HM_reason <- paste0("There are fewer than 2 samples within", subsetname, ", impossible to plot a heatmap.")
            proceed <- FALSE
        }
        curr_pt <- colData(obj)[samplesToKeep_sp, ]

        if (any(c(!is.null(compareby), (hmtype != "exploratory")))) {
            #investigate if there are at least two of each class of things to compare to go ahead
            if (length(unique(curr_pt[ , compareby])) < 2){
                skip_HM_reason <- paste0("There are fewer than 2 classes within variable ", compareby, ", impossible to obtain a p-value. Try setting hmtype = 'exploratory'.")
                proceed <- FALSE
            }

            if (length(unique(curr_pt[ , compareby])) == 2){
                #Comparison is binary. Is there at least 2 of each class to get a p-value?
                if ((min(table(curr_pt[ , compareby]))) < 2){
                    flog.warn(paste0("There are less than 2 samples within at least one class of variable ", compareby, ", impossible to obtain a p-value. Setting hmtype = 'exploratory'."))
                    hmtype <- "exploratory"
                }
            }
        }

        if (proceed){

            if (hmtype == "exploratory"){
                stattype <- "variance"
                if (is.null(ntop)){
                    ntop <- max_rows_in_heatmap
                }
            } else if (hmtype == "PA"){
                stattype <- "PA"
            } else {
                stattype <- "auto"
            }

            if (assay_for_matrix == "BaseCounts"){

                currobj <- filter_experiment(SEobj = obj, featcutoff = c(0,0), samplesToKeep = samplesToKeep_sp, featuresToKeep = global_FTK, normalization = normalization, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = c(0,0), discard_SDoverMean_below = NULL, flush_out_empty_samples = FALSE)

            } else if (assay_for_matrix == "GeneCounts"){
                currobj <- obj[global_FTK, samplesToKeep_sp]
            }

            flog.info(paste("Plotting", analysisname, "heatmap."))
            fullheatmap_column_order <- NULL
            fullheatmap_column_dend <- NULL
            fullheatmap_row_order <- NULL
            fullheatmap_row_dend <- NULL

            #Compose an appropriate title for the plot
            if (length(unique(subset_points)) > 1){
                maintit <- paste(hmtypemsg, analysisname, paste("within", subset_points[sp]), sep = " | ")
            } else {
                maintit <- paste(hmtypemsg, analysisname, sep = " | ")
            }
            if (!is.null(addtit)) {
                maintit <- paste(addtit, maintit, sep = "\n")
            }

            #Get counts matrix
            if (normalization == "relabund"){
                countmat <- as.matrix(assays(currobj)[["PPM"]])
            }

            if (ignoreunclassified == TRUE){
               dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified", "hypothetical protein", "putative protein")
               rowsToKeep <- which(!(rownames(countmat) %in% dunno))
               countmat <- countmat[rowsToKeep, ]
            }

            #Rename rows to include description if not taxonomic data or MetaCyc which has enormous descriptions
            #if (!(analysis %in% c("LKT", "MetaCyc"))){
            if (!analysis %in% taxonomic_spaces){
                feattable <- rowData(currobj)
                feattable$Feature <- paste(feattable$Accession, feattable$Description, sep = "-")
                rownames(countmat) <- feattable$Feature[match(rownames(countmat), feattable$Accession)]
                #Paste descriptions to row_order vector if not LKT space.
                if (!is.null(row_order)){
                    row_order <- row_order[row_order %in% rownames(feattable)]
                    row_order <- paste(feattable[row_order, "Accession"], feattable[row_order, "Description"], sep = "-")
                }
            }

            matrixSamples <- colnames(countmat)
            matrixRows <- rownames(countmat)

            topcats <- nrow(countmat)
            if (!(is.null(ntop))) {
                topcats <- min(topcats, ntop)
            }

            #Calculate matrix stats and get new matrix.
            if (stattype == "variance"){
                matstats <- calculate_matrix_stats(countmatrix = countmat, stattype = "variance")
                #Bank raw stats
                matstatsbank <- as.data.frame(matstats)
                if ("GenomeCompleteness" %in% names(assays(currobj))){
                    matstatsbank$Taxa <- rownames(matstatsbank)
                    if (!use_checkM2_style_GenomeCompleteness){
                        genomecompletenessmat <- as.matrix(assays(currobj)$GenomeCompleteness) + as.matrix(assays(currobj)$GenomeContamination)
                    } else {
                        genomecompletenessmat <- as.matrix(assays(currobj)$GenomeCompleteness)
                    }
                    genomecompletenessdf <- as.data.frame(genomecompletenessmat)
                    genomecompletenessdf$MedianGenomeComp <- rowMedians(genomecompletenessmat)
                    genomecompletenessdf$SDGenomeComp <- rowSds(genomecompletenessmat)
                    genomecompletenessdf$Taxa <- rownames(genomecompletenessdf)
                    genomecompletenessdfmedians <- genomecompletenessdf[, c("Taxa", "MedianGenomeComp", "SDGenomeComp")]
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdfmedians, by = "Taxa")
                    rownames(matstatsbank) <- matstatsbank$Taxa
                    matstatsbank$Taxa <- NULL
                }
                svec[[s]] <- matstatsbank
                ffeatmsg <- paste0("Number of features assessed = ", nrow(matstats))
                matstats$Colour <- rep("black", nrow(matstats))
                countmat2 <- countmat[rownames(matstats), ]

                #Transform to log2 space
                if (assay_for_matrix == "BaseCounts"){
                    countmat2 <- convert_matrix_log2(mat = countmat2, transformation = "to_log2")
                }

                #Turn off row clustering if explicitly setting row_order
                if (!is.null(row_order)){
                    cluster_rows <- FALSE
                }

                if (all(c(cluster_rows, (any(!(c(cluster_samples_per_heatmap, cluster_features_per_heatmap))))))){
                    flog.info("Clustering samples and features using entire matrix to obtain sample and feature order for all heatmaps.")
                    #create a heatmap from the entire count matrix for getting column order.
                    htfull <- suppressMessages(Heatmap(countmat2, name = "FullHM", column_title = "FullHM", column_names_gp = gpar(fontsize = 1), cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 1, col = "black"), use_raster = TRUE))

                    fullheatmap_column_order <- suppressWarnings(column_order(htfull))
                    fullheatmap_column_dend <- suppressWarnings(column_dend(htfull))
                    fullheatmap_row_order <- suppressWarnings(row_order(htfull))
                    fullheatmap_row_dend <- suppressWarnings(row_dend(htfull))
                }

                if (all(c(!cluster_features_per_heatmap, !is.null(fullheatmap_row_order)))){
                    countmat2 <- countmat2[fullheatmap_row_order, ]
                    matstats <- matstats[rownames(countmat2), ]
                }

                if (!is.null(row_order)){
                    #Ensure all elements in the vector are actually in the matrix
                    row_order_have <- row_order[row_order %in% rownames(countmat2)]
                    #warn if any are missing
                    if (length(row_order) != length(row_order_have)) {
                        flog.warn("Some elements passed to row_order were not found. Check filtration criteria and/or input SEobj")
                    }
                    if (length(row_order_have) < 2){
                        flog.warn("There are fewer than 2 features available. Impossible to draw a heatmap. Aborting now.")
                        return(NULL)
                    }
                    countmat2 <- countmat2[row_order_have, ]
                    matstats <- matstats[rownames(countmat2), ]
                }

                #Create a list of matrices each of maximum 50 rows
                rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / max_rows_in_heatmap))
                matlist <- lapply(1:length(rowlist), function(x){countmat2[rowlist[[x]], ]})
                statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                rowlblcol_list <- lapply(1:length(rowlist), function(x){ rep("black", length(rowlist[[x]])) })
                stattit <- paste("Top", topcats, "most variant features across samples")
                statmsg <- stattype

            } else {
                if (!is.null(wilcox_paired_by)){
                    flog.info(paste("Will attempt to pair samples by", wilcox_paired_by, "for Mann-Whitney-Wilcoxon test"))
                }
                classesdf <- make_classes_df(curr_pt = colData(currobj), compareby = compareby, wilcox_paired_by = wilcox_paired_by)

                discretenames <- sort(unique(classesdf$cl))
                if ("wilcox_pairs" %in% colnames(classesdf)){
                    flog.info(paste("Mann-Whitney-Wilcoxon test between", discretenames[1], "and", discretenames[2], "will be paired by", wilcox_paired_by))
                }

                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = FALSE, statsonlog = statsonlog, stattype = stattype, classesdf = classesdf, invertbinaryorder = invertbinaryorder, numthreads = numthreads, threshPA = threshPA, fun_for_l2fc = fun_for_l2fc)

                ffeatmsg <- paste0("Number of features assessed = ", nrow(matstats))

                #Bank raw stats
                matstatsbank <- as.data.frame(matstats)
                if ("GenomeCompleteness" %in% names(assays(currobj))){
                    matstatsbank$Taxa <- rownames(matstatsbank)
                    if (!use_checkM2_style_GenomeCompleteness){
                        genomecompletenessmat <- as.matrix(assays(currobj)$GenomeCompleteness) + as.matrix(assays(currobj)$GenomeContamination)
                    } else {
                        genomecompletenessmat <- as.matrix(assays(currobj)$GenomeCompleteness)
                    }
                    genomecompletenessdf <- as.data.frame(genomecompletenessmat)
                    genomecompletenessdf$MedianGenomeComp <- rowMedians(genomecompletenessmat)
                    genomecompletenessdf$SDGenomeComp <- rowSds(genomecompletenessmat)
                    genomecompletenessdf$Taxa <- rownames(genomecompletenessdf)
                    genomecompletenessdfmedians <- genomecompletenessdf[, c("Taxa", "MedianGenomeComp", "SDGenomeComp")]
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdfmedians, by = "Taxa")
                    rownames(matstatsbank) <- matstatsbank$Taxa
                    matstatsbank$Taxa <- NULL
                }
                svec[[s]] <- matstatsbank

                #Account for the fact that if there is only a single class within compareby variable, stats may have been coerced to variance.
                if (matstats$Method[1] == "variance") {
                    matstats$Colour <- rep("black", nrow(matstats))
                    countmat2 <- countmat[rownames(matstats), ]
                    if (assay_for_matrix != "GeneCount"){
                        #Transform to log2 space
                        countmat2 <- convert_matrix_log2(mat = countmat2, transformation = "to_log2")
                    }

                    if (all(c(cluster_rows, (any(!(c(cluster_samples_per_heatmap, cluster_features_per_heatmap))))))){
                        flog.info("Clustering samples and features using entire matrix to obtain sample and feature order for all heatmaps.")
                        #create a heatmap from the entire count matrix for getting column order.
                        htfull <- Heatmap(countmat2, name = "FullHM", column_title = "FullHM", column_names_gp = gpar(fontsize = 1), cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 1, col = "black"), use_raster = TRUE)

                        fullheatmap_column_order <- column_order(htfull)
                        fullheatmap_column_dend <- column_dend(htfull)
                        fullheatmap_row_order <- row_order(htfull)
                        fullheatmap_row_dend <- row_dend(htfull)
                    }

                    if (all(c(cluster_rows, (!cluster_features_per_heatmap)))){
                        countmat2 <- countmat2[fullheatmap_row_order, ]
                        matstats <- matstats[rownames(countmat2), ]
                    }

                    if (!is.null(row_order)){
                        #Ensure all elements in the vector are actually in the matrix
                        row_order_have <- row_order[row_order %in% rownames(countmat2)]
                        #warn if any are missing
                        if (length(row_order) != length(row_order_have)) {
                            flog.warn("Some elements passed to row_order were not found. Check filtration criteria and/or input SEobj")
                        }
                        if (length(row_order_have) < 2){
                            flog.warn("There are fewer than 2 features available. Impossible to draw a heatmap. Aborting now.")
                            return(NULL)
                        }
                        countmat2 <- countmat2[row_order_have, ]
                        matstats <- matstats[rownames(countmat2), ]
                    }

                    #Create a list of matrices each of maximum 50 rows
                    topcats <- min(nrow(countmat2), max_rows_in_heatmap)
                    rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / max_rows_in_heatmap))
                    matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                    statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                    rowlblcol_list <- lapply(1:length(rowlist), function(x){rep("black", length(rowlist[[x]]))})
                    stattit <- paste("Top", topcats, "most variant features across samples")
                    statmsg <- paste("Var", compareby, sep = "_")

                } else {
                    #Stats are not variance
                    if (all(c((!(is.null(presetlist$minl2fc))), ("l2fc" %in% colnames(matstats)), (hmtype != "PA")))) {
                        #Check if there is enough leftover after filter.
                        if (length(which(matstats$absl2fc > presetlist$minl2fc)) < 2){
                            flog.info(paste("There are fewer than 2 features which have > ", presetlist$minl2fc, "l2fc."))
                            #Redefine min l2fc to whatever is the lowest available keeping between 2 and 40 features to plot.
                            presetlist$minl2fc <- matstats[order(matstats$absl2fc, decreasing = TRUE),][min(max(length(matstats$absl2fc), 2), 40), ]$absl2fc
                            flog.info(paste("Resetting minimum l2fc to", round(presetlist$minl2fc, 2)))
                        }
                        minl2fcmsg <- paste("log2foldchange >", round(presetlist$minl2fc, 2))
                        matstats <- subset(matstats, absl2fc > presetlist$minl2fc)
                    } else {
                        minl2fcmsg <- "log2foldchange > 0"
                    }

                    if (all(c((!(is.null(presetlist$maxl2fc))), ("l2fc" %in% colnames(matstats)), (hmtype != "PA")))) {
                        #Check if there is enough leftover after filter.
                        if (length(which(matstats$absl2fc < presetlist$maxl2fc)) < 2){
                            flog.info(paste("There are fewer than 2 features which have < ", presetlist$maxl2fc, "l2fc."))
                            #Redefine max l2fc to whatever is the lowest available keeping between 2 and 40 features to plot.
                            presetlist$maxl2fc <- matstats[order(matstats$absl2fc, decreasing = FALSE),][min(max(length(matstats$absl2fc), 2), 40), ]$absl2fc
                            flog.info(paste("Resetting maximum l2fc to", round(presetlist$maxl2fc, 2)))
                        }
                        maxl2fcmsg <- paste("log2foldchange <", round(presetlist$maxl2fc, 2))
                        matstats <- subset(matstats, absl2fc < presetlist$maxl2fc)
                    } else {
                        maxl2fcmsg <- "log2foldchange > Inf"
                    }

                    if ("l2fc" %in% colnames(matstats)){
                        matstats$Colour <- ifelse(matstats$l2fc < 0, "#900000", "#000000")
                        statmsg <- paste("MWW", compareby, sep="_")
                    } else if ("OddsRatio" %in% colnames(matstats)){
                        matstats$Colour <- ifelse(matstats$OddsRatio < 1, "#900000", "#000000")
                        statmsg <- paste("Fisher", compareby, sep="_")
                    } else {
                        matstats$Colour <- rep("#000000", nrow(matstats))

                        if (matstats$Method[1] == "permanova"){
                            statmsg <- paste("PERMANOVA", compareby, sep="_")
                        } else if (matstats$Method[1] == "anova"){
                            statmsg <- paste("ANOVA", compareby, sep="_")
                        } else if (matstats$Method[1] == "variance"){
                            statmsg <- paste("Var", compareby, sep="_")
                        }
                    }

                    binary_directionality <- (discretenames[1:2][c(!invertbinaryorder, invertbinaryorder)])
                    if (matstats$Method[1] == "fisher"){
                        binary_directionality_msg <- paste("Odds Ratio > 1 means enriched in", binary_directionality)
                    } else {
                        binary_directionality_msg <- paste("Positive l2fc means increased in", binary_directionality)
                    }

                    #Obtain a matrix that represents cells in the heatmap
                    countmat2 <- as.matrix(countmat[rownames(matstats), ])
                    if (all(c((hmtype == "PA"), hmasPA))){
                        #Transform to Presence/Absence according to threshold space
                        countmat2 <- convert_matrix_PA(mat = countmat2, threshPA = threshPA)
                    } else {
                        if (assay_for_matrix != "GeneCounts"){
                            #Transform to log2 space
                            countmat2 <- convert_matrix_log2(mat = countmat2, transformation = "to_log2")
                        }
                    }

                    #Decide row and sample ordering
                    if (all(c(cluster_rows, (any(!(c(cluster_samples_per_heatmap, cluster_features_per_heatmap))))))){
                        flog.info("Clustering samples and features using entire matrix to obtain sample and feature order for all heatmaps.")
                        #create a heatmap from the entire count matrix for getting column order.
                        htfull <- Heatmap(countmat2, name = "FullHM", column_title = "FullHM", column_names_gp = gpar(fontsize = 1), cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 1, col = "black"), use_raster = TRUE)

                        fullheatmap_column_order <- column_order(htfull)
                        fullheatmap_column_dend <- column_dend(htfull)
                        fullheatmap_row_order <- row_order(htfull)
                        fullheatmap_row_dend <- row_dend(htfull)
                    }

                    if (all(c(cluster_rows, (!cluster_features_per_heatmap)))) {
                        countmat2 <- countmat2[fullheatmap_row_order, ]
                        matstats <- matstats[rownames(countmat2), ]
                    }

                    if (!is.null(row_order)){
                        #Ensure all elements in the vector are actually in the matrix
                        row_order_have <- row_order[row_order %in% rownames(countmat2)]
                        #warn if any are missing
                        if (length(row_order) != length(row_order_have)) {
                            flog.warn("Some elements passed to row_order were not found. Check filtration criteria and/or input SEobj")
                        }
                        if (length(row_order_have) < 2){
                            flog.warn("There are fewer than 2 features available. Impossible to draw a heatmap. Aborting now.")
                            return(NULL)
                        }
                        countmat2 <- countmat2[row_order_have, ]
                        matstats <- matstats[rownames(countmat2), ]
                    }

                    #If adj_pval_for_threshold is set to auto, then find out which is best and re-set it to either TRUE or FALSE

                    if (("pval" %in% colnames(matstats)) && adj_pval_for_threshold == "auto"){
                        propsigadj <- length(which(matstats$padj_fdr < 0.05)) / length(matstats$padj_fdr)
                        propsignonadj <- length(which(matstats$pval < 0.05)) / length(matstats$pval)
                        fracsigadj <- propsigadj / propsignonadj
                        if ((!is.na(fracsigadj)) && (fracsigadj > 0.2)){
                            adj_pval_for_threshold <- TRUE
                        } else {
                            adj_pval_for_threshold <- FALSE
                        }
                        flog.info(paste("P-value adjustment set to auto. Upon adjustment using FDR, the proportion still below 0.05 after adjustment is", round((fracsigadj * 100), 1), "% of unadjusted p-values, so p-value adjustment is set to", as.character(adj_pval_for_threshold)))
                    }

                    if (("pval" %in% colnames(matstats)) && adj_pval_for_threshold != TRUE){
                        sigmeas <- "pval"
                    } else if (("pval" %in% colnames(matstats)) && adj_pval_for_threshold == TRUE){
                        sigmeas <- "padj_fdr"
                    }

                    if (is.null(showonlypbelow)){
                        if (nrow(countmat2) < topcats){
                            topcats <- nrow(countmat2)
                            flog.info("There are not enough features matching the criteria imposed.")
                            flog.info(paste("Showing top", topcats, "features."))
                        }
                        countmat2 <- countmat2[1:topcats, ]
                        #Create a list of matrices each of maximum ~50 rows
                        if((topcats %% max_rows_in_heatmap) == 1){
                            chunksize <- (max_rows_in_heatmap - 2)
                        } else {
                            chunksize <- max_rows_in_heatmap
                        }

                        rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / chunksize))
                        matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                        statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                        rowlblcol_list <- lapply(1:length(rowlist), function(x) { matstats$Colour[rowlist[[x]]] })

                        if (matstats$Method[1] == "MannWhitneyWilcoxon") {
                            stattit <- paste("Top", topcats, "different between", compareby, "using Mann-Whitney-Wilcoxon U-test")
                            if ("wilcox_pairs" %in% colnames(classesdf)){
                                stattit <- paste(stattit, paste0("Mann-Whitney-Wilcoxon U-test paired by ", wilcox_paired_by), sep = "\n")
                            }
                        } else if (matstats$Method[1] == "permanova"){
                            stattit <- paste("Top", topcats, "different between", compareby, "using PERMANOVA")
                        } else if (matstats$Method[1] == "anova"){
                            stattit <- paste("Top", topcats, "different between", compareby, "using ANOVA")
                        } else if (matstats$Method[1] == "fisher"){
                            stattit <- paste("Top", topcats, "Present/Absent between", compareby, "using Fishers test")
                        }

                    } else {
                        if (adj_pval_for_threshold != TRUE){
                            rowcutoff <- which(matstats$pval <= showonlypbelow)
                        } else {
                            rowcutoff <- which(matstats$padj_fdr <= showonlypbelow)
                        }

                        #Limit number of features to requested number or number available
                        rowcutoff <- rowcutoff[1:(min(topcats, length(rowcutoff)))]

                        #Must have at least two rows in a matrix to plot a heatmap
                        if (length(rowcutoff) > 1){
                            #Account for the fact that if the remainder of a chunk of 50 is 1 then a heatmap with a single feature cannot be drawn. Oh my, I have seen everything havent I...
                            if((length(rowcutoff) %% max_rows_in_heatmap) == 1){
                                chunksize <- (max_rows_in_heatmap - 2)
                            } else {
                                chunksize <- max_rows_in_heatmap
                            }

                            rowlist <- split(rowcutoff, ceiling(seq_along(rowcutoff) / chunksize))
                            matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                            statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                            if (matstats$Method[1] == "MannWhitneyWilcoxon") {
                                stattit <- paste(sigmeas, "<", showonlypbelow, "different between", compareby, "using Mann-Whitney-Wilcoxon")
                                if ("wilcox_pairs" %in% colnames(classesdf)){
                                    stattit <- paste(stattit, paste0("Mann-Whitney-Wilcoxon U-test paired by ", wilcox_paired_by), sep = "\n")
                                }
                            } else if (matstats$Method[1] == "permanova"){
                                stattit <- paste(sigmeas, "<", showonlypbelow, "different between", compareby, "using PERMANOVA")
                            } else if (matstats$Method[1] == "anova"){
                                stattit <- paste(sigmeas, "<", showonlypbelow, "different between", compareby, "using ANOVA")
                            } else if (matstats$Method[1] == "fisher"){
                                stattit <- paste(sigmeas, "<", showonlypbelow, "Present/Absent between", compareby, "using Fishers test")
                            }
                        } else {
                            #allIhave<-min(nrow(countmat2), 30)
                            #stattit<-paste("No values", "p <", showonlypbelow, "between", compareby,"\nShowing top", allIhave, "with any p")
                            #matlist<-list(countmat2[1:allIhave, ])
                            #rowlblcol_list<-list(matstats$Colour[1:allIhave])
                            #Report nothing was found
                            matlist <- NULL
                            if (!(is.na(matstats$Method[1]))) {
                                if (matstats$Method[1] == "MannWhitneyWilcoxon") {
                                    stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using MannWhitneyWilcoxon"))
                                    if ("wilcox_pairs" %in% colnames(classesdf)){
                                        stattit <- paste(stattit, paste0("Mann-Whitney-Wilcoxon U-test paired by ", wilcox_paired_by), sep = "\n")
                                    }
                                } else if (matstats$Method[1] == "permanova"){
                                    stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using PERMANOVA"))
                                } else if (matstats$Method[1] == "anova"){
                                    stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using ANOVA"))
                                } else if (matstats$Method[1] == "fisher"){
                                    stattit <- c(paste(sigmeas, "<", showonlypbelow, "Present/Absent between"), paste(compareby, "using Fishers test"))
                                    if (threshPA != 0){
                                        stattit <- c(stattit, paste("Presence is considered PPM >=", threshPA))
                                    }
                                } #End conditional for getting stats title
                            } else {
                                #Stats failed for some reason, like there are no stat features fulfilling the filtering criteria
                                stattit <- NULL
                            }
                        } #End conditional that there are two rows or more to plot
                    } #End conditional for only showing feats with certain pval
                } #End conditional that stats method is not variance
            } #End conditional for choosing stats method

            #Name stats in svec
            stattitle <- paste(analysisname, statmsg, subsetname, sep = "_")
            names(svec)[s] <- stattitle
            s <- s + 1

            #If the matrix list is not empty, plot the heatmap. If there is nothing to plot, say so.
            if (length(matlist) > 0){
                #Cycle through list of matrices to transform into heatmaps.
                if (!(is.null(maxnumheatmaps))){
                    #Prune graphics list to largest number allowed
                    mhm <- min(maxnumheatmaps, length(matlist))
                } else {
                    mhm <- length(matlist)
                }

                for (hm in 1:mhm){
                    #Regenerate the current matrix being plot from matrix list
                    mathm <- matlist[[hm]]

                    stathm <- statslist[[hm]]
                    rowlblcol <- stathm$Colour
                    ht1fs <- 10

                    #Plot the heatmap
                    fontsizexy <- hm_fontsize_computer(mat_rownames = rownames(mathm), mat_colnames = colnames(mathm), upper_n = 400, upper_fs = 0.1, lower_n = 10, lower_fs = 9, cex = 0.42)

                    fontsizex <- fontsizexy[1]
                    fontsizey <- fontsizexy[2]

                    #Fix row names, add a carriage return (\n) if over 40 characters and not LKT
                    if (!analysis %in% taxonomic_spaces){
                        rownames(mathm) <- sapply(rownames(mathm), function(x) { split_featname(featname = x, thresh_featname_split = 60) } )
                    }

                    if (label_samples == FALSE){
                        fontsizex <- 0
                    }

                    hmdf <- as.data.frame(matrix(data = 0, nrow = nrow(colData(currobj)), ncol = length(colcategories)))
                    cores <- vector("list", length = length(colcategories))

                    for (g in 1:length(colcategories)){
                        hmdf[ , g] <- colData(currobj)[ , which(colnames(colData(currobj)) == colcategories[g])]
                        colnames(hmdf)[g] <- colcategories[g]

                        #Test if variable can be coerced to numeric
                        if (!(can_be_made_numeric(hmdf[ , g], cats_to_ignore = class_to_ignore))){
                            if (is.null(cdict)){
                                cores[[g]] <- as.vector(rainbow(length(unique(hmdf[ ,g]))))
                                names(cores[[g]]) <- sort(unique(hmdf[ ,g]))
                                #Use colour table if available
                                if ("ctable" %in% names(metadata(currobj))){
                                    colourshave <- names(cores[[g]])[names(cores[[g]]) %in% rownames(metadata(currobj)$ctable)]
                                    cores[[g]][colourshave]<- metadata(currobj)$ctable[colourshave, "Hex"]
                                }
                            } else {
                                ct <- cdict[[colcategories[g]]]
                                ct <- subset(ct, Name %in% hmdf[ , g])
                                cores[[g]] <- as.vector(ct$Hex)
                                names(cores[[g]]) <- as.vector(ct$Name)
                            }
                        } else {
                            #Variable is (or can be made) numeric, but check for variance in the numbers
                            numvals <- as.numeric(hmdf[, g][which(!(hmdf[, g] %in% class_to_ignore))])
                            hmdf[, g] <- as.numeric(hmdf[, g])
                            if ((max(numvals) - min(numvals)) > 0 ){
                                #If values contain class to ignore, make them black else, only span white and dark blue
                                cores[[g]] <- colorRamp2(c(min(numvals), max(numvals)), c("#3a8aa7", "#780078"), space = "RGB")
                            } else {
                                #Not enough variance, so make them discrete
                                cores[[g]] <- as.vector(rainbow(length(unique(hmdf[, g]))))
                                names(cores[[g]]) <- sort(unique(hmdf[, g]))
                            }
                        }
                        names(cores)[g] <- colcategories[g]
                    }

                    #Make colour scale for relabund heatmap
                    if (hmasPA == FALSE) {
                        if (scaled == TRUE) {
                            #Scale the matrix
                            sampordernames <- colnames(mathm)
                            mathm <- t(apply(mathm, 1, scale))
                            colnames(mathm) <- sampordernames

                            #This is the colour spectrum we are aiming to span
                            PctHmColours <- c("#1307FC", "#FFFFFF", "#F70C00")
                            #Let us see what the distribution looks like to fit it to the colour spectrum
                            #Transform to PPM
                            #RelabundBreakPts <- c(min(mathm), ((min(mathm) + max(mathm)) / 2), max(mathm))
                            RelabundBreakPts <- c(min(mathm), 0, max(mathm))
                            relabundscalename <- "scaling"
                            RelabundBreakPtsLbls <- round(RelabundBreakPts, 2)
                            HMrelabundBreaks <- RelabundBreakPtsLbls
                            relabundheatmapCols <- colorRamp2(HMrelabundBreaks, PctHmColours)
                        } else {
                            #This is the colour spectrum we are aiming to span
                            PctHmColours <- c("blue4", "blue", "slategray1", "khaki", "orange", "tomato", "red", "magenta2", "magenta4")
                            if (analysis %in% taxonomic_spaces){
                                PctBreakPts <- c(0.0001, 0.001, 0.1, 1, 2.5, 5, 10, 50, 100)
                                RelabundBreakPts <- signif(PctBreakPts, digits = 3)
                                relabundscalename <- "Relative Abundance (%)"
                                RelabundBreakPtsLbls <- as.character(paste0(RelabundBreakPts, "%"))
                                HMrelabundBreaks <- Pct2log2PPM(PctBreakPts)
                                relabundheatmapCols <- colorRamp2(HMrelabundBreaks, PctHmColours)
                            } else {
                                if (assay_for_matrix == "GeneCounts"){

                                    if (length(sort(unique(as.vector(mathm)))) == 1){
                                        if (sort(unique(as.vector(mathm))) == 0){
                                            relabundheatmapCols <- "#000000"
                                            RelabundBreakPtsLbls <- "0"
                                            HMrelabundBreaks <- 0
                                        } else {
                                            relabundheatmapCols <- "#FF0000"
                                            RelabundBreakPtsLbls <- as.character(sort(unique(as.vector(mathm))))
                                            HMrelabundBreaks <- as.numeric(sort(unique(as.vector(mathm))))
                                        }
                                    } else {
                                        relabundheatmapCols <- colorRamp2(sort(unique(as.vector(mathm))), c("#000000", rainbow( length(sort(unique(as.vector(mathm)))) - 1 )))
                                        RelabundBreakPtsLbls <- as.character(sort(unique(as.vector(mathm))))
                                        HMrelabundBreaks <- sort(unique(as.vector(mathm)))
                                    }
                                    relabundscalename <- "Number of genes"
                                } else {
                                    #matrix cells are BaseCounts in PPM
                                    #Let us see what the distribution looks like to fit it to the colour spectrum
                                    #Transform to PPM
                                    quantprobs <- round(((2:(length(PctHmColours) - 1)) * (1 / length(PctHmColours))), 1)
                                    #quantprobs <- c(0.10, 0.20, 0.35, 0.50, 0.85, 0.95, 0.99)
                                    nonlogmat <- convert_matrix_log2(mat = mathm, transformation = "from_log2")
                                    #Transform to get quantiles
                                    countmatdistrib <- apply(nonlogmat, function (x) { quantile(x, probs = quantprobs) }, MARGIN = 2)
                                    #take median values of these
                                    medpoints <- rowMedians(countmatdistrib)
                                    PPMrelabundBreaks <- c(0, medpoints, max(nonlogmat))
                                    HMrelabundBreaks <- log2(PPMrelabundBreaks + 1)
                                    RelabundBreakPts <- PPMrelabundBreaks
                                    relabundscalename <- "Parts Per Million"
                                    RelabundBreakPtsLbls <- as.character(paste(RelabundBreakPts, "PPM"))
                                    relabundheatmapCols <- colorRamp2(HMrelabundBreaks, PctHmColours)
                                }
                            }
                        }
                    } else {
                        #Plot cells within heatmap as present/absent
                        relabundheatmapCols <- colorRamp2(c(0, 1), c("blue4", "red"))
                        RelabundBreakPtsLbls <- c("Absent", paste("Present >=", threshPA, "PPM") )
                        HMrelabundBreaks <- c(0, 1)
                        relabundscalename <- "Pres/Abs"
                    }

                    #if (all(c((!is.null(textby)), (fontsizex > 0.2)))){
                    if (!is.null(textby)){
                        hmdf_txt <- as.character(colData(currobj)[ , which(colnames(colData(currobj)) == textby[1])])
                        ha_column <- HeatmapAnnotation(which = "column", df = hmdf, col = cores, textby = anno_text(hmdf_txt, gp = gpar(fontsize = fontsizex, col = "black")), annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 7, col = "black"))
                    } else {
                        ha_column <- HeatmapAnnotation(which = "column", df = hmdf, col = cores, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 7, col = "black"))
                    }

                    #Build plot title
                    nsampmsg <- paste0("Number of samples in heatmap = ", ncol(mathm))
                    msgs <- c(maintit, stattit, presetlist$filtermsg, paste(nsampmsg, ffeatmsg, sep = " | "))
                    plotit <- paste0(msgs, collapse = "\n")

                    if (all(c((any(!(c(is.null(presetlist$minl2fc), is.null(presetlist$maxl2fc))))), ("l2fc" %in% colnames(stathm))))){
                        l2fcmsg <- paste(minl2fcmsg, maxl2fcmsg, sep = " | ")
                        plotit <- paste(plotit, l2fcmsg, sep = "\n")
                    }

                    if (any(c("l2fc", "OddsRatio") %in% colnames(stathm))){
                        plotit <- paste(plotit, binary_directionality_msg, sep = "\n")
                    }

                    #Add plot number if there is more than one heatmap matrix.
                    if (length(matlist) > 1){
                        hmcounter <- paste(hm, length(matlist), sep = "/")
                        plotit <- paste(plotit, paste("Heatmap", hmcounter), sep = "\n")
                    }

                    #Switch plot title if doing a dual heatmap with either genome completeness or percentage from contigs
                    #Get genome completeness hmdf if in taxonomic space
                    ht1fs <- 8
                    hm1tit <- plotit

                    if (all(c(show_GenomeCompleteness, analysis %in% taxonomic_spaces))) {
                        secondaryheatmap <- "GenomeCompleteness"
                    }

                    if (secondaryheatmap == "GenomeCompleteness") {
                        if ("GenomeCompleteness" %in% names(assays(currobj))){
                            gchmdf <- genomecompletenessdf[rownames(mathm), colnames(mathm)]
                            gchmdf <- as.matrix(gchmdf)
                            #gchmdf <- gchmdf * 100
                            #Cap genome completeness to 400% max, to preserve scale across heatmaps
                            gchmdf[which(gchmdf[] > 400)] <- 400
                            if (ncol(mathm) < threshold_for_double_plot){
                                ht1fs <- 7
                                dualHMtit <- plotit
                                hm1tit <- "Relative Abundance"
                            }
                        }
                    } else if (secondaryheatmap == "PctFromCtgs"){
                        if ("PctFromCtgs" %in% names(assays(currobj))){
                            pctctgsdf <- as.data.frame(assays(currobj)$PctFromCtgs)
                            gchmdf <- pctctgsdf[rownames(mathm), colnames(mathm)]
                            gchmdf <- as.matrix(gchmdf)
                            if (ncol(mathm) < threshold_for_double_plot){
                                ht1fs <- 7
                                dualHMtit <- plotit
                                hm1tit <- "Relative Abundance"
                            }
                        }
                    }

                    #Include row annotations for pvalues and l2fc if required.
                    if (("l2fc" %in% colnames(stathm))) {
                        #Get l2fc amplitude for the current heatmap
                        l2fcamplitude <- round(stathm$l2fc, 1)
                        maximuml2fctoshow = 15
                        #Take into account that l2fc may contain only infinites
                        l2fcamplitudeculled <- l2fcamplitude
                        l2fcamplitudeculled[which(l2fcamplitudeculled == Inf)] <- maximuml2fctoshow
                        l2fcamplitudeculled[which(l2fcamplitudeculled == -Inf)] <- (maximuml2fctoshow * -1)
                        l2fcamplitudeshifted <- l2fcamplitudeculled + 15 #Offset amplitude to contain only positive values, maintaining their relativity
                    }

                    #Include annotations for Odds Ratio if PA
                    if ("OddsRatio" %in% colnames(stathm)){
                        ORamplitude <- paste0(round(stathm$OddsRatio, 2), paste0("(", paste(round(stathm$OR95lwr, 2), round(stathm$OR95upr, 2), sep = "-"), ")"))
                    }

                    if (showl2fc == TRUE){
                        showl2fc = "text"
                    }

                    #This is not the most eficient way of doing it, but will keep it like this for now as it is working.
                    if (all(c(("pval" %in% colnames(stathm)), showpval))){
                        if (showroundedpval == TRUE) {
                            statannotnonadj <- as.character(format(round(stathm[, "padj_none"], 3), nsmall = 3))
                            statannotadj <- as.character(format(round(stathm[, "padj_fdr"], 3), nsmall = 3))
                        } else {
                            statannotnonadj <- as.character(signif(stathm[, "padj_none"], digits = 3))
                            statannotadj <- as.character(signif(stathm[, "padj_fdr"], digits = 3))
                        }
                        statannot <- paste(statannotnonadj, statannotadj, sep = " | ")

                        #If showing Log2FC then do so
                        if (all(c(("l2fc" %in% colnames(stathm)), (showl2fc %in% c("text", "plot"))))) {
                            if (showl2fc == "text"){
                                #Show pvalue AND l2fc
                                row_ha <- HeatmapAnnotation(which = "row", Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)), Log2FC = anno_text(l2fcamplitude, gp = gpar(fontsize = fontsizey)), gap = unit(3, "mm"))
                            } else {
                                row_ha <- HeatmapAnnotation(which = "row", Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)), Log2FC = anno_points(l2fcamplitudeshifted, ylim = c(0, 30), width = unit(0.8, "cm"), axis_param = list(side = "bottom", at = c(0, 15, 30), labels = c("<-15", "0", ">15"), labels_rot = 90)), annotation_name_gp = gpar(fontsize = 6, col = "black"))
                            }
                        } else if ("OddsRatio" %in% colnames(stathm)) {
                            #Show Pval and Odds Ratio
                            row_ha <- rowAnnotation(Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)), OR = anno_text(ORamplitude, gp = gpar(fontsize = fontsizey)))
                        } else {
                            #Show only pval
                            row_ha <- rowAnnotation(Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)))
                        }

                    } else {
                        #OK, not showing p-values
                        if (all(c(("l2fc" %in% colnames(stathm)), (showl2fc %in% c("text", "plot"))))) {
                            #If showing Log2FC then do so
                            if (showl2fc == "text"){
                                row_ha <- rowAnnotation(Log2FC = anno_text(l2fcamplitude, gp = gpar(fontsize = fontsizey)))
                            } else {
                                row_ha <- rowAnnotation(Log2FC = anno_points(l2fcamplitudeshifted, ylim = c(0, 30), width = unit(0.8, "cm"), axis_param = list(side = "bottom", at = c(0, 15, 30), labels = c("<-15", "0", ">15"), labels_rot = 90)), annotation_name_gp = gpar(fontsize = 6, col = "black"))
                            }
                        } else if ("OddsRatio" %in% colnames(stathm)) {
                            #Show only Odds Ratio
                            row_ha <- rowAnnotation(OR = anno_text(ORamplitude, gp = gpar(fontsize = fontsizey)))
                        } else if ("Rank" %in% colnames(stathm)){
                            #Include variance rank if Present
                            VarRank <- as.character(stathm$Rank)
                            row_ha <- rowAnnotation(Rank = anno_text(VarRank, gp = gpar(fontsize = fontsizey)))
                        } else {
                            row_ha <- NULL
                        }
                    }

                    #Plot Gram and phyla, if applicable
                    if (all(c(showGram, (analysisname %in% c(taxonomic_spaces, "Species", "Genus", "Family", "Order", "Class"))))){
                        data(Gram)#only for colours
                        tt <- as.data.frame(rowData(currobj))
                        #Make backwards compatible
                        if (!("Gram" %in% colnames(tt))){
                            #There is no Gram information on SummarizedExperiment feature table
                            tt <- tt[rownames(mathm), c(analysisname, "Phylum")]
                            tt <- left_join(tt, Gram, by = "Phylum")
                        } else {
                            #Gram information is present on SummarizedExperiment feature table
                            tt <- tt[rownames(mathm), c(analysisname, "Phylum", "Gram")]
                        }
                        #Fill in missing colours
                        tt <- left_join(tt, Gram[ , c("Phylum", "PhylumColour")], by = "Phylum")
                        tt[is.na(tt$PhylumColour), "PhylumColour"] <- "#BCC2C2"
                        tt$GramColour <- "#BCC2C2"
                        tt$GramColour[which(tt$Gram == "negative")] <- "#FC0345"
                        tt$GramColour[which(tt$Gram == "positive")] <- "#7D00C4"
                        phycols <- setNames(as.character(tt[!duplicated(tt$Phylum), "PhylumColour"]), as.character(tt[!duplicated(tt$Phylum), "Phylum"]))
                        gramcols <- setNames(as.character(tt[!duplicated(tt$Gram), "GramColour"]), as.character(tt[!duplicated(tt$Gram), "Gram"]))
                        hatax <- rowAnnotation(Phylum = tt$Phylum, Gram = tt$Gram, col = list(Phylum = phycols, Gram = gramcols),  annotation_name_gp = gpar(fontsize = 6, col = "black"), show_legend = TRUE)
                    } else {
                        hatax <- NULL
                    }

                    #Determine column order explicitly if required and draw heatmap
                    if (!is.null(ordercolsby)) {
                        if (can_be_made_numeric(colData(currobj)[ , which(colnames(colData(currobj)) == ordercolsby)])){
                            column_order <- order(as.numeric(colData(currobj)[ , which(colnames(colData(currobj)) == ordercolsby)]))
                        } else {
                            column_order <- order(colData(currobj)[ , which(colnames(colData(currobj)) == ordercolsby)])
                        }
                        cluster_column_slices <- FALSE
                        fullheatmap_column_dend <- FALSE
                    } else {
                        column_order <- NULL
                    }

                    if (all(c(no_underscores, (analysis %in% taxonomic_spaces)))) {
                        #rownames(mathm) <- gsub("__", " ", rownames(mathm))
                        rownames(mathm) <- gsub("_", " ", rownames(mathm))
                    }
                    if (analysis %in% taxonomic_spaces){
                        hmfontface <- "italic"
                    } else {
                        hmfontface <- "plain"
                    }
                    #rownames(mathm) <- strtrim(rownames(mathm), 60)

                    #Define if groups will be split
                    if (!(is.null(splitcolsby))){

                        splitcol <- as.character(colData(currobj)[ , which(colnames(colData(currobj)) == splitcolsby)])

                        if (!is.null(column_split_group_order)) {
                            colgroups <- column_split_group_order[column_split_group_order %in% unique(splitcol)]
                            cluster_column_slices <- FALSE
                        } else {
                            colgroups <- sort(unique(splitcol))
                        }

                        column_split <- factor(splitcol, levels = colgroups)

                        if (is.null(column_order)){
                            column_clustering <- TRUE
                        } else {
                            column_clustering <- FALSE
                        }

                        ht1 <- Heatmap(mathm, name = relabundscalename, cluster_columns = column_clustering, column_split = column_split, column_order = column_order, cluster_column_slices = cluster_column_slices, show_column_dend = TRUE, column_dend_side = "top", column_gap = unit(3, "mm"), column_title = hm1tit, column_title_gp = gpar(fontsize = ht1fs), top_annotation = ha_column, col = relabundheatmapCols, column_names_gp = gpar(fontsize = fontsizex), right_annotation = row_ha, left_annotation = hatax, cluster_rows = cluster_rows, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol, fontface = hmfontface),  row_names_max_width = unit(6, "cm"), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)), ...)

                    } else {

                        column_split <- NULL
                        if (all(is.null(ordercolsby), any(cluster_samples_per_heatmap, !cluster_rows))){

                            ht1 <- Heatmap(mathm, name = relabundscalename, column_title = hm1tit, column_title_gp = gpar(fontsize = ht1fs), top_annotation = ha_column, col = relabundheatmapCols, column_names_gp = gpar(fontsize = fontsizex), column_dend_height = unit(5, "mm"), right_annotation = row_ha, left_annotation = hatax, cluster_rows = cluster_rows, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol, fontface = hmfontface), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 4)), row_names_max_width = unit(6, "cm"), ...)

                        } else {

                            #Coerce column order to the order obtained using the full countmatrix, or explicitly by ordercolsby
                            ht1 <- Heatmap(mathm, name = relabundscalename, column_title = hm1tit, column_order = column_order, column_title_gp = gpar(fontsize = ht1fs), top_annotation = ha_column, col = relabundheatmapCols, column_names_gp = gpar(fontsize = fontsizex), cluster_columns = fullheatmap_column_dend, column_dend_height = unit(5, "mm"), right_annotation = row_ha, left_annotation = hatax, cluster_rows = cluster_rows, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol, fontface = hmfontface), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 4)), row_names_max_width = unit(6, "cm"), ...)

                        }
                    }

                    #Make a genome completeness heatmap if in taxonomic space
                    if (all(c((secondaryheatmap %in% c("GenomeCompleteness", "PctFromCtgs")), (any(c(("GenomeCompleteness" %in% names(assays(currobj))), ("PctFromCtgs" %in% names(assays(currobj))))))))) {

                        if (ncol(mathm) < threshold_for_double_plot){
                            SHM_ha_column <- HeatmapAnnotation(df = hmdf, col = cores, show_annotation_name = FALSE)
                        } else {
                            SHM_ha_column <- ha_column
                        }

                        RelabundRowOrder <- suppressWarnings(row_order(ht1))
                        HT1ColumnOrder <- suppressWarnings(column_order(ht1))

                        if (secondaryheatmap == "GenomeCompleteness"){
                            #Draw heatmap with completeness
                            GCheatmapCols <- colorRamp2(c(0, 100, 200, 300, 400), c("white", "forestgreen", "blue", "firebrick1", "black"))
                            ht2 <- Heatmap(gchmdf, name = "GenComp", column_split = column_split, column_title = "% Genome completeness", column_title_gp = gpar(fontsize = ht1fs), top_annotation = SHM_ha_column, col = GCheatmapCols, column_names_gp = gpar(fontsize = fontsizex), right_annotation = NULL, left_annotation = hatax, cluster_rows = FALSE, column_order = unlist(HT1ColumnOrder), row_order = RelabundRowOrder, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol, fontface = hmfontface), heatmap_legend_param = list(direction = "horizontal", title = "% GenComp", labels = c("0%", "100%", "200%", "300%", "> 400%"), title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 4)), row_names_max_width = unit(6, "cm"), ...)
                        } else if (secondaryheatmap == "PctFromCtgs"){
                            #Draw heatmap with percentage from contigs
                            GCheatmapCols <- colorRamp2(c(0, 100), c("white", "midnightblue"))

                            ht2 <- Heatmap(gchmdf, name = "PctFromCtgs", column_split = column_split, column_title = "% Taxonomic info from Contigs", column_title_gp = gpar(fontsize = ht1fs), top_annotation = SHM_ha_column, col = GCheatmapCols, column_names_gp = gpar(fontsize = fontsizex), right_annotation = NULL, left_annotation = hatax, cluster_rows = FALSE, column_order = unlist(HT1ColumnOrder), row_order = RelabundRowOrder, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol, fontface = hmfontface), heatmap_legend_param = list(direction = "horizontal", title = "PctFromCtgs", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 4)), row_names_max_width = unit(6, "cm"), ...)
                        }

                        #Plot heatmaps side by side if there are fewer than the threshold number of samples or less. Else plot one on each page.
                        if (ncol(mathm) < threshold_for_double_plot){
                            ht_list = ht1 + ht2
                            draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "right", ht_gap = unit(0.2, "cm"),padding = unit(c(2, 2, 2, 2), "mm"), column_title = dualHMtit, column_title_gp = gpar(fontsize = ht1fs))
                            drawseparateGChm <- FALSE
                        } else {
                            par(oma = c(10, 7, 3, 10) + 0.1, xpd = TRUE)
                            draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
                            drawseparateGChm <- TRUE
                        }

                    } else {
                        #Draw the heatmap
                        par(oma = c(10, 7, 3, 10) + 0.1, xpd = TRUE)
                        draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
                        drawseparateGChm <- FALSE
                    }

                    #Print what the column annotations are
                    #for(an in colnames(hmdf)) {
                    #    decorate_annotation(an, {
                    #        grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp = gpar(fontsize = 10, col = "black"))
                    #    } )
                    #}
                    if (all(c(("pval" %in% colnames(stathm)), showpval))) {
                        decorate_annotation("Pval", {
                            grid.text("p | FDR p", y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 4, col = "black"))
                        })
                    }
                    if (all(c(("l2fc" %in% colnames(stathm)), (showl2fc %in% c("text", "plot", TRUE))))) {
                        decorate_annotation("Log2FC", {
                            grid.text(paste(fun_for_l2fc, "Log2FC", sep = "\n"), y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 4, col = "black"))
                        })
                    }
                    if ("OddsRatio" %in% colnames(stathm)) {
                        decorate_annotation("OR", {
                            grid.text("Odds Ratio", y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 5, col = "black"))
                        })
                    }
                    if ("Rank" %in% colnames(stathm)) {
                        decorate_annotation("Rank", {
                            grid.text("Rank", y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 5, col = "black"))
                        })
                    }

                    #gvec[[n]]<-recordPlot()
                    n <- n + 1

                    if (drawseparateGChm == TRUE){
                        draw(ht2, heatmap_legend_side = "bottom", annotation_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
                        n <- n + 1
                    }

                    if(n > 100){
                        stop("There are too many combinations. I think you have had enough plots. I am stopping here.")
                    }
                }
            } else {
                #There is no valid matrix, so print out the conditions which led the matrix to be empty.
                plotit <- c(stattit, paste(presetlist$filtermsg, sep=", "))
                plot.new()
                grid.table(c(maintit, "No features fulfilling", plotit), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 10))
                #gvec[[n]]<-recordPlot()
                n <- n + 1
            } #End conditional if there are heatmaps to plot from stats matrix

        } else {
            flog.warn(skip_HM_reason)
            n <- n + 1
        } #End conditional if there is more than a two samples and or two features in subset
    } #End subsetby loop

    #Redefine stats list as ones only containing data
    svec <- svec[sapply(svec, function(x){ !(is.null(x)) } )]

    if (returnstats == TRUE){
        return(svec)
    } else {
        flog.info("Heatmaps generated.")
    }
}
