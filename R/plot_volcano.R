#' Plots volcano plots of differential relative abundance between two classes of a metadata variable, using as input a JAMS-style SummarizedExperiment object.
#'
#' plot_volcano is the volcano-plot sibling of plot_relabund_heatmap. It shares the *entire*
#' upstream pipeline with that function (object vetting, agglomeration, sample and feature
#' selection, subsetting, filtration presets, normalization, and the statistics computed by
#' calculate_matrix_stats), and differs only in how the result is drawn: instead of a
#' ComplexHeatmap, each comparison is rendered as a ggplot2 volcano with -log10(p) on the
#' Y axis and the log2 fold change on the X axis.
#'
#' There are two modes, chosen automatically from the nature of compareby. When compareby
#' is categorical, each pair of classes is contrasted with a Mann-Whitney-Wilcoxon U-test and
#' the X axis is the log2 fold change. When compareby is continuous, every feature is instead
#' correlated against it (Spearman by default) and the X axis is the correlation coefficient;
#' the Y axis, the thresholds, the labelling, the refusal pages and everything else are shared.
#' This mirrors plot_relabund_features, which likewise switches to Spearman for a continuous
#' compareby. To compare High against Low groups of a continuous variable instead of
#' correlating against it, set dichotomize_continuous_compareby_threshold.
#'
#' Because a volcano plot is, by construction, a two-group figure, this function is far
#' stricter than plot_relabund_heatmap about what it will agree to draw. Every refusal is
#' both logged (flog.warn) and, by default, written as a full page into the current graphics
#' device, so that a PDF produced in an unattended loop always records *why* a given
#' comparison is absent rather than silently skipping it. Note that this means an
#' interactive call which refuses will open a device and draw a text page; set
#' report_failures_on_device = FALSE if you would rather it only logged.
#'
#' The one deliberate departure from plot_relabund_heatmap is that nothing is ever discarded
#' by effect size. A minl2fc arriving from applyfilters, or passed explicitly, becomes the
#' position of the vertical significance lines rather than a filter, because removing
#' features by fold change would destroy the very shape a volcano exists to show, and would
#' silently change the multiple-testing denominator along with it.

#' @param ExpObj JAMS-style SummarizedExperiment object.

#' @param glomby String giving the taxonomic level at which to agglomerate counts. This argument should only be used with taxonomic SummarizedExperiment objects. When NULL (the default), there is no agglomeration.

#' @param compareby String specifying the metadata variable name defining the two groups being contrasted. This argument is **mandatory**: unlike plot_relabund_heatmap, there is no "exploratory" mode, because a volcano plot without a comparison is meaningless. If compareby is NULL, missing from the metadata, or resolves to fewer than two usable classes, nothing is plotted and the reason is reported.

#' @param classes_to_compare Optional character vector of exactly two class names within compareby. When supplied, only those two classes are contrasted and all other samples are dropped from the comparison. When NULL (the default) and compareby holds exactly two classes, those two are used; when compareby holds three or more classes, every pairwise combination is plotted on its own page (see max_pairwise_comparisons). Note that the order in which you pass the two names does not set the direction: as in plot_relabund_heatmap, the two classes are sorted alphabetically and a positive fold change means increased in the first of them. Use invertbinaryorder to flip it. Has no meaning in correlation mode and is reported as ignored there.

#' @param dichotomize_continuous_compareby_threshold Controls what happens when compareby holds continuous (numeric) data, for which a fold change is undefined. Accepts NULL (the default, which switches the function into correlation mode: see corr_method), FALSE (treat each distinct numeric value as its own category, which is what you want for numerically-coded groups such as timepoints 0/7/14, bearing in mind that such classes are then ordered as text, so 0, 14, 7), a single number (the cut point), the name of a function such as "median" or "mean", or a function returning a single number. When a cut point is resolved, samples are split into High (>= threshold) and Low (< threshold) using dichotomize_variable_in_md, exactly as in plot_Ordination. The threshold is computed **once**, on the whole vetted object before any subsetting, so that High and Low mean the same thing on every page. The input object is never modified.

#' @param subsetby String (or vector of up to three strings) specifying metadata variable name(s) for subsetting samples. If passed, the full set of comparisons is repeated within each subset. If NULL, data is not subset. Default is NULL.

#' @param samplesToKeep Vector with sample names to keep. If NULL, all samples within the SummarizedExperiment object are kept. Default is NULL. Note that aggressive sample selection is the single most common cause of an unplottable volcano; see min_samples_per_class and refuse_if_underpowered.

#' @param featuresToKeep Vector with feature names to keep. If NULL, all features within the SummarizedExperiment object are kept. Default is NULL. Please note that when agglomerating features with the glomby argument (see above), feature names passed to featuresToKeep must be post-agglomeration feature names. Restricting to a handful of features also shrinks the multiple-testing denominator, which changes every FDR-adjusted p-value relative to an unrestricted run; this is reported on the plot.

#' @param only_allow_CSBs Requires a logical value. If set to TRUE, and the ExpObj passed as input is a taxonomic SummarizedExperiment object of JAMS analysis class "ConsolidatedGenomeBin", the only features considered are Consolidated Species Bins (CSBs). Default is FALSE.

#' @param highlight_features Optional vector of feature names which should always be ringed in black on the volcano, and also named when label_features is TRUE, whether or not they reach the significance thresholds. They do not consume any of the max_features_to_label slots. Useful for tracking a handful of taxa or genes of interest across many comparisons. Names are matched after agglomeration, and for non-taxonomic objects either the bare Accession or the "Accession-Description" form is accepted. Default is NULL.

#' @param wilcox_paired_by String specifying the metadata variable name for pairing samples. Has no meaning in correlation mode and is reported as ignored there. If passed and not NULL, the Mann-Whitney-Wilcoxon U-test will be paired, with samples paired by duplicate values in the wilcox_paired_by column. If the samples cannot be unambiguously paired, the test falls back to unpaired mode with a warning. Default is NULL.

#' @param invertbinaryorder Requires a logical value. If set to TRUE, the sign of the log2 fold change is inverted, which mirrors the volcano left-to-right and swaps which class counts as "up". Default is FALSE. Has no meaning in correlation mode, where the sign of the coefficient is a property of the data, and is reported as ignored there.

#' @param applyfilters Optional string specifying filtration setting "combos", used as a shorthand for setting featcutoff and GenomeCompletenessCutoff. Permissible values are "light", "moderate" or "stringent". Default is "light". Note an important difference from plot_relabund_heatmap: the minl2fc implied by a preset is **not** used to discard features here, because removing features by effect size would gut the very shape a volcano exists to show. It is instead used as the *default* value of l2fc_threshold, i.e. as the vertical significance cut-off. Set l2fc_threshold explicitly to override.

#' @param featcutoff Requires a numeric vector of length 2 for filtering out features by relative abundance: minimum PPM, and the percentage of samples which must reach it. Default is NULL, which defers to applyfilters. As in plot_relabund_heatmap, when subsetby is used the cut-off is applied within the whole object to define a common feature universe across subsets.

#' @param GenomeCompletenessCutoff Requires a numeric vector of length 2 for filtering out taxonomic features by genome completeness: minimum percentage completeness, and the percentage of samples which must reach it. Ignored for non-taxonomic objects. Default is NULL, which defers to applyfilters.

#' @param use_checkM2_style_GenomeCompleteness Requires a logical value. If TRUE (the default), the genome completeness banked alongside the statistics (see returnstats) is the CheckM2-style completeness on its own. If FALSE, the contamination is added to it, as older JAMS reporting did. Only applicable to taxonomic objects which carry a GenomeCompleteness assay, and it affects only the banked statistics, never the plot.

#' @param discard_SDoverMean_below Numeric value setting the minimum standard deviation over mean (sd/mean) for a feature to be kept. When NULL (the default), this filtration is not applied. As in plot_relabund_heatmap, it is applied when the common feature universe is established over the whole object, not again within each subset, and it has no effect when assay_for_matrix = "GeneCounts".

#' @param pval_threshold Numeric value setting the p-value at or below which a feature is considered significant, i.e. the height of the horizontal dashed line. Default is 0.05. Whether this is applied to the raw or to the FDR-adjusted p-value is governed by adj_pval_for_threshold.

#' @param adj_pval_for_threshold Requires TRUE, FALSE or the string "auto". Governs **both** which p-value is plotted on the Y axis and which one pval_threshold is applied to, so that the dashed line on the figure always means what it says. When "auto", the proportion of features still below 0.05 after FDR adjustment is compared with the proportion below 0.05 before adjustment, and the adjusted p-value is used if that ratio exceeds 0.2 (the same rule plot_relabund_heatmap uses). Default is FALSE, i.e. plot and threshold the raw p-value. Whichever is chosen, the count of features passing the *other* one is still reported in the subtitle.

#' @param l2fc_threshold Numeric value, used in binary mode only, setting the minimum absolute log2 fold change for a feature to be considered significant, i.e. the position of the two vertical dashed lines. When NULL (the default), it is taken from the minl2fc implied by applyfilters (1 for "light" and "moderate", 2 for taxonomic / 2.5 for functional objects under "stringent"), falling back to 1 when no preset applies. If a call written for plot_relabund_heatmap passes minl2fc instead, that value is picked up and used here, with a warning. Set to 0 to colour on the p-value alone, except that features whose effect size is exactly zero always stay grey.

#' @param corr_method String naming the correlation test used when compareby is continuous. Permissible values are "spearman" (the default) and "pearson". Spearman is the sensible default for abundance data: it is rank based, as the Mann-Whitney-Wilcoxon test used in binary mode is, and it does not assume normality of zero-inflated relative abundances. It is also what plot_correlation_heatmap defaults to. Ignored when compareby is categorical, or when dichotomize_continuous_compareby_threshold has turned a continuous variable into High/Low groups.

#' @param corr_threshold Numeric value between 0 and 1, used in correlation mode only, setting the minimum absolute correlation coefficient for a feature to be considered significant, i.e. the position of the two vertical dashed lines. When NULL (the default), it is taken from the minabscorrcoeff implied by applyfilters (0.4 for "light", 0.6 for "moderate", 0.8 for "stringent"), falling back to 0.4 when no preset applies. If a call written for another JAMS function passes minabscorrcoeff instead, that value is picked up and used here, with a warning. As with l2fc_threshold, features below it are never discarded, only drawn in grey.

#' @param label_features Requires a logical value. If TRUE (the default), feature names are drawn next to the dots which pass the labelling thresholds, using ggrepel so that labels do not overlap each other or the points.

#' @param label_pval_threshold Numeric value setting the p-value at or below which a feature becomes a candidate for being named on the plot. When NULL (the default), pval_threshold is used. Set this lower than pval_threshold to colour generously but name sparingly.

#' @param label_l2fc_threshold Numeric value setting the minimum absolute log2 fold change for a feature to become a candidate for being named on the plot, in binary mode. When NULL (the default), l2fc_threshold is used.

#' @param label_corr_threshold Numeric value setting the minimum absolute correlation coefficient for a feature to become a candidate for being named on the plot, in correlation mode. When NULL (the default), corr_threshold is used.

#' @param max_features_to_label Numeric value capping how many feature names are drawn. When more features qualify, the ones shown are those with the largest product of -log10(p) and absolute log2 fold change, i.e. the points furthest out into the upper corners of the volcano. Features named in highlight_features are always kept and do not count against the cap. Default is 30.

#' @param label_fontsize Numeric value giving the ggplot2 size (in mm) of the feature labels. When NULL (the default), a size between about 1.3 and 3.5 is computed from how many labels are being drawn, and shrunk by a further tenth when labels wrap over more than two lines, so that a crowded plot automatically gets smaller text.

#' @param label_wrap_width Numeric value giving the number of characters after which a feature name is wrapped onto a new line. Wrapping breaks on spaces where there are any, and hard-breaks over-long unbroken tokens (such as taxonomic names) otherwise. Default is 24.

#' @param max_label_nchar Numeric value giving the number of characters beyond which a feature name is truncated (with a trailing ellipsis) before wrapping. This stops a single MetaCyc-style description from swallowing the panel. Default is 80.

#' @param no_underscores Requires a logical value. If set to TRUE, removes underscores from taxonomic feature names on the plot, so that "LKT__s__Staphylococcus_aureus" is drawn as "LKT s Staphylococcus aureus". Default is FALSE.

#' @param point_size Numeric value giving the size of the dots. When NULL (the default), a size between 0.6 and 2.2 is computed from the number of features being plotted.

#' @param point_alpha Numeric value between 0 and 1 giving the opacity of the dots. Default is 0.7, which keeps dense clouds legible.

#' @param use_class_colours Requires a logical value. Applies in binary mode only, since correlation mode has no classes to look up and always uses the fixed red/blue pair (or up_colour / down_colour). If TRUE (the default), the "up" and "down" dots are coloured with the colours already assigned to those two class names in the object's colour table (metadata(ExpObj)$ctable) or in a cdict passed to the cdict argument, so that a volcano matches the colours used for the same groups in every other JAMS plot. If either class has no colour, or both resolve to the same colour, a fixed red/blue pair is used instead and this is noted in the log.

#' @param cdict Optional named list of per-variable colour dictionaries (as returned by make_colour_dictionary with object_to_return = "cdict"). Consulted before the object's ctable. Default is NULL.

#' @param up_colour Optional colour for the features on the positive side of the X axis: in binary mode those increased in the class named by the directionality message (the alphabetically first of the two unless invertbinaryorder = TRUE), and in correlation mode those positively correlated with compareby. Overrides use_class_colours for that side only. Default is NULL.

#' @param down_colour Optional colour for the features on the negative side of the X axis. Overrides use_class_colours for that side only. Default is NULL.

#' @param nonsignificant_colour Colour for features which do not pass both thresholds. Default is "grey75".

#' @param show_threshold_lines Requires a logical value. If TRUE (the default), dashed guide lines are drawn at the p-value and effect-size thresholds, together with a faint solid line at zero effect size.

#' @param cap_l2fc_at Numeric value at which to clamp the X axis. Infinite log2 fold changes (which arise when a feature is wholly absent from one group and fun_for_l2fc is not "geom_mean") are drawn at the cap as triangles rather than being dropped. When NULL (the default), the cap is the largest finite value present widened by 5 percent, and it is never allowed to fall inside the threshold lines. Ignored in correlation mode, where the X axis is fixed to the full [-1, 1] range of the coefficient so that every correlation volcano is directly comparable with every other.

#' @param cap_neglog10p_at Numeric value at which to clamp the Y axis, used in the same way for p-values that underflow to zero. When NULL (the default), the cap is the largest finite value present widened by 5 percent, and it is never allowed to fall below the horizontal threshold line.

#' @param fun_for_l2fc String specifying the function used to aggregate samples within each group when computing the log2 fold change. Permissible entries are "sum", "mean", "median" or "geom_mean". Default is "geom_mean", for which the geometric mean within each group is exp(mean(log(x + 1))) and which therefore never produces an infinite fold change. Ignored when normalization = "clr", where the effect size is the difference between the arithmetic group means of the centred log-ratio (see normalization), and ignored in correlation mode, which has no groups to aggregate.

#' @param assay_for_matrix String specifying the SummarizedExperiment assay to use. Permissible values are "BaseCounts" (the default, normalized according to normalization) and "GeneCounts" (only available in non-taxonomic objects, and never normalized).

#' @param normalization String specifying how BaseCounts should be normalized. Permissible values are "relabund" (Parts Per Million, the default) and "clr" (centred log-ratio; "compositions" is accepted as a synonym). Ignored when assay_for_matrix is "GeneCounts". Note that a ratio between two centred log-ratios is not a fold change, so under "clr" the X axis instead shows the difference between the two group means of the CLR, converted to log2 units, and is labelled as such. The p-values are still Mann-Whitney-Wilcoxon, computed on the CLR values, and will therefore not be identical to those obtained under "relabund".

#' @param PPM_normalize_to_bases_sequenced Requires a logical value. When TRUE, non-taxonomic objects are relative-abundance-normalized against the total sequenced bases rather than against the coding fraction. Default is FALSE.

#' @param statsonlog Requires a logical value, retained for signature compatibility with plot_relabund_heatmap. Note that the Mann-Whitney-Wilcoxon U-test is rank based, so a monotonic log2 transformation cannot change its p-values. Default is FALSE.

#' @param ignoreunclassified Requires a logical value. If TRUE (the default), the wholly unclassified feature (for example "LKT__Unclassified", or "<analysis>_none") is dropped before the statistics are computed, so that it does not consume a slot in the multiple-testing correction for a term nobody would report. This matches what plot_relabund_heatmap does in practice.

#' @param drop_features_absent_in_all_samples Requires a logical value. If TRUE, features whose counts are zero in every sample of the comparison are removed before testing. Note that in correlation mode features which are *constant* across all samples (which includes, but is wider than, those that are everywhere zero) are always removed regardless of this setting, because a constant feature has no correlation with anything; the count is reported on the plot. This is arguably more correct, since such features carry no information, but it shrinks the FDR denominator and therefore makes the adjusted p-values differ from those plot_relabund_heatmap would report for the same comparison. Default is FALSE, which preserves heatmap parity.

#' @param numthreads Numeric value accepted for signature compatibility with plot_relabund_heatmap and passed through to calculate_matrix_stats, which does not use it on the Mann-Whitney-Wilcoxon path. Default is 1.

#' @param min_samples_per_class Numeric value giving the minimum number of samples which must be present in each of the two classes for the comparison to be attempted at all. Values below 2 are raised to 2. Default is 2, which is the arithmetic minimum for a Mann-Whitney-Wilcoxon U-test, but note that under the other defaults nothing will actually be drawn below 4 per group, because that is where the test first becomes able to return a p-value of 0.05 at all. See refuse_if_underpowered.

#' @param min_samples_for_correlation Numeric value giving the minimum number of samples carrying a usable numeric value for compareby before a correlation is attempted at all. Values below 3 are raised to 3. Default is 5, which is where the exact two-sided Spearman p-value first becomes able to reach 0.05. Note that a coefficient estimated from fewer than about ten points is very unstable, and the plot says so when n is below ten.

#' @param refuse_if_underpowered Requires a logical value. If TRUE (the default), before computing anything the function works out the smallest two-sided p-value the test could possibly return for the group sizes available: 2 / choose(n1 + n2, n1) for the unpaired Mann-Whitney-Wilcoxon U-test, 2 / 2^n for the paired signed-rank test on n pairs, and, in correlation mode, 2 / n! for the exact two-sided Spearman test on n points, which is why five samples are the minimum there. Pearson has no such floor, because a perfect straight line through three points gives an arbitrarily small p-value; in that case only the min_samples_for_correlation guard and a small-n note on the plot apply. If even that floor lies above pval_threshold then no feature can ever be significant, the volcano would be a uniformly grey cloud, and the comparison is refused with a page explaining the arithmetic and naming the group size that would be needed. The floor is deliberately not inflated for multiple testing, because the smallest attainable Benjamini-Hochberg adjusted p-value is the raw floor itself; FDR adjustment can only make matters worse in practice, never better. Set to FALSE to plot it anyway and read the figure purely as a picture of effect sizes.

#' @param min_features_for_volcano Numeric value giving the minimum number of surviving features below which a volcano is not drawn. Values below 2 are raised to 2. Default is 3.

#' @param max_pairwise_comparisons Numeric value capping how many pairwise volcanoes are drawn when compareby holds three or more classes. Pairs are taken in the order returned by combn on the sorted class names, and a page is written recording that the cap was hit. Default is 10.

#' @param max_total_pages Numeric value capping the total number of pages (volcanoes plus explanatory pages) this single call may commit to the device, so that a careless combination of subsetby and a many-class compareby cannot produce an unmanageable document. Default is 200.

#' @param plot_when_nothing_significant Requires a logical value. If TRUE (the default), a volcano in which no feature passes the thresholds is still drawn, with a prominent note, because the absence of signal is itself a result. If FALSE, a text page is written instead, in the manner of plot_relabund_heatmap.

#' @param report_failures_on_device Requires a logical value. If TRUE (the default), every refusal to plot is written as a full text page to the current graphics device in addition to being logged, so that a PDF built in a loop is self-documenting. Set to FALSE when calling interactively and you would rather not have a device opened.

#' @param addtit Optional string with text to prepend to the plot title. Default is NULL.

#' @param returnstats Requires a logical value. If TRUE, returns a named list holding the full statistics data frame for every comparison for which statistics were actually computed (all features, not only the significant ones; comparisons refused before the test was run are absent, though each is reported on its own page). Each data frame carries everything calculate_matrix_stats produced, and, when the comparison went on to be plotted, also neglog10_plotted_pval, Volcano_direction ("Up"/"Down" in binary mode, "Positive"/"Negative" in correlation mode, "NS", or "Not plotted" for features whose statistic was NA), Volcano_status (the human-readable legend label, such as "Up in Responder"), and Volcano_up_class and Volcano_down_class naming the two sides. List names are analysisname_MWW_compareby_subset_classA_vs_classB in binary mode, and analysisname_SPEA (or _PEAR) _compareby_subset_all_samples in correlation mode. Default is FALSE, in which case the function only draws.

#' @param class_to_ignore String or vector specifying classes whose presence in compareby (or in subsetby) causes a sample to be dropped from that comparison. Default is "N_A".

#' @param ... Further arguments are accepted so that a call written for plot_relabund_heatmap can be redirected here without editing. Three of them are honoured, because they have a direct volcano meaning: minl2fc sets l2fc_threshold, minabscorrcoeff sets corr_threshold, and showonlypbelow sets pval_threshold. Everything else is ignored, and every ignored argument name is named in a warning.

#' @return Invisibly NULL, or a named list of statistics data frames when returnstats is TRUE. Plots are drawn to the current graphics device.

#' @export

plot_volcano <- function(ExpObj = NULL, glomby = NULL, compareby = NULL, classes_to_compare = NULL, dichotomize_continuous_compareby_threshold = NULL, subsetby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, only_allow_CSBs = FALSE, highlight_features = NULL, wilcox_paired_by = NULL, invertbinaryorder = FALSE, applyfilters = "light", featcutoff = NULL, GenomeCompletenessCutoff = NULL, use_checkM2_style_GenomeCompleteness = TRUE, discard_SDoverMean_below = NULL, pval_threshold = 0.05, adj_pval_for_threshold = FALSE, l2fc_threshold = NULL, corr_method = "spearman", corr_threshold = NULL, label_features = TRUE, label_pval_threshold = NULL, label_l2fc_threshold = NULL, label_corr_threshold = NULL, max_features_to_label = 30, label_fontsize = NULL, label_wrap_width = 24, max_label_nchar = 80, no_underscores = FALSE, point_size = NULL, point_alpha = 0.7, use_class_colours = TRUE, cdict = NULL, up_colour = NULL, down_colour = NULL, nonsignificant_colour = "grey75", show_threshold_lines = TRUE, cap_l2fc_at = NULL, cap_neglog10p_at = NULL, fun_for_l2fc = "geom_mean", assay_for_matrix = "BaseCounts", normalization = "relabund", PPM_normalize_to_bases_sequenced = FALSE, statsonlog = FALSE, ignoreunclassified = TRUE, drop_features_absent_in_all_samples = FALSE, numthreads = 1, min_samples_per_class = 2, min_samples_for_correlation = 5, refuse_if_underpowered = TRUE, min_features_for_volcano = 3, max_pairwise_comparisons = 10, max_total_pages = 200, plot_when_nothing_significant = TRUE, report_failures_on_device = TRUE, addtit = NULL, returnstats = FALSE, class_to_ignore = "N_A", ...){

    set.seed(2138)

    #Account for JAMS2 spaces
    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")

    #Fixed fallback colours. These are the same pair plot_Ordination uses for dichotomized groups.
    fallback_up_colour <- "#D7301F"
    fallback_down_colour <- "#2B8CBE"

    #Which of the two statistical models this call will use. Decided once, below, from the
    #nature of compareby and of dichotomize_continuous_compareby_threshold.
    #  "binary"      -> Mann-Whitney-Wilcoxon, X axis = log2 fold change between two classes
    #  "correlation" -> Spearman/Pearson against the continuous variable, X axis = the coefficient
    volcano_mode <- "binary"

    #Collect the stats of every comparison attempted
    svec <- list()
    #Count how many pages have been committed to the device, and how many of those are actual volcanoes
    pages_drawn <- 0
    volcanoes_drawn <- 0

    ####################################################################
    #### Internal helpers                                           ####
    ####################################################################

    #Single place where a refusal is both logged and, optionally, committed to the device.
    #Every early exit in this function goes through here, which is what makes an unattended
    #PDF self-documenting rather than mysteriously short.
    report_refusal <- function(headline = NULL, reasons = NULL, advice = NULL){
        #c() already drops NULL components, so only NA and empty strings need weeding out here
        msgs <- as.character(c(headline, reasons, advice))
        msgs <- msgs[!is.na(msgs)]
        msgs <- msgs[nzchar(msgs)]
        for (m in msgs){
            flog.warn(m)
        }
        if (report_failures_on_device){
            #Wrap long lines so that grid.table does not run off the page
            wrapped <- unlist(lapply(msgs, function(m){ strwrap(m, width = 95) }))
            #Keep the page from overflowing vertically
            if (length(wrapped) > 42){
                wrapped <- c(wrapped[1:41], "[...]")
            }
            basefs <- if (length(wrapped) > 24) 7 else 9
            #grid.newpage() rather than plot.new(), because this function also emits ggplot
            #pages onto the same device and mixing base and grid page breaks can overlay them.
            grid.newpage()
            grid.table(wrapped, rows = NULL, cols = NULL, theme = ttheme_default(base_size = basefs, core = list(fg_params = list(hjust = 0, x = 0.02))))
            pages_drawn <<- pages_drawn + 1
        }

        return(invisible(NULL))
    }

    #Hard-wrap a feature name so that it fits a volcano label. strwrap alone is useless for
    #taxonomic names, which contain no spaces at all, so over-long unbroken tokens are split
    #first and only then re-flowed.
    wrap_feature_label <- function(x = NULL, width = 24, maxchar = 80){
        out <- vapply(as.character(x), function(s){
            if (is.na(s)){
                return(NA_character_)
            }
            if (nchar(s) > maxchar){
                s <- paste0(substr(s, 1, (maxchar - 3)), "...")
            }
            if (nchar(s) <= width){
                return(s)
            }
            toks <- unlist(strsplit(s, split = " ", fixed = TRUE))
            toks <- toks[nzchar(toks)]
            if (length(toks) < 1){
                return(s)
            }
            toks <- unlist(lapply(toks, function(tk){
                if (nchar(tk) <= width){
                    return(tk)
                }
                starts <- seq(1, nchar(tk), by = width)
                return(substring(tk, starts, pmin((starts + width - 1), nchar(tk))))
            }))
            return(paste0(strwrap(paste0(toks, collapse = " "), width = width), collapse = "\n"))
        }, FUN.VALUE = character(1), USE.NAMES = FALSE)

        return(out)
    }

    #Smallest two-sided p-value a Mann-Whitney-Wilcoxon U-test can return for these group
    #sizes with no ties. This is the single most useful sanity check for a volcano, because
    #small samplesToKeep vectors routinely produce comparisons in which no feature could
    #possibly clear 0.05 no matter how clean the separation.
    min_attainable_wilcox_p <- function(n1 = NULL, n2 = NULL, paired = FALSE){
        if (any(is.na(c(n1, n2))) || (n1 < 1) || (n2 < 1)){
            return(1)
        }
        if (paired){
            #Signed rank test on n pairs: smallest two-sided p is 2 / 2^n
            npairs <- min(n1, n2)
            minp <- 2 / (2 ^ npairs)
        } else {
            minp <- 2 / choose((n1 + n2), n1)
        }

        return(min(1, minp))
    }

    #Resolve the colour of a class name from a cdict, then from the object's ctable.
    resolve_class_colour <- function(classname = NULL, variable = NULL, SEobj = NULL){
        thecol <- NA_character_
        if (!is.null(cdict) && !is.null(variable)){
            for (v in variable){
                if (!is.na(thecol)){
                    break
                }
                ct <- cdict[[v]]
                if (!is.null(ct) && all(c("Name") %in% colnames(ct))){
                    colcol <- intersect(c("Hex", "Colour"), colnames(ct))[1]
                    if (!is.na(colcol)){
                        hit <- which(as.character(ct$Name) == classname)
                        if (length(hit) > 0){
                            thecol <- as.character(ct[[colcol]])[hit[1]]
                        }
                    }
                }
            }
        }
        if (is.na(thecol) && !is.null(SEobj) && ("ctable" %in% names(metadata(SEobj)))){
            ctab <- metadata(SEobj)$ctable
            colcol <- intersect(c("Hex", "Colour"), colnames(ctab))[1]
            if (!is.na(colcol)){
                #The ctable is rowname-indexed by class name, but do not assume it
                namevec <- if ("Name" %in% colnames(ctab)) as.character(ctab$Name) else rownames(ctab)
                hit <- which(namevec == classname)
                if (length(hit) > 0){
                    thecol <- as.character(ctab[[colcol]])[hit[1]]
                }
            }
        }
        if (!is.na(thecol)){
            #Reject anything that is not a colour R can actually use
            valid <- tryCatch({ grDevices::col2rgb(thecol); TRUE }, error = function(e){ FALSE })
            if (!valid){
                thecol <- NA_character_
            }
        }

        return(thecol)
    }

    ####################################################################
    #### Argument vetting, before any expensive work                ####
    ####################################################################

    #Coerce the logical switches. A string such as "yes" would otherwise reach if() and abort,
    #and for report_failures_on_device that would turn a graceful refusal into a hard failure.
    #Arguments which default to TRUE keep that default on unreadable input; the rest default to FALSE.
    for (lgl_true in c("label_features", "use_class_colours", "show_threshold_lines", "ignoreunclassified", "refuse_if_underpowered", "plot_when_nothing_significant", "report_failures_on_device", "use_checkM2_style_GenomeCompleteness")){
        assign(lgl_true, !isFALSE(suppressWarnings(as.logical(get(lgl_true))[1])))
    }
    for (lgl_false in c("only_allow_CSBs", "invertbinaryorder", "no_underscores", "PPM_normalize_to_bases_sequenced", "statsonlog", "drop_features_absent_in_all_samples", "returnstats")){
        assign(lgl_false, isTRUE(suppressWarnings(as.logical(get(lgl_false))[1])))
    }

    #### Honour the plot_relabund_heatmap arguments that have a volcano equivalent ####
    #This is what makes a heatmap call redirectable here without editing: minl2fc and
    #showonlypbelow mean something on a volcano, they just mean "where the line goes"
    #rather than "which features are discarded".
    dots <- list(...)
    dotargs <- names(dots)
    if (is.null(dotargs)){
        dotargs <- character(0)
    }
    dotargs <- dotargs[nzchar(dotargs)]

    if (("minl2fc" %in% dotargs) && !is.null(dots$minl2fc)){
        if (is.null(l2fc_threshold)){
            l2fc_threshold <- as.numeric(dots$minl2fc)
            flog.warn(paste0("minl2fc = ", l2fc_threshold, " was passed. On a volcano, features are never discarded by effect size, so this is being used as l2fc_threshold, i.e. the position of the vertical significance lines."))
        } else {
            flog.warn("Both minl2fc and l2fc_threshold were passed. l2fc_threshold wins.")
        }
        dotargs <- setdiff(dotargs, "minl2fc")
    }
    if (("minabscorrcoeff" %in% dotargs) && !is.null(dots$minabscorrcoeff)){
        if (is.null(corr_threshold)){
            corr_threshold <- as.numeric(dots$minabscorrcoeff)
            flog.warn(paste0("minabscorrcoeff = ", corr_threshold, " was passed. On a volcano, features are never discarded by effect size, so this is being used as corr_threshold, i.e. the position of the vertical significance lines in correlation mode."))
        } else {
            flog.warn("Both minabscorrcoeff and corr_threshold were passed. corr_threshold wins.")
        }
        dotargs <- setdiff(dotargs, "minabscorrcoeff")
    }
    if (("showonlypbelow" %in% dotargs) && !is.null(dots$showonlypbelow)){
        pval_threshold <- as.numeric(dots$showonlypbelow)
        flog.warn(paste0("showonlypbelow = ", pval_threshold, " was passed. On a volcano, features are never hidden by p-value, so this is being used as pval_threshold, i.e. the height of the horizontal significance line."))
        dotargs <- setdiff(dotargs, "showonlypbelow")
    }
    if (length(dotargs) > 0){
        flog.warn(paste0("plot_volcano: the following argument(s) have no volcano equivalent and are being IGNORED: ", paste0(dotargs, collapse = ", "), "."))
    }

    if (!(is.numeric(pval_threshold) && (length(pval_threshold) == 1) && is.finite(pval_threshold) && (pval_threshold > 0) && (pval_threshold <= 1))){
        stop("pval_threshold must be a single number greater than 0 and no greater than 1.")
    }
    if (!is.null(l2fc_threshold)){
        if (!(is.numeric(l2fc_threshold) && (length(l2fc_threshold) == 1) && is.finite(l2fc_threshold) && (l2fc_threshold >= 0))){
            stop("l2fc_threshold must be NULL or a single non-negative number.")
        }
    }
    if (!(fun_for_l2fc %in% c("sum", "mean", "median", "geom_mean"))){
        stop("fun_for_l2fc must be one of \"sum\", \"mean\", \"median\" or \"geom_mean\".")
    }
    #Unvetted, a character value here would make max() return a string and turn every
    #subsequent effect-size comparison into a silent lexicographic one, labelling nothing.
    for (lab_thr in c("label_pval_threshold", "label_l2fc_threshold", "label_corr_threshold")){
        lab_val <- get(lab_thr)
        if (!is.null(lab_val) && !(is.numeric(lab_val) && (length(lab_val) == 1) && is.finite(lab_val) && (lab_val >= 0))){
            stop(paste0(lab_thr, " must be NULL or a single non-negative number."))
        }
    }
    if (!(length(corr_method) == 1) || !(corr_method %in% c("spearman", "pearson"))){
        stop("corr_method must be either \"spearman\" or \"pearson\".")
    }
    if (!is.null(corr_threshold)){
        if (!(is.numeric(corr_threshold) && (length(corr_threshold) == 1) && is.finite(corr_threshold) && (corr_threshold >= 0) && (corr_threshold <= 1))){
            stop("corr_threshold must be NULL or a single number between 0 and 1, since it is an absolute correlation coefficient.")
        }
    }
    if (!(identical(adj_pval_for_threshold, TRUE) || identical(adj_pval_for_threshold, FALSE) || identical(adj_pval_for_threshold, "auto"))){
        stop("adj_pval_for_threshold must be TRUE, FALSE or \"auto\".")
    }
    if (!(assay_for_matrix %in% c("BaseCounts", "GeneCounts"))){
        stop("assay_for_matrix must be either \"BaseCounts\" or \"GeneCounts\".")
    }
    #plot_relabund_heatmap documents "compositions" for the CLR transform while filter_experiment expects "clr". Accept both.
    if (identical(normalization, "compositions")){
        normalization <- "clr"
    }
    if (!(normalization %in% c("relabund", "clr"))){
        stop("normalization must be either \"relabund\" or \"clr\" (\"compositions\" is accepted as a synonym of \"clr\").")
    }
    #Clamp the numeric knobs. as.numeric() on rubbish yields NA, which would later reach an
    #if() as a missing value and abort before any explanatory page could be drawn, so fall
    #back to the documented default rather than propagating the NA.
    vet_number <- function(x = NULL, lo = NULL, default = NULL){
        x <- suppressWarnings(as.numeric(x))[1]
        if (length(x) < 1 || is.na(x) || !is.finite(x)){
            x <- default
        }

        return(max(lo, x))
    }
    min_samples_per_class <- vet_number(min_samples_per_class, lo = 2, default = 2)
    min_samples_for_correlation <- vet_number(min_samples_for_correlation, lo = 3, default = 5)
    min_features_for_volcano <- vet_number(min_features_for_volcano, lo = 2, default = 3)
    label_wrap_width <- vet_number(label_wrap_width, lo = 4, default = 24)
    max_label_nchar <- vet_number(max_label_nchar, lo = label_wrap_width, default = 80)
    max_features_to_label <- vet_number(max_features_to_label, lo = 0, default = 30)
    max_total_pages <- vet_number(max_total_pages, lo = 1, default = 200)
    max_pairwise_comparisons <- vet_number(max_pairwise_comparisons, lo = 1, default = 10)
    point_alpha <- min(1, vet_number(point_alpha, lo = 0, default = 0.7))

    #Vet the dichotomization argument the same way plot_Ordination does, up front.
    dichot_spec <- dichotomize_continuous_compareby_threshold
    treat_compareby_as_is <- isFALSE(dichot_spec)
    if (treat_compareby_as_is){
        dichot_spec <- NULL
    }
    if (!is.null(dichot_spec)){
        #A number passed as a string (easy to do, given that all JAMS metadata is character) is a number, not a function name
        if (is.character(dichot_spec) && (length(dichot_spec) == 1) && !is.na(dichot_spec) && can_be_made_numeric(dichot_spec)){
            dichot_spec <- as.numeric(dichot_spec)
        }
        dichot_spec_ok <- is.function(dichot_spec) || (is.numeric(dichot_spec) && (length(dichot_spec) == 1) && is.finite(dichot_spec)) || (is.character(dichot_spec) && (length(dichot_spec) == 1) && !is.na(dichot_spec))
        if (!dichot_spec_ok){
            stop("dichotomize_continuous_compareby_threshold must be NULL, FALSE, a single number, the name of a function such as \"median\", or a function returning a single number.")
        }
    }

    if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
        stop("plot_volcano can only take a SummarizedExperiment object as input.")
    }

    #### compareby is mandatory and must exist ####
    if (is.null(compareby) || (length(compareby) < 1) || all(is.na(compareby)) || !any(nzchar(as.character(compareby)))){
        report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = c("No value was passed to the compareby argument.", "A volcano plot shows the log2 fold change between exactly two groups of samples against the p-value of the difference, so, unlike plot_relabund_heatmap, it has no exploratory mode and cannot be drawn without a comparison."), advice = "Pass the name of a metadata variable holding two classes, for example compareby = \"Treatment\". Run colnames(colData(ExpObj)) to see what is available.")

        return(invisible(NULL))
    }
    if (length(compareby) > 1){
        flog.warn(paste0("plot_volcano: compareby must be a single variable; only the first (\"", compareby[1], "\") will be used."))
        compareby <- compareby[1]
    }
    compareby <- as.character(compareby)

    available_vars <- colnames(colData(ExpObj))
    missing_vars <- character(0)
    for (argname in c("compareby", "subsetby", "wilcox_paired_by")){
        thisvar <- get(argname)
        if (!is.null(thisvar) && is.character(thisvar) && (length(thisvar) > 0)){
            absent <- setdiff(thisvar[!is.na(thisvar)], available_vars)
            if (length(absent) > 0){
                missing_vars <- c(missing_vars, paste0(argname, " = \"", paste0(absent, collapse = "\", \""), "\""))
            }
        }
    }
    if (length(missing_vars) > 0){
        report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = c(paste0("The following metadata variable(s) named in the call are not present in the SummarizedExperiment object: ", paste0(missing_vars, collapse = "; "), "."), paste0("Variables available are: ", paste0(available_vars, collapse = ", "), ".")), advice = "Check the spelling against colnames(colData(ExpObj)). Note that JAMS metadata column names may contain spaces, so a name copied from a make.names()-mangled data frame will not match.")

        return(invisible(NULL))
    }

    #### Pre-flight on samplesToKeep and featuresToKeep ####
    #ExpObjVetting stops() on some of these, which would abort a whole PDF loop. Catch the
    #common, recoverable cases here so they become a page rather than an error.
    if (!is.null(samplesToKeep)){
        samplesToKeep <- unique(as.character(samplesToKeep))
        stk_have <- samplesToKeep[samplesToKeep %in% colnames(ExpObj)]
        if (length(stk_have) < 1){
            report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = c(paste0("None of the ", length(samplesToKeep), " sample name(s) passed to samplesToKeep are present in this object."), paste0("Examples of what was asked for: ", paste0(utils::head(samplesToKeep, 5), collapse = ", "), "."), paste0("Examples of what is present: ", paste0(utils::head(colnames(ExpObj), 5), collapse = ", "), ".")), advice = "Compare your vector against colnames(ExpObj).")

            return(invisible(NULL))
        }
        if (length(stk_have) < min((min_samples_per_class * 2), min_samples_for_correlation)){
            report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = c(paste0("Only ", length(stk_have), " of the ", length(samplesToKeep), " sample(s) passed to samplesToKeep are present in this object."), paste0("A volcano needs at least ", min((min_samples_per_class * 2), min_samples_for_correlation), " samples, so there is nothing to test.")), advice = "Widen samplesToKeep, or lower min_samples_per_class / min_samples_for_correlation if you accept that the test will be badly underpowered.")

            return(invisible(NULL))
        }
        if (length(stk_have) < length(samplesToKeep)){
            flog.warn(paste0("plot_volcano: ", (length(samplesToKeep) - length(stk_have)), " of the ", length(samplesToKeep), " sample(s) passed to samplesToKeep are not present in this object and will be ignored."))
        }
    }
    if (!is.null(featuresToKeep) && is.null(glomby)){
        featuresToKeep <- unique(as.character(featuresToKeep))
        ftk_have <- featuresToKeep[featuresToKeep %in% rownames(ExpObj)]
        if (length(ftk_have) < min_features_for_volcano){
            report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = c(paste0("Only ", length(ftk_have), " of the ", length(featuresToKeep), " feature name(s) passed to featuresToKeep are present in this object."), paste0("A volcano needs at least ", min_features_for_volcano, " features to be worth drawing; with fewer than that, plot_relabund_features is the better tool."), if (length(ftk_have) > 0) paste0("Found: ", paste0(utils::head(ftk_have, 10), collapse = ", "), ".") else NULL), advice = "Check the names against rownames(ExpObj), remembering that featuresToKeep is interpreted in the post-agglomeration namespace when glomby is set.")

            return(invisible(NULL))
        }
    }

    ####################################################################
    #### Vet the object                                             ####
    ####################################################################

    variables_to_fix <- c(compareby, subsetby)

    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, only_allow_CSBs = only_allow_CSBs, glomby = glomby, variables_to_fix = variables_to_fix, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    #This runs before the mode is known, so it can only enforce the looser of the two
    #requirements. The real, mode-specific gates are applied per subset further down.
    min_samples_anymode <- min((min_samples_per_class * 2), min_samples_for_correlation)
    if (ncol(obj) < min_samples_anymode){
        report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = c(paste0("Only ", ncol(obj), " sample(s) remain after applying samplesToKeep and removing samples whose value for \"", compareby, "\" is in class_to_ignore (", paste0(class_to_ignore, collapse = ", "), ")."), paste0("A volcano needs at least ", min_samples_anymode, " samples, whether it is comparing two groups or correlating against a continuous variable.")), advice = "Check that the samples you want are not carrying a class_to_ignore value in the comparison variable.")

        return(invisible(NULL))
    }
    if (nrow(obj) < min_features_for_volcano){
        report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = paste0("Only ", nrow(obj), " feature(s) remain in the object, which is fewer than the ", min_features_for_volcano, " needed to draw a volcano."), advice = "Relax featuresToKeep, or use plot_relabund_features to show a handful of features individually.")

        return(invisible(NULL))
    }

    ####################################################################
    #### Continuous compareby: refuse, or dichotomize into High/Low ####
    ####################################################################
    #JAMS metadata is all-character, so "continuous" is decided by whether the usable values
    #can be read as numbers. This is done ONCE, on the whole vetted object and BEFORE any
    #subsetting, so that the threshold, and therefore the meaning of High and Low, is
    #identical on every page. The input object is never modified.

    dichotomized_compareby <- FALSE
    compareby_original <- compareby
    dichot_labels <- c(High = "High", Low = "Low")
    dichot_msg <- NULL

    cb_vals <- as.character(colData(obj)[[compareby]])
    cb_usable <- cb_vals[!is.na(cb_vals) & !(cb_vals %in% class_to_ignore)]
    if (length(cb_usable) < 1){
        report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = paste0("Comparison variable \"", compareby, "\" has no usable values among the samples being plotted: every value is either NA or in class_to_ignore (", paste0(class_to_ignore, collapse = ", "), ")."), advice = "Check the metadata column, or pass a different class_to_ignore.")

        return(invisible(NULL))
    }
    cb_is_numeric <- can_be_made_numeric(cb_usable)
    cb_n_distinct <- if (cb_is_numeric) length(unique(as.numeric(cb_usable))) else length(unique(cb_usable))

    #Respect an explicit declaration in phenolabels, if there is one for this variable
    cb_declared_type <- NA_character_
    pl <- metadata(obj)$phenolabels
    if (!is.null(pl) && all(c("Var_label", "Var_type") %in% colnames(pl)) && (compareby %in% pl$Var_label)){
        cb_declared_type <- as.character(pl$Var_type[match(compareby, pl$Var_label)])
    }

    cb_is_continuous <- cb_is_numeric && (cb_n_distinct > 2) && (is.na(cb_declared_type) || (cb_declared_type == "continuous"))

    if (cb_is_continuous && treat_compareby_as_is){
        flog.info(paste0("Comparison variable \"", compareby, "\" looks continuous, but dichotomize_continuous_compareby_threshold = FALSE, so each of its ", cb_n_distinct, " distinct values will be treated as a separate category."))

    } else if (cb_is_continuous && is.null(dichot_spec)){
        #### Correlation mode ####
        #There is no fold change to compute against a continuous variable, but there is a
        #perfectly good per-feature correlation, which calculate_matrix_stats already knows how
        #to produce. This is also what plot_relabund_features does with a continuous compareby,
        #so correlating (rather than refusing) is the behaviour consistent with the rest of JAMS.
        volcano_mode <- "correlation"
        flog.info(paste0("Comparison variable \"", compareby, "\" holds continuous data: its ", length(cb_usable), " usable values are all numeric and span ", cb_n_distinct, " distinct values (from ", signif(min(as.numeric(cb_usable)), 4), " to ", signif(max(as.numeric(cb_usable)), 4), ")."))
        flog.info(paste0("Switching to correlation mode: each feature will be correlated against \"", compareby, "\" with a two-sided ", corr_method, " test, and the X axis of the volcano will be the correlation coefficient rather than a log2 fold change. To compare High against Low groups instead, set dichotomize_continuous_compareby_threshold (a number, \"median\", or any function returning one number); to treat the distinct numeric values as categories, set it to FALSE."))

        #Several arguments only mean something for a two-group contrast. Say so rather than
        #letting the user believe they took effect.
        inert <- character(0)
        if (!is.null(classes_to_compare)){
            inert <- c(inert, "classes_to_compare")
            classes_to_compare <- NULL
        }
        if (!is.null(wilcox_paired_by)){
            inert <- c(inert, "wilcox_paired_by")
            wilcox_paired_by <- NULL
        }
        if (isTRUE(invertbinaryorder)){
            inert <- c(inert, "invertbinaryorder")
        }
        if (!identical(fun_for_l2fc, "geom_mean")){
            inert <- c(inert, "fun_for_l2fc")
        }
        if (!is.null(l2fc_threshold)){
            inert <- c(inert, "l2fc_threshold (use corr_threshold in correlation mode)")
        }
        if (length(inert) > 0){
            flog.warn(paste0("In correlation mode the following argument(s) have no meaning and are being ignored: ", paste0(inert, collapse = ", "), "."))
        }

    } else if (cb_is_continuous){
        cb_num <- as.numeric(cb_usable)
        if (is.numeric(dichot_spec)){
            thrinfo <- resolve_dichotomization_threshold(x = cb_num, threshold = dichot_spec)
        } else {
            thrinfo <- resolve_dichotomization_threshold(x = cb_num, threshold = NULL, func_for_threshold = dichot_spec)
        }

        dichot_vals <- dichotomize_variable_in_md(md = data.frame(cbvar = cb_vals, stringsAsFactors = FALSE), column_to_dichotomize = "cbvar", threshold = thrinfo$threshold, class_to_ignore = class_to_ignore, high_label = dichot_labels[["High"]], low_label = dichot_labels[["Low"]])

        dichot_colname <- paste0(compareby, "_dichotomized")
        if (dichot_colname %in% colnames(colData(obj))){
            flog.warn(paste0("A metadata variable called \"", dichot_colname, "\" already exists. It will be overwritten (within this plot only; your SummarizedExperiment object is not modified) by the dichotomized values of \"", compareby, "\"."))
        }
        colData(obj)[[dichot_colname]] <- dichot_vals

        #Defensive: should not happen, as ExpObjVetting already removed class_to_ignore samples for compareby
        dichot_missing <- which(!(dichot_vals %in% dichot_labels))
        if (length(dichot_missing) > 0){
            flog.warn(paste0(length(dichot_missing), " sample(s) could not be assigned to High or Low for \"", compareby, "\" and will be omitted."))
            obj <- obj[ , -dichot_missing]
        }

        n_high <- sum(colData(obj)[[dichot_colname]] == dichot_labels[["High"]])
        n_low <- sum(colData(obj)[[dichot_colname]] == dichot_labels[["Low"]])
        dichot_msg <- paste0(compareby, " dichotomized at ", thrinfo$description, ": High >= ", thrinfo$formatted_threshold, " (n=", n_high, "); Low < ", thrinfo$formatted_threshold, " (n=", n_low, ")")
        flog.info(dichot_msg)

        if (min(n_high, n_low) < min_samples_per_class){
            report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = c(dichot_msg, paste0("Dichotomizing \"", compareby, "\" at ", thrinfo$description, " left ", n_high, " sample(s) in High and ", n_low, " in Low, and at least ", min_samples_per_class, " are needed in each group.")), advice = "Choose a threshold which splits the samples more evenly, or widen samplesToKeep.")

            return(invisible(NULL))
        }

        compareby_original <- compareby
        compareby <- dichot_colname
        dichotomized_compareby <- TRUE
        #classes_to_compare, if given, must now be expressed in High/Low terms
        if (!is.null(classes_to_compare)){
            if (!all(as.character(classes_to_compare) %in% dichot_labels)){
                flog.warn(paste0("classes_to_compare was given as (", paste0(classes_to_compare, collapse = ", "), "), but \"", compareby_original, "\" has been dichotomized into \"", dichot_labels[["High"]], "\" and \"", dichot_labels[["Low"]], "\". Ignoring classes_to_compare."))
                classes_to_compare <- NULL
            }
        }

    } else if (!is.null(dichot_spec)){
        flog.info(paste0("Comparison variable \"", compareby, "\" is categorical; ignoring dichotomize_continuous_compareby_threshold."))
    }

    ####################################################################
    #### Filtration presets                                         ####
    ####################################################################

    #Set applyfilters to null if featuresToKeep is set, exactly as plot_relabund_heatmap does:
    #the user has already said precisely which features they want.
    if (!is.null(featuresToKeep)){
        applyfilters <- NULL
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, is16S = ("asv2lkt" %in% names(metadata(obj))))

    #A volcano must NOT drop features by effect size, or its characteristic shape disappears
    #and the FDR denominator silently changes. The minl2fc (or minabscorrcoeff) implied by a
    #preset is therefore reinterpreted here as the significance cut-off, i.e. where the
    #vertical lines go. Both modes resolve their own threshold the same way.
    if (is.null(l2fc_threshold)){
        if (!is.null(presetlist$minl2fc) && is.finite(presetlist$minl2fc)){
            l2fc_threshold <- as.numeric(presetlist$minl2fc)
            flog.info(paste0("l2fc_threshold was not set; taking it from applyfilters = \"", applyfilters, "\", which implies a minimum absolute log2 fold change of ", round(l2fc_threshold, 2), ". Note that, unlike in plot_relabund_heatmap, features below this are NOT discarded: they are drawn in grey as non-significant."))
        } else {
            l2fc_threshold <- 1
            flog.info("l2fc_threshold was not set and no filtration preset implies one; defaulting to 1.")
        }
    }
    if (is.null(corr_threshold)){
        if (!is.null(presetlist$minabscorrcoeff) && is.finite(presetlist$minabscorrcoeff)){
            corr_threshold <- as.numeric(presetlist$minabscorrcoeff)
            if (volcano_mode == "correlation"){
                flog.info(paste0("corr_threshold was not set; taking it from applyfilters = \"", applyfilters, "\", which implies a minimum absolute correlation coefficient of ", round(corr_threshold, 2), ". As with the fold change, features below this are NOT discarded, only drawn in grey."))
            }
        } else {
            corr_threshold <- 0.4
            if (volcano_mode == "correlation"){
                flog.info("corr_threshold was not set and no filtration preset implies one; defaulting to 0.4, which is the value the \"light\" preset uses for minabscorrcoeff.")
            }
        }
    }

    #From here on, a single pair of names carries whichever effect size this call is using,
    #so that everything downstream of the statistics is shared between the two modes.
    if (volcano_mode == "correlation"){
        effect_threshold <- corr_threshold
        effect_axis_name <- if (corr_method == "spearman") "Spearman rho" else "Pearson r"
        effect_size_name <- effect_axis_name
    } else {
        effect_threshold <- l2fc_threshold
        effect_size_name <- "log2 fold change"
        effect_axis_name <- "log2 fold change"
    }

    if (is.null(label_pval_threshold)){
        label_pval_threshold <- pval_threshold
    }
    if (is.null(label_l2fc_threshold)){
        label_l2fc_threshold <- l2fc_threshold
    }
    if (is.null(label_corr_threshold)){
        label_corr_threshold <- corr_threshold
    }
    label_effect_threshold <- if (volcano_mode == "correlation") label_corr_threshold else label_l2fc_threshold

    if (assay_for_matrix == "GeneCounts"){
        flog.warn("Counts matrix used for the volcano will represent the *number of genes* for each feature, rather than its relative abundance. For using relative abundance (the default), set assay_for_matrix = \"BaseCounts\"")
        hmtypemsg <- if (volcano_mode == "correlation") "Gene count Correlation Volcano Plot" else "Gene count Volcano Plot"
        global_FTK <- rownames(obj)
    } else {
        hmtypemsg <- if (normalization == "clr") "CLR-transformed Volcano Plot" else "Relative Abundance Volcano Plot"
        if (volcano_mode == "correlation"){
            hmtypemsg <- sub("Volcano Plot", "Correlation Volcano Plot", hmtypemsg)
        }

        #Obtain a "global" list of features so that the same feature universe is used within every subset.
        #Guard each component against NULL: `all(NULL != c(0,0))` is TRUE, which would make this always fire.
        needs_global_filter <- any(c((!is.null(applyfilters)), (!is.null(presetlist$featcutoff) && any(presetlist$featcutoff != 0)), (!is.null(presetlist$GenomeCompletenessCutoff) && any(presetlist$GenomeCompletenessCutoff != 0)), (!is.null(discard_SDoverMean_below))))
        if (needs_global_filter){
            currobj <- filter_experiment(SEobj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = NULL, featuresToKeep = NULL, normalization = "relabund", PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, only_allow_CSBs = FALSE, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, discard_SDoverMean_below = discard_SDoverMean_below)
            global_FTK <- rownames(currobj)
        } else {
            global_FTK <- rownames(obj)
        }
    }

    if (length(global_FTK) < min_features_for_volcano){
        report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = c(paste0("Only ", length(global_FTK), " feature(s) survived filtration, which is fewer than the ", min_features_for_volcano, " needed to draw a volcano."), paste0("Filtration applied: ", paste0(presetlist$filtermsg, collapse = " | "))), advice = "Relax applyfilters (try \"light\", or NULL), or lower featcutoff and GenomeCompletenessCutoff.")

        return(invisible(NULL))
    }

    ####################################################################
    #### Subsetting                                                 ####
    ####################################################################

    if (!(is.null(subsetby))){
        #Do NOT hand a continuous compareby to multiple_subsetting_sample_selector: it runs
        #shapiro.test() per subset for the QC column, and shapiro.test ERRORS (it does not warn)
        #when a group has fewer than 3 samples or when every value in it is identical. Both are
        #routine in a tier-2 subset, and the error would abort the whole call before any
        #explanatory page could be drawn. None of the columns plot_volcano consumes from
        #Subsets_stats depend on compareby, so it is simply withheld in correlation mode; the
        #sample counts are then slightly generous (class_to_ignore samples are no longer
        #dropped here), which the per-subset guards further down re-check anyway.
        subset_compareby <- if (volcano_mode == "correlation") NULL else compareby
        subset_list <- multiple_subsetting_sample_selector(SEobj = obj, phenotable = NULL, subsetby = subsetby, compareby = subset_compareby, cats_to_ignore = class_to_ignore)
        if (is.null(subset_list)){
            report_refusal(headline = "NO VOLCANO PLOT WAS DRAWN", reasons = paste0("The variable(s) passed to subsetby (", paste0(subsetby, collapse = ", "), ") could not be used for subsetting, because none of them is a discrete metadata variable."), advice = "Pass a categorical variable to subsetby, or set subsetby = NULL.")

            return(invisible(NULL))
        }
        subset_df <- subset_list$Subsets_stats
        subset_df <- subset_df[which(subset_df$Subset_Tier_Level != 0), , drop = FALSE]
        min_samples_in_subset <- if (volcano_mode == "correlation") min_samples_for_correlation else (min_samples_per_class * 2)
        if (any(subset_df$Num_samples_in_subset < min_samples_in_subset)){
            LowSampSubsets <- subset_df[which(subset_df$Num_samples_in_subset < min_samples_in_subset), "Subset_Tier_Class_Name"]
            flog.warn(paste0("Subsets ", paste0(LowSampSubsets, collapse = ", "), " contain fewer than ", min_samples_in_subset, " samples, and will thus be omitted."))
            subset_df <- subset_df[which(subset_df$Num_samples_in_subset >= min_samples_in_subset), , drop = FALSE]
        }
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

    ####################################################################
    #### Main loop: one or more comparisons within each subset      ####
    ####################################################################

    nothing_was_drawn_reasons <- character(0)

    for (sp in 1:length(subset_points)){

        #The cap has to be tested here as well as inside the pair loop, because the many `next`
        #statements below can each commit a refusal page without ever reaching the pair loop.
        if (pages_drawn >= max_total_pages){
            flog.warn(paste0("plot_volcano: ", pages_drawn, " pages have already been committed to the device, which is the cap set by max_total_pages. Stopping here rather than producing an unmanageable document."))
            break
        }

        if (!(is.null(subsetby))){
            samplesToKeep_sp <- subset_list[[subset_points[sp]]]
            flog.info(paste("Plotting within", subset_points[sp]))
            subsetname <- subset_points[sp]
            subsettit <- paste("within", subset_points[sp])
        } else {
            samplesToKeep_sp <- rownames(colData(obj))
            subsetname <- "no_sub"
            subsettit <- NULL
        }
        #Defend against a subset vector carrying names which are no longer in the object
        samplesToKeep_sp <- samplesToKeep_sp[samplesToKeep_sp %in% colnames(obj)]

        subset_label <- if (is.null(subsetby)) "this dataset" else paste0("subset \"", subsetname, "\"")

        min_samples_here <- if (volcano_mode == "correlation") min_samples_for_correlation else (min_samples_per_class * 2)
        if (length(samplesToKeep_sp) < min_samples_here){
            thisreason <- paste0("Only ", length(samplesToKeep_sp), " sample(s) are available in ", subset_label, ", and at least ", min_samples_here, " are needed for a ", if (volcano_mode == "correlation") "correlation" else "two-group comparison", ".")
            report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(subset_label)), reasons = thisreason, advice = "Widen samplesToKeep or choose a coarser subsetby variable.")
            nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
            next
        }

        curr_pt <- colData(obj)[samplesToKeep_sp, , drop = FALSE]
        cb_here <- as.character(curr_pt[[compareby]])
        names(cb_here) <- samplesToKeep_sp
        usable_here <- !is.na(cb_here) & !(cb_here %in% class_to_ignore)
        classes_here <- sort(unique(cb_here[usable_here]))

        if (volcano_mode == "correlation"){
            #A correlation needs spread in the variable, not classes in it. A variable which is
            #constant within this subset has an undefined correlation with anything.
            cb_num_here <- suppressWarnings(as.numeric(cb_here[usable_here]))
            cb_num_here <- cb_num_here[!is.na(cb_num_here)]
            if (length(cb_num_here) < min_samples_for_correlation){
                thisreason <- paste0("Only ", length(cb_num_here), " sample(s) in ", subset_label, " carry a usable numeric value for \"", compareby, "\", and at least ", min_samples_for_correlation, " are needed to correlate against it.")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(subset_label)), reasons = thisreason, advice = "Widen samplesToKeep, or choose a coarser subsetby variable.")
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }
            if (length(unique(cb_num_here)) < 2){
                thisreason <- paste0("Comparison variable \"", compareby, "\" takes a single value (", signif(cb_num_here[1], 4), ") across every sample in ", subset_label, ", so it has no variance and cannot be correlated with anything.")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(subset_label)), reasons = thisreason, advice = "Choose a subsetby variable which does not hold the comparison variable constant, or correlate against a different variable.")
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }
            #One comparison per subset. The NA is a placeholder so the loop below can be shared.
            pairs_to_plot <- list(NA_character_)

        } else if (length(classes_here) < 2){
            thisreason <- paste0("Comparison variable \"", compareby, "\" resolves to ", length(classes_here), " usable class", if (length(classes_here) == 1) paste0(" (\"", classes_here, "\")") else "es", " within ", subset_label, ", so there is nothing to contrast.")
            report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(subset_label)), reasons = c(thisreason, "A volcano plot needs exactly two groups of samples."), advice = "Use plot_relabund_heatmap with hmtype = \"exploratory\" to rank these features by variance instead, or choose a comparison variable which varies within this subset.")
            nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
            next

        } else if (!is.null(classes_to_compare)){
            ctc <- unique(as.character(classes_to_compare))
            ctc_have <- ctc[ctc %in% classes_here]
            if (length(ctc_have) != 2){
                thisreason <- paste0("classes_to_compare was given as (", paste0(ctc, collapse = ", "), "), but ", length(ctc_have), " of those class(es) are present in ", subset_label, ". Classes present here are: ", paste0(classes_here, collapse = ", "), ".")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(subset_label)), reasons = thisreason, advice = "Pass exactly two class names which both occur in every subset, or set classes_to_compare = NULL to let every pair be plotted automatically.")
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }
            pairs_to_plot <- list(sort(ctc_have))
        } else if (length(classes_here) == 2){
            pairs_to_plot <- list(classes_here)
        } else {
            allpairs <- utils::combn(classes_here, 2, simplify = FALSE)
            if (length(allpairs) > max_pairwise_comparisons){
                report_refusal(headline = paste0("PAIRWISE COMPARISONS CAPPED IN ", toupper(subset_label)), reasons = c(paste0("Comparison variable \"", compareby, "\" holds ", length(classes_here), " classes here: ", paste0(classes_here, collapse = ", "), "."), paste0("That is ", length(allpairs), " pairwise volcano plots, which exceeds max_pairwise_comparisons = ", max_pairwise_comparisons, "."), paste0("Only the first ", max_pairwise_comparisons, " pairs (in alphabetical order) follow.")), advice = "Raise max_pairwise_comparisons, or name the two classes you want with classes_to_compare.")
                allpairs <- utils::head(allpairs, max_pairwise_comparisons)
            } else {
                flog.info(paste0("Comparison variable \"", compareby, "\" has ", length(classes_here), " classes within ", subset_label, ". Drawing all ", length(allpairs), " pairwise volcano plots."))
            }
            pairs_to_plot <- allpairs
        }

        #### Build the counts matrix for this subset, once ####
        if (assay_for_matrix == "BaseCounts"){
            #Only pass a GenomeCompletenessCutoff for taxonomic objects: declare_filtering_presets,
            #which filter_experiment calls again internally, raises a base R warning() otherwise, once per subset.
            gcc_here <- if (analysis %in% taxonomic_spaces) c(0, 0) else NULL
            currobj <- filter_experiment(SEobj = obj, featcutoff = c(0, 0), samplesToKeep = samplesToKeep_sp, featuresToKeep = global_FTK, only_allow_CSBs = FALSE, normalization = normalization, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = gcc_here, discard_SDoverMean_below = NULL, flush_out_empty_samples = FALSE)
            wanted_assay <- if (normalization == "clr") "CLR" else "PPM"
        } else {
            currobj <- obj[global_FTK, samplesToKeep_sp]
            wanted_assay <- "GeneCounts"
        }

        if (!(wanted_assay %in% names(assays(currobj)))){
            thisreason <- paste0("The \"", wanted_assay, "\" assay is not present in this SummarizedExperiment object (assays available: ", paste0(names(assays(currobj)), collapse = ", "), ").")
            report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(subset_label)), reasons = thisreason, advice = if (wanted_assay == "GeneCounts") "GeneCounts only exists in non-taxonomic (functional) objects. Set assay_for_matrix = \"BaseCounts\"." else "Check that the object was built by a current version of JAMSbeta.")
            nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
            next
        }
        countmat <- as.matrix(assays(currobj)[[wanted_assay]])
        #A raw (un-normalized) copy is kept purely so that "absent from every sample" can be
        #judged on counts rather than on a CLR value, for which zero means nothing of the sort.
        can_judge_absence <- ("BaseCounts" %in% names(assays(currobj)))
        if (can_judge_absence){
            rawcountmat <- as.matrix(assays(currobj)[["BaseCounts"]])[rownames(countmat), colnames(countmat), drop = FALSE]
        } else {
            rawcountmat <- countmat
            if (drop_features_absent_in_all_samples){
                flog.warn("drop_features_absent_in_all_samples = TRUE, but this object carries no BaseCounts assay, so absence cannot be judged on raw counts. The option is being ignored for this comparison rather than guessing from normalized values.")
            }
        }

        #Drop the wholly unclassified / unattributed feature before testing, so that it does
        #not occupy a slot in the multiple-testing correction for a term nobody would report.
        if (ignoreunclassified == TRUE){
            dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified", "hypothetical protein", "putative protein")
            rowsToKeep <- which(!(rownames(countmat) %in% dunno))
            countmat <- countmat[rowsToKeep, , drop = FALSE]
            rawcountmat <- rawcountmat[rownames(countmat), , drop = FALSE]
        }

        #Rename rows to include the description, if not taxonomic data and not MetaCyc Product,
        #whose descriptions are enormous. This matches plot_relabund_heatmap exactly, so that a
        #feature named on a volcano can be found by name on the corresponding heatmap.
        renamed_with_description <- FALSE
        if (!analysis %in% c(taxonomic_spaces, "Product")){
            feattable <- rowData(currobj)
            if (all(c("Accession", "Description") %in% colnames(feattable))){
                newnames <- paste(feattable$Accession, feattable$Description, sep = "-")[match(rownames(countmat), feattable$Accession)]
                #Never let an unmatched accession become NA, which would break every lookup downstream
                newnames[is.na(newnames)] <- rownames(countmat)[is.na(newnames)]
                #Two accessions can share a description, and every lookup from here on is by name,
                #so duplicates would silently resolve to the wrong row.
                if (any(duplicated(newnames))){
                    flog.warn(paste0(sum(duplicated(newnames)), " feature name(s) would be duplicated once the Description is appended to the Accession. They are being made unique with a numeric suffix."))
                    newnames <- make.unique(newnames, sep = "_")
                }
                rownames(countmat) <- newnames
                rownames(rawcountmat) <- newnames
                renamed_with_description <- TRUE
            }
        }

        #Resolve highlight_features within this subset's namespace. Accept either the bare
        #accession or the "Accession-Description" form.
        highlight_here <- character(0)
        if (!is.null(highlight_features)){
            hf <- unique(as.character(highlight_features))
            highlight_here <- hf[hf %in% rownames(countmat)]
            hf_left <- setdiff(hf, highlight_here)
            #Only try the bare-accession form when the rownames really are "Accession-Description".
            #Doing it unconditionally would let a taxonomic name containing a hyphen match the wrong bin.
            if ((length(hf_left) > 0) && renamed_with_description){
                bare <- sub("-.*$", "", rownames(countmat))
                extra <- rownames(countmat)[bare %in% hf_left]
                highlight_here <- unique(c(highlight_here, extra))
            }
            still_missing <- length(setdiff(hf, c(highlight_here, sub("-.*$", "", highlight_here))))
            if (still_missing > 0){
                flog.info(paste0("highlight_features: ", still_missing, " of ", length(hf), " requested feature(s) are not present in ", subset_label, " and will simply not be highlighted."))
            }
        }

        #### One volcano per pair of classes ####
        for (pr in 1:length(pairs_to_plot)){

            if (pages_drawn >= max_total_pages){
                flog.warn(paste0("plot_volcano: ", pages_drawn, " pages have already been committed to the device, which is the cap set by max_total_pages. Stopping here rather than producing an unmanageable document."))
                break
            }

            thispair <- pairs_to_plot[[pr]]

            if (volcano_mode == "correlation"){
                #### Correlation: one comparison, every usable sample, no classes ####
                pairname <- "all_samples"
                pair_label <- paste0(effect_axis_name, " vs ", compareby)
                where_label <- if (is.null(subsetby)) pair_label else paste0(pair_label, " within ", subsetname)

                pairsamples <- samplesToKeep_sp[usable_here]
                pairsamples <- pairsamples[pairsamples %in% colnames(countmat)]
                cbvals <- suppressWarnings(as.numeric(as.character(colData(currobj)[pairsamples, compareby])))
                usable_num <- !is.na(cbvals)
                if (any(!usable_num)){
                    flog.warn(paste0(sum(!usable_num), " sample(s) in ", where_label, " carry a value for \"", compareby, "\" which could not be read as a number and are being dropped from the correlation."))
                }
                pairsamples <- pairsamples[usable_num]
                cbvals <- cbvals[usable_num]
                n_used <- length(pairsamples)

                #Hand-built, because make_classes_df exists to map samples to discrete classes.
                #rownames matter: calculate_matrix_stats reorders the counts matrix by them.
                classesdf <- data.frame(Sample = pairsamples, cl = cbvals, stringsAsFactors = FALSE)
                rownames(classesdf) <- classesdf$Sample

                #Inert in this mode, but referenced by shared code further down
                n1 <- n_used
                n2 <- n_used
                paired_possible <- FALSE
                paired_in_use <- FALSE
                numeric_class_labels <- FALSE
                discretenames <- c(NA_character_, NA_character_)
                display_names <- c(NA_character_, NA_character_)

                if (n_used < min_samples_for_correlation){
                    thisreason <- paste0("Only ", n_used, " sample(s) in ", where_label, " carry a usable numeric value for \"", compareby, "\", and at least ", min_samples_for_correlation, " are needed to correlate against it.")
                    report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = thisreason, advice = "Widen samplesToKeep, or lower min_samples_for_correlation if you accept a correlation on very few points.")
                    nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                    next
                }
                if (length(unique(cbvals)) < 2){
                    thisreason <- paste0("Comparison variable \"", compareby, "\" is constant across every sample in ", where_label, ", so it has no variance and cannot be correlated with anything.")
                    report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = thisreason, advice = "Choose a subsetby variable which does not hold the comparison variable constant.")
                    nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                    next
                }
                flog.info(paste0("Correlating ", nrow(countmat), " feature(s) against \"", compareby, "\" over ", n_used, " sample(s) with a two-sided ", corr_method, " test."))

            } else {

            pairname <- paste(thispair[1], "vs", thispair[2], sep = "_")
            pair_label <- paste0(thispair[1], " vs ", thispair[2])
            where_label <- if (is.null(subsetby)) pair_label else paste0(pair_label, " within ", subsetname)

            pairsamples <- samplesToKeep_sp[usable_here & (cb_here %in% thispair)]
            pairsamples <- pairsamples[pairsamples %in% colnames(countmat)]
            classvec <- as.character(colData(currobj)[pairsamples, compareby])
            n1 <- sum(classvec == thispair[1])
            n2 <- sum(classvec == thispair[2])

            #### Sanity check: enough samples in each class ####
            if (min(n1, n2) < min_samples_per_class){
                thisreason <- paste0("For ", where_label, " there are ", n1, " sample(s) in \"", thispair[1], "\" and ", n2, " in \"", thispair[2], "\", and at least ", min_samples_per_class, " are needed in each class for a Mann-Whitney-Wilcoxon U-test.")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = thisreason, advice = "Widen samplesToKeep, merge classes, or choose a different comparison.")
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }

            #### Map samples to classes ####
            #This is done before the power check, because make_classes_df is what decides
            #whether the test will actually be paired, and the paired and unpaired tests have
            #very different p-value floors.
            pair_pt <- colData(currobj)[pairsamples, , drop = FALSE]
            classesdf <- make_classes_df(curr_pt = pair_pt, compareby = compareby, wilcox_paired_by = wilcox_paired_by)
            paired_possible <- !is.null(wilcox_paired_by)
            paired_in_use <- ("wilcox_pairs" %in% colnames(classesdf))

            #make_classes_df reorders rows when pairing, and can emit NA rows if the pairing
            #is degenerate. Trust it as the record of what will be tested, rather than the
            #vector we handed it, so that the reported group sizes are the ones actually used.
            if (any(is.na(classesdf$Sample)) || any(is.na(classesdf$cl))){
                thisreason <- paste0("Pairing the samples of ", where_label, " by \"", wilcox_paired_by, "\" produced an incomplete sample-to-class mapping, so the comparison cannot be run.")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = c(thisreason, "This happens when the values of the pairing variable do not line up one-for-one across the two classes."), advice = paste0("Check that every value of \"", wilcox_paired_by, "\" occurs exactly twice, once in each class, or set wilcox_paired_by = NULL."))
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }
            pairsamples <- as.character(classesdf$Sample)
            n1 <- sum(as.character(classesdf$cl) == thispair[1])
            n2 <- sum(as.character(classesdf$cl) == thispair[2])

            #Re-apply the size guard to the numbers make_classes_df actually produced, which
            #its pairing logic can change, rather than trusting the pre-pairing count.
            if (min(n1, n2) < min_samples_per_class){
                thisreason <- paste0("Once the samples of ", where_label, " had been mapped to classes", if (paired_in_use) paste0(" and paired by \"", wilcox_paired_by, "\"") else "", ", there were ", n1, " sample(s) in \"", thispair[1], "\" and ", n2, " in \"", thispair[2], "\", and at least ", min_samples_per_class, " are needed in each class.")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = thisreason, advice = "Widen samplesToKeep, merge classes, or drop wilcox_paired_by.")
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }

            #calculate_matrix_stats treats numeric-looking class labels as a continuous variable
            #and then refuses to run a binary test on them. Numerically-coded groups (timepoints
            #0 and 7, dose 0 and 10) are common in JAMS metadata, so prefix the two labels into
            #unambiguously non-numeric strings and map them back for display. A constant prefix
            #cannot change the alphabetical order of two strings, so directionality is preserved.
            numeric_class_labels <- can_be_made_numeric(as.character(classesdf$cl))
            if (numeric_class_labels){
                flog.info(paste0("Class labels \"", thispair[1], "\" and \"", thispair[2], "\" of \"", compareby, "\" are numeric. They are being handled as two discrete categories, ordered as text, for the purposes of this comparison."))
                classesdf$cl <- paste0("JAMSclass_", as.character(classesdf$cl))
            }
            discretenames <- sort(unique(as.character(classesdf$cl)))
            display_names <- if (numeric_class_labels) sub("^JAMSclass_", "", discretenames) else discretenames

            if (length(discretenames) != 2){
                thisreason <- paste0("After assembling the comparison for ", where_label, ", ", length(discretenames), " class(es) were found rather than exactly two. This should not happen and suggests the metadata changed underneath the comparison.")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = thisreason, advice = "Check colData(ExpObj) for stray values in the comparison variable.")
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }

            if (paired_in_use){
                flog.info(paste0("Mann-Whitney-Wilcoxon test between ", display_names[1], " and ", display_names[2], " will be paired by ", wilcox_paired_by))
            } else if (paired_possible){
                flog.warn(paste0("Pairing by \"", wilcox_paired_by, "\" was requested but could not be established for ", where_label, ". Running the UNPAIRED test."))
            }

            n_used <- length(pairsamples)

            } #End of the binary-mode comparison setup

            #### Sanity check: can the test possibly return a significant p-value at all? ####
            #With small groups the Mann-Whitney-Wilcoxon U-test has a hard floor on the p-value
            #it can return, no matter how perfect the separation: 2 / choose(n1 + n2, n1) for
            #the unpaired test, and 2 / 2^n for the signed-rank test on n pairs. Four samples
            #per group are needed before an unpaired two-sided test can reach 0.05 at all, and
            #six pairs before a paired one can. Working this out before computing anything turns
            #a baffling all-grey volcano into an explicit statement.
            #Note that this floor is deliberately NOT inflated for multiple testing. The smallest
            #attainable Benjamini-Hochberg adjusted p-value equals the raw floor, since BH divides
            #by the rank and every feature could in principle tie at that floor.
            #Spearman has the same kind of floor: with n points and no ties the smallest
            #two-sided exact p is 2 / n!, so 5 points are needed before 0.05 is reachable.
            #Pearson has no such floor (a perfect straight line through 3 points gives an
            #arbitrarily small p), so there it is only worth warning that n is small.
            if (volcano_mode == "correlation"){
                if (corr_method == "spearman"){
                    #lfactorial, not factorial: gamma() warns "value out of range" from n = 171 up
                    minp_floor <- min(1, exp(log(2) - lfactorial(n_used)))
                } else {
                    minp_floor <- 0
                }
            } else {
                minp_floor <- min_attainable_wilcox_p(n1 = n1, n2 = n2, paired = paired_in_use)
            }
            floor_name <- "p-value"

            underpowered_msg <- NULL
            if (minp_floor > pval_threshold){
                if (volcano_mode == "correlation"){
                    testname <- "Spearman rank correlation"
                    arrangements <- paste0(n_used, "! = ", format(factorial(n_used), big.mark = ",", scientific = FALSE), " possible orderings")
                    underpowered_msg <- paste0("UNDERPOWERED: with ", n_used, " samples the smallest two-sided exact ", testname, " ", floor_name, " attainable is ", signif(minp_floor, 3), ", which is above pval_threshold = ", pval_threshold, ". No feature can be significant.")
                    reasons <- c(paste0("For ", where_label, " there are ", n_used, " sample(s) carrying a usable value for \"", compareby, "\"."), paste0("Spearman compares ranks, so with ", n_used, " points there are only ", arrangements, ", and the two most extreme of them give the smallest attainable two-sided p of ", signif(minp_floor, 3), "."), paste0("The threshold asked for is ", pval_threshold, ", which is below that floor, so every single feature would necessarily be drawn as non-significant."), "Note that when there are ties (and zero-inflated abundance data are full of them) cor.test falls back on an asymptotic approximation which CAN report smaller p-values than this exact floor. Those are an artefact of the tie correction at this sample size and should not be trusted.")
                } else {
                    testname <- if (paired_in_use) "paired Wilcoxon signed-rank" else "unpaired Mann-Whitney-Wilcoxon"
                    arrangements <- if (paired_in_use) paste0("2^", min(n1, n2), " = ", format(2 ^ min(n1, n2), big.mark = ","), " possible sign patterns") else paste0("choose(", (n1 + n2), ", ", n1, ") = ", format(choose((n1 + n2), n1), big.mark = ","), " possible arrangements")
                    underpowered_msg <- paste0("UNDERPOWERED: with ", n1, " vs ", n2, " samples the smallest two-sided ", testname, " ", floor_name, " attainable is ", signif(minp_floor, 3), ", which is above pval_threshold = ", pval_threshold, ". No feature can be significant.")
                    reasons <- c(paste0("For ", where_label, " there are ", n1, " sample(s) in \"", thispair[1], "\" and ", n2, " in \"", thispair[2], "\"."), paste0("The two-sided ", testname, " test compares ranks, so for these group sizes there are only ", arrangements, ". Even perfect separation therefore cannot give a p-value below ", signif(minp_floor, 3), "."), paste0("The threshold asked for is ", pval_threshold, ", which is below that floor, so every single feature would necessarily be drawn as non-significant."))
                }
                if (isTRUE(adj_pval_for_threshold) || identical(adj_pval_for_threshold, "auto")){
                    reasons <- c(reasons, "FDR adjustment can only make this worse, never better.")
                }
                if (refuse_if_underpowered){
                    #Work out exactly how many samples would be needed, rather than hand-waving
                    need_n <- NA
                    for (kk in 2:60){
                        kk_floor <- if (volcano_mode == "correlation") min(1, exp(log(2) - lfactorial(kk))) else min_attainable_wilcox_p(n1 = kk, n2 = kk, paired = paired_in_use)
                        if (kk_floor <= pval_threshold){
                            need_n <- kk
                            break
                        }
                    }
                    per_what <- if (volcano_mode == "correlation") " samples in total." else " samples per group."
                    if (is.na(need_n)){
                        need_msg <- paste0("No sample size up to 60 would bring the smallest attainable ", floor_name, " down to ", pval_threshold, ".")
                    } else {
                        need_msg <- paste0("The smallest size whose attainable ", floor_name, " reaches ", pval_threshold, " is ", need_n, per_what)
                    }
                    report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label), " - THE COMPARISON CANNOT REACH SIGNIFICANCE"), reasons = reasons, advice = c(need_msg, "Options: widen samplesToKeep; raise pval_threshold; or set refuse_if_underpowered = FALSE to draw the plot anyway and read it purely as a picture of effect sizes."))
                    nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, paste0(where_label, ": underpowered (floor ", signif(minp_floor, 3), " > ", pval_threshold, ")"))
                    next
                } else {
                    flog.warn(underpowered_msg)
                    flog.warn("refuse_if_underpowered = FALSE, so the volcano will be drawn anyway. Read it as a picture of effect sizes, not of significance.")
                }
            }

            smalln_msg <- NULL
            if (volcano_mode == "correlation" && (n_used < 10)){
                smalln_msg <- paste0("NOTE: only ", n_used, " samples. A correlation coefficient estimated from so few points is extremely unstable, and a single outlier can produce an apparently strong, significant correlation.")
                flog.warn(smalln_msg)
            }

            #### Assemble the counts matrix for this pair ####
            countmat_pair <- countmat[ , pairsamples, drop = FALSE]
            raw_pair <- rawcountmat[rownames(countmat_pair), pairsamples, drop = FALSE]
            dropped_absent <- 0
            if (drop_features_absent_in_all_samples && can_judge_absence){
                keeprows <- which(rowSums(raw_pair) > 0)
                dropped_absent <- (nrow(countmat_pair) - length(keeprows))
                if (dropped_absent > 0){
                    flog.info(paste0("Dropping ", dropped_absent, " feature(s) which are absent from every sample of ", where_label, ". This shrinks the multiple-testing denominator, so the adjusted p-values will not match those plot_relabund_heatmap reports for the same comparison."))
                }
                countmat_pair <- countmat_pair[keeprows, , drop = FALSE]
            }

            #A feature with no variance has an undefined correlation with anything: cor.test
            #returns NA, warns about a zero standard deviation, and on some inputs errors
            #outright, which would take the whole comparison down with it. In correlation mode
            #these are therefore removed before testing, not afterwards. This is unconditional,
            #unlike drop_features_absent_in_all_samples, because there is no meaningful answer
            #to keep. Note that it also shrinks the multiple-testing denominator.
            zerovar_dropped <- 0
            if (volcano_mode == "correlation"){
                featsd <- suppressWarnings(apply(countmat_pair, 1, function(x){ stats::sd(as.numeric(x)) }))
                keepvar <- which(is.finite(featsd) & (featsd > 0))
                zerovar_dropped <- (nrow(countmat_pair) - length(keepvar))
                if (zerovar_dropped > 0){
                    flog.info(paste0("Dropping ", zerovar_dropped, " feature(s) which are constant across every sample of ", where_label, ". A constant feature has an undefined correlation, so it cannot be tested; this does reduce the number of features the FDR correction is applied over."))
                }
                countmat_pair <- countmat_pair[keepvar, , drop = FALSE]
            }

            if (nrow(countmat_pair) < min_features_for_volcano){
                thisreason <- paste0("Only ", nrow(countmat_pair), " feature(s) remain for ", where_label, ", which is fewer than the ", min_features_for_volcano, " needed to draw a volcano.")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = c(thisreason, paste0("Filtration applied: ", paste0(presetlist$filtermsg, collapse = " | "))), advice = "Relax applyfilters and featuresToKeep, or use plot_relabund_features for a handful of features.")
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }
            if (nrow(countmat_pair) < 10){
                flog.warn(paste0("Only ", nrow(countmat_pair), " feature(s) are being tested for ", where_label, ". A volcano of so few points conveys little, and any FDR adjustment over so small a denominator is barely a correction at all."))
            }

            #### Compute the statistics ####

            stattype_here <- if (volcano_mode == "correlation") corr_method else "binary"

            #Rank tests on zero-inflated abundance data emit one "cannot compute exact p-value
            #with ties" warning PER FEATURE. Left alone that is hundreds of warnings for a single
            #plot. Count them and say so once instead, and let anything unexpected through.
            tie_warnings <- 0
            matstats <- tryCatch({
                withCallingHandlers({
                    calculate_matrix_stats(countmatrix = countmat_pair, uselog = FALSE, statsonlog = statsonlog, stattype = stattype_here, classesdf = classesdf, invertbinaryorder = invertbinaryorder, numthreads = numthreads, fun_for_l2fc = fun_for_l2fc)
                }, warning = function(w){
                    if (grepl("exact p-value|ties|standard deviation is zero", conditionMessage(w))){
                        tie_warnings <<- (tie_warnings + 1)
                        invokeRestart("muffleWarning")
                    }
                })
            }, error = function(e){
                return(simpleError(conditionMessage(e)))
            })

            ties_msg <- NULL
            if (tie_warnings > 0){
                ties_msg <- paste0("NOTE: ", tie_warnings, " feature(s) contained tied values, so their p-values come from the asymptotic approximation rather than from an exact test. This is normal for zero-inflated abundance data.")
                flog.info(ties_msg)
            }

            if (inherits(matstats, "error") || !is.data.frame(matstats) || (nrow(matstats) < 1)){
                errmsg <- if (inherits(matstats, "error")) conditionMessage(matstats) else "the statistics routine returned nothing usable"
                thisreason <- paste0("The statistics could not be computed for ", where_label, ": ", errmsg)
                scope_msg <- if (volcano_mode == "correlation") paste0("The correlation used ", n_used, " sample(s), over ", nrow(countmat_pair), " features.") else paste0("Groups were ", n1, " x \"", thispair[1], "\" and ", n2, " x \"", thispair[2], "\", over ", nrow(countmat_pair), " features.")
                scope_advice <- if (volcano_mode == "correlation") "This usually means degenerate data. Constant features are already removed in correlation mode, so suspect the comparison variable itself: check that it really varies across these samples." else "This usually means degenerate data, such as a feature which is constant or missing across every sample. Try relaxing the filters, or set drop_features_absent_in_all_samples = TRUE."
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = c(thisreason, scope_msg), advice = scope_advice)
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }

            #Restore the user-facing class labels on the directionality column
            if (numeric_class_labels && ("Class_where_increased" %in% colnames(matstats))){
                matstats$Class_where_increased <- sub("^JAMSclass_", "", matstats$Class_where_increased)
            }

            #### A centred log-ratio is already a log-ratio, so a *ratio* of CLR values is not ####
            #### a fold change. Replace the effect size with the difference of group means of  ####
            #### CLR, converted from natural-log to log2 units.                                ####
            if (volcano_mode == "binary" && normalization == "clr" && assay_for_matrix == "BaseCounts"){
                g1 <- classesdf$Sample[which(as.character(classesdf$cl) == discretenames[1])]
                g2 <- classesdf$Sample[which(as.character(classesdf$cl) == discretenames[2])]
                clrdiff <- (rowMeans(countmat_pair[ , g1, drop = FALSE]) - rowMeans(countmat_pair[ , g2, drop = FALSE])) / log(2)
                if (invertbinaryorder == TRUE){
                    clrdiff <- (clrdiff * -1)
                }
                matstats$l2fc <- as.numeric(clrdiff[rownames(matstats)])
                matstats$absl2fc <- abs(matstats$l2fc)
                matstats$Class_where_increased <- NA_character_
                matstats$Class_where_increased[which(matstats$l2fc == 0)] <- "neither"
                matstats$Class_where_increased[which(matstats$l2fc > 0)] <- display_names[1:2][c(!invertbinaryorder, invertbinaryorder)]
                matstats$Class_where_increased[which(matstats$l2fc < 0)] <- display_names[1:2][c(invertbinaryorder, !invertbinaryorder)]
                effect_size_name <- "CLR difference (log2 units)"
                effect_axis_name <- "CLR difference (log2 units)"
                flog.info("Normalization is CLR, so the effect size plotted is the difference between the group means of the centred log-ratio, expressed in log2 units, rather than a ratio-based fold change.")
            }

            #From here on everything is shared. In correlation mode the coefficient takes the
            #place the fold change occupies in binary mode, so the plotting code below needs
            #no further knowledge of which test produced it.
            if (volcano_mode == "correlation"){
                if (!("correl" %in% colnames(matstats))){
                    thisreason <- paste0("The ", corr_method, " test returned no correlation coefficient for ", where_label, ".")
                    report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = thisreason, advice = "This should not happen; check that compareby really is numeric for these samples.")
                    nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                    next
                }
                matstats$l2fc <- as.numeric(matstats$correl)
                matstats$absl2fc <- abs(matstats$l2fc)
            }

            #### Bank the raw statistics for every feature, exactly as plot_relabund_heatmap does ####
            matstatsbank <- as.data.frame(matstats)
            if ("GenomeCompleteness" %in% names(assays(currobj))){
                matstatsbank$Taxa <- rownames(matstatsbank)
                gc_rows <- rownames(matstatsbank)[rownames(matstatsbank) %in% rownames(currobj)]
                if (length(gc_rows) > 0){
                    if (!use_checkM2_style_GenomeCompleteness && ("GenomeContamination" %in% names(assays(currobj)))){
                        genomecompletenessmat <- as.matrix(assays(currobj)$GenomeCompleteness)[gc_rows, pairsamples, drop = FALSE] + as.matrix(assays(currobj)$GenomeContamination)[gc_rows, pairsamples, drop = FALSE]
                    } else {
                        genomecompletenessmat <- as.matrix(assays(currobj)$GenomeCompleteness)[gc_rows, pairsamples, drop = FALSE]
                    }
                    genomecompletenessdf <- data.frame(Taxa = gc_rows, MedianGenomeComp = rowMedians(genomecompletenessmat), SDGenomeComp = jams_rowSds(genomecompletenessmat), stringsAsFactors = FALSE)
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdf, by = "Taxa")
                    rownames(matstatsbank) <- matstatsbank$Taxa
                }
                matstatsbank$Taxa <- NULL
            }

            #Name and bank immediately. The verdict columns are added, and this entry
            #overwritten, further down if the comparison turns out to be plottable.
            stattitle <- paste(analysisname, paste(if (volcano_mode == "correlation") toupper(substr(corr_method, 1, 4)) else "MWW", compareby, sep = "_"), subsetname, pairname, sep = "_")
            if (stattitle %in% names(svec)){
                stattitle <- paste(stattitle, length(svec), sep = "_")
            }
            svec[[stattitle]] <- matstatsbank

            #### Decide which p-value is being plotted and thresholded ####
            adj_here <- adj_pval_for_threshold
            if (identical(adj_here, "auto")){
                propsigadj <- length(which(matstats$padj_fdr < 0.05)) / length(matstats$padj_fdr)
                propsignonadj <- length(which(matstats$pval < 0.05)) / length(matstats$pval)
                fracsigadj <- propsigadj / propsignonadj
                adj_here <- ((!is.na(fracsigadj)) && (fracsigadj > 0.2))
                flog.info(paste0("P-value adjustment set to auto. Upon FDR adjustment, the proportion still below 0.05 is ", round((fracsigadj * 100), 1), "% of the unadjusted proportion, so p-value adjustment is set to ", as.character(adj_here), "."))
            }
            sigmeas <- if (isTRUE(adj_here)) "padj_fdr" else "pval"
            sigmeas_pretty <- if (isTRUE(adj_here)) "FDR-adjusted p" else "p"
            othermeas <- if (isTRUE(adj_here)) "pval" else "padj_fdr"
            othermeas_pretty <- if (isTRUE(adj_here)) "raw p" else "FDR-adjusted p"

            #### Assemble the plotting data frame ####
            vdf <- data.frame(Feature = rownames(matstats), l2fc = as.numeric(matstats$l2fc), pval = as.numeric(matstats$pval), padj_fdr = as.numeric(matstats$padj_fdr), stringsAsFactors = FALSE)
            vdf$pplot <- vdf[[sigmeas]]

            n_tested <- nrow(vdf)
            badrows <- which(is.na(vdf$l2fc) | is.na(vdf$pplot) | is.nan(vdf$l2fc) | is.nan(vdf$pplot))
            na_msg <- NULL
            if (length(badrows) > 0){
                na_msg <- paste0("NOTE: ", length(badrows), " of ", n_tested, " feature(s) yielded an unusable statistic (NA or NaN) and are not shown.")
                flog.warn(na_msg)
                vdf <- vdf[-badrows, , drop = FALSE]
            }
            if (nrow(vdf) < min_features_for_volcano){
                thisreason <- paste0("After discarding features with unusable statistics, only ", nrow(vdf), " feature(s) remain for ", where_label, ".")
                report_refusal(headline = paste0("NO VOLCANO PLOT FOR ", toupper(where_label)), reasons = c(thisreason, na_msg), advice = "This normally means the data in this comparison are degenerate. Check the counts matrix for constant or empty features.")
                nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, thisreason)
                next
            }

            vdf$neglog10p <- -log10(vdf$pplot)

            #### Clamp infinities so that they are shown at the edge rather than silently dropped ####
            finite_l2fc <- vdf$l2fc[is.finite(vdf$l2fc)]
            n_inf_x <- sum(!is.finite(vdf$l2fc))
            if (volcano_mode == "correlation"){
                #A correlation coefficient lives in [-1, 1] by construction, so fixing the axis
                #there makes every correlation volcano directly comparable with every other,
                #and no point can ever need capping. cap_l2fc_at has nothing to do here.
                xcap <- 1
                if (!is.null(cap_l2fc_at)){
                    flog.info("cap_l2fc_at has no effect in correlation mode: the X axis is fixed to the full [-1, 1] range of the coefficient.")
                }
            } else if (!is.null(cap_l2fc_at)){
                #Never let an explicit cap hide the threshold lines themselves
                xcap <- max(abs(as.numeric(cap_l2fc_at)), (l2fc_threshold * 1.05), 0.1)
            } else if (length(finite_l2fc) > 0){
                xcap <- max(c(abs(finite_l2fc), (l2fc_threshold * 1.2), 0.5)) * 1.05
            } else {
                xcap <- max((l2fc_threshold * 1.5), 1)
            }
            vdf$l2fc_plot <- pmin(pmax(vdf$l2fc, (xcap * -1)), xcap)
            vdf$l2fc_plot[is.na(vdf$l2fc_plot)] <- 0

            finite_y <- vdf$neglog10p[is.finite(vdf$neglog10p)]
            n_inf_y <- sum(!is.finite(vdf$neglog10p))
            if (!is.null(cap_neglog10p_at)){
                ycap <- max(abs(as.numeric(cap_neglog10p_at)), (-log10(pval_threshold) * 1.05))
            } else if (length(finite_y) > 0){
                ycap <- max(c(finite_y, -log10(pval_threshold), 1)) * 1.05
            } else {
                ycap <- max(-log10(pval_threshold) * 2, 2)
            }
            vdf$neglog10p_plot <- pmin(vdf$neglog10p, ycap)
            vdf$Capped <- ((abs(vdf$l2fc) > xcap) | !is.finite(vdf$l2fc) | (vdf$neglog10p > ycap) | !is.finite(vdf$neglog10p))
            vdf$Capped[is.na(vdf$Capped)] <- FALSE
            cap_msg <- NULL
            if (any(vdf$Capped)){
                capbits <- character(0)
                if (n_inf_x > 0){
                    capbits <- c(capbits, paste0(n_inf_x, " with an infinite fold change (feature wholly absent from one group; use fun_for_l2fc = \"geom_mean\" to avoid this)"))
                }
                if (n_inf_y > 0){
                    capbits <- c(capbits, paste0(n_inf_y, " with a p-value that underflowed to zero"))
                }
                cap_msg <- paste0("NOTE: ", sum(vdf$Capped), " feature(s) are drawn as triangles at the plot edge because they lie beyond the axis limits", if (length(capbits) > 0) paste0(", including ", paste0(capbits, collapse = " and ")) else "", ".")
            }

            #### Direction and significance ####
            if (volcano_mode == "correlation"){
                up_class <- paste0("increasing ", compareby)
                down_class <- paste0("decreasing ", compareby)
                lvl_up <- paste0("Positively correlated with ", compareby)
                lvl_dn <- paste0("Negatively correlated with ", compareby)
            } else {
                #Matches plot_relabund_heatmap: positive l2fc means increased in binary_directionality.
                up_class <- display_names[1:2][c(!invertbinaryorder, invertbinaryorder)]
                down_class <- display_names[1:2][c(invertbinaryorder, !invertbinaryorder)]
                lvl_up <- paste0("Up in ", up_class)
                lvl_dn <- paste0("Up in ", down_class)
            }
            eff_thresh <- max(effect_threshold, .Machine$double.eps)

            vdf$Significant <- (vdf$pplot <= pval_threshold)
            lvl_ns <- "Not significant"
            statusvec <- rep(lvl_ns, nrow(vdf))
            statusvec[which(vdf$Significant & (vdf$l2fc >= eff_thresh))] <- lvl_up
            statusvec[which(vdf$Significant & (vdf$l2fc <= (eff_thresh * -1)))] <- lvl_dn
            vdf$Status <- factor(statusvec, levels = c(lvl_up, lvl_ns, lvl_dn))

            n_up <- sum(statusvec == lvl_up)
            n_dn <- sum(statusvec == lvl_dn)
            n_ns <- sum(statusvec == lvl_ns)
            n_other <- sum(vdf[[othermeas]] <= pval_threshold, na.rm = TRUE)

            #Record the verdict on the banked stats, so that returnstats gives back exactly what was drawn.
            #Volcano_status carries the human-readable label used in the legend; Volcano_direction
            #carries the plain Up / Down / NS token, which is what you want for subsetting.
            #Features whose statistic came back NA were tested, and did count towards the
            #multiple-testing correction, but could not be placed on the plot. They are marked
            #"Not plotted" rather than being conflated with features that were never tested.
            matstatsbank$neglog10_plotted_pval <- -log10(matstatsbank[[sigmeas]])
            matstatsbank$Volcano_status <- "Not plotted"
            matstatsbank$Volcano_direction <- "Not plotted"
            matstatsbank[vdf$Feature, "Volcano_status"] <- as.character(vdf$Status)
            dirtoken <- rep("NS", nrow(vdf))
            dirtoken[statusvec == lvl_up] <- if (volcano_mode == "correlation") "Positive" else "Up"
            dirtoken[statusvec == lvl_dn] <- if (volcano_mode == "correlation") "Negative" else "Down"
            matstatsbank[vdf$Feature, "Volcano_direction"] <- dirtoken
            matstatsbank$Volcano_up_class <- up_class
            matstatsbank$Volcano_down_class <- down_class
            svec[[stattitle]] <- matstatsbank

            if ((n_up + n_dn) == 0){
                flog.warn(paste0("No feature reached ", sigmeas_pretty, " <= ", pval_threshold, " with an absolute ", effect_axis_name, " >= ", round(effect_threshold, 2), " for ", where_label, "."))
                if (!plot_when_nothing_significant){
                    report_refusal(headline = paste0("NO SIGNIFICANT FEATURES FOR ", toupper(where_label)), reasons = c(paste0("All ", nrow(vdf), " feature(s) tested fell short of ", sigmeas_pretty, " <= ", pval_threshold, " combined with an absolute ", effect_axis_name, " >= ", round(effect_threshold, 2), "."), if (volcano_mode == "correlation") paste0("The correlation used ", n_used, " samples.") else paste0("Groups were ", n1, " x \"", thispair[1], "\" and ", n2, " x \"", thispair[2], "\"."), paste0("Of those features, ", sum(vdf$Significant), " passed the p-value threshold on its own, and ", sum(abs(vdf$l2fc) >= eff_thresh, na.rm = TRUE), " passed the effect-size threshold on its own."), paste0("Filtration applied: ", paste0(presetlist$filtermsg, collapse = " | "))), advice = "Set plot_when_nothing_significant = TRUE (the default) to draw the volcano anyway; an empty volcano is itself a result. Otherwise relax pval_threshold, or the effect-size threshold (l2fc_threshold in binary mode, corr_threshold in correlation mode).")
                    nothing_was_drawn_reasons <- c(nothing_was_drawn_reasons, paste0(where_label, ": no feature passed both thresholds"))
                    next
                }
            }

            #### Colours: reuse the colours these two classes already have elsewhere in JAMS ####
            upcol <- up_colour
            dncol <- down_colour
            colour_source <- "fixed JAMS defaults"
            if ((is.null(upcol) || is.null(dncol)) && use_class_colours && (volcano_mode == "binary")){
                lookup_vars <- unique(c(compareby, compareby_original))
                c_up <- resolve_class_colour(classname = up_class, variable = lookup_vars, SEobj = currobj)
                c_dn <- resolve_class_colour(classname = down_class, variable = lookup_vars, SEobj = currobj)
                if (!is.na(c_up) && !is.na(c_dn) && !identical(tolower(c_up), tolower(c_dn))){
                    if (is.null(upcol)){
                        upcol <- c_up
                    }
                    if (is.null(dncol)){
                        dncol <- c_dn
                    }
                    colour_source <- "the colours already assigned to these classes in this project"
                } else {
                    if (is.na(c_up) || is.na(c_dn)){
                        flog.info(paste0("No colour was found for \"", if (is.na(c_up)) up_class else down_class, "\" in the colour dictionary or in the object's ctable, so the fixed red/blue pair will be used for this volcano."))
                    } else {
                        flog.info(paste0("Classes \"", up_class, "\" and \"", down_class, "\" resolve to the same colour (", c_up, "), which would make the two sides of the volcano indistinguishable. Using the fixed red/blue pair instead."))
                    }
                }
            }
            if (is.null(upcol)){
                upcol <- fallback_up_colour
            }
            if (is.null(dncol)){
                dncol <- fallback_down_colour
            }

            statuscols <- c(upcol, nonsignificant_colour, dncol)
            names(statuscols) <- c(lvl_up, lvl_ns, lvl_dn)
            statuslabels <- c(paste0(lvl_up, " (n = ", n_up, ")"), paste0(lvl_ns, " (n = ", n_ns, ")"), paste0(lvl_dn, " (n = ", n_dn, ")"))
            names(statuslabels) <- c(lvl_up, lvl_ns, lvl_dn)

            #### Which features get named ####
            labdf <- NULL
            lab_msg <- NULL
            if (label_features){
                eff_lab_thresh <- max(label_effect_threshold, .Machine$double.eps)
                hl_idx <- which(vdf$Feature %in% highlight_here)
                cand <- which((vdf$pplot <= label_pval_threshold) & (abs(vdf$l2fc) >= eff_lab_thresh))
                #Highlighted features are labelled unconditionally, so they must not eat into the cap
                cand <- setdiff(cand, hl_idx)
                #Rank by how far out into the corners of the volcano a feature sits
                corner_score <- (vdf$neglog10p_plot * pmin(abs(vdf$l2fc_plot), xcap))
                corner_score[is.na(corner_score)] <- 0
                n_candidates <- length(cand)
                if (n_candidates > 0){
                    cand <- cand[order(corner_score[cand], decreasing = TRUE)]
                }
                if (n_candidates > max_features_to_label){
                    cand <- utils::head(cand, max_features_to_label)
                    lab_msg <- paste0("Naming the ", max_features_to_label, " of ", n_candidates, " qualifying feature(s) which sit furthest into the corners of the volcano; raise max_features_to_label to name more.")
                }
                labidx <- unique(c(cand, hl_idx))
                if (length(labidx) > 0){
                    labdf <- vdf[labidx, , drop = FALSE]
                    labtxt <- labdf$Feature
                    if (no_underscores && (analysis %in% taxonomic_spaces)){
                        labtxt <- gsub("_", " ", labtxt)
                    }
                    labdf$Label <- wrap_feature_label(x = labtxt, width = label_wrap_width, maxchar = max_label_nchar)
                    if (is.null(label_fontsize)){
                        nlab <- nrow(labdf)
                        lab_fs <- max(1.5, min(3.5, (11 / sqrt(max(nlab, 1)))))
                        maxlines <- max(vapply(strsplit(labdf$Label, "\n", fixed = TRUE), length, integer(1)))
                        if (maxlines > 2){
                            lab_fs <- (lab_fs * 0.88)
                        }
                        lab_fs <- round(lab_fs, 2)
                    } else {
                        lab_fs <- as.numeric(label_fontsize)
                    }
                }
            }

            #### Point size ####
            if (is.null(point_size)){
                pt_sz <- round(max(0.6, min(2.2, (30 / sqrt(max(nrow(vdf), 1))))), 2)
            } else {
                pt_sz <- as.numeric(point_size)
            }

            #### Titles, built to match plot_relabund_heatmap ####
            maintit <- paste(c(hmtypemsg, analysisname, pair_label, subsettit), collapse = " | ")
            if (!is.null(addtit)){
                maintit <- paste(addtit, maintit, sep = "\n")
            }

            if (volcano_mode == "correlation"){
                stattit <- paste0("Two-sided ", if (corr_method == "spearman") "Spearman rank" else "Pearson", " correlation of each feature against ", compareby, " | n = ", n_used, " samples")
            } else {
                stattit <- paste0("Mann-Whitney-Wilcoxon U-test | ", compareby, ": ", thispair[1], " (n = ", n1, ") vs ", thispair[2], " (n = ", n2, ")")
                if (paired_in_use){
                    stattit <- paste0(stattit, " | paired by ", wilcox_paired_by)
                }
            }
            nsampmsg <- paste0("Number of samples in comparison = ", length(pairsamples))
            #n_tested, not nrow(vdf), is the denominator p.adjust actually used: features whose
            #statistic came back NA were counted by the correction but cannot be plotted.
            ffeatmsg <- paste0("Number of features assessed = ", n_tested, if (nrow(vdf) != n_tested) paste0(" (", nrow(vdf), " plottable)") else "")
            sigrule <- paste0("Significance: ", sigmeas_pretty, " <= ", pval_threshold, " AND absolute ", effect_axis_name, " >= ", round(effect_threshold, 2))
            if (volcano_mode == "correlation"){
                binary_directionality_msg <- paste0("A positive ", effect_axis_name, " means the feature becomes more abundant as ", compareby, " increases")
            } else {
                binary_directionality_msg <- paste0("Positive ", effect_size_name, " means increased in ", up_class)
            }
            countmsg <- paste0(lvl_up, ": ", n_up, " | ", lvl_dn, ": ", n_dn, " | not significant: ", n_ns, " | features passing ", othermeas_pretty, " <= ", pval_threshold, ": ", n_other)

            submsgs <- c(stattit)
            if (dichotomized_compareby && !is.null(dichot_msg)){
                submsgs <- c(submsgs, dichot_msg)
            }
            if (!is.null(presetlist$filtermsg)){
                submsgs <- c(submsgs, presetlist$filtermsg)
            }
            submsgs <- c(submsgs, paste(nsampmsg, ffeatmsg, sep = " | "), sigrule, binary_directionality_msg, countmsg)
            if (!is.null(featuresToKeep)){
                submsgs <- c(submsgs, paste0("NOTE: featuresToKeep restricted this comparison to ", n_tested, " feature(s), so every FDR-adjusted p-value is corrected over that reduced set."))
            }
            if (zerovar_dropped > 0){
                submsgs <- c(submsgs, paste0("NOTE: ", zerovar_dropped, " feature(s) were constant across all samples and were dropped before testing, because a constant feature has no correlation. This reduces the number of features the FDR correction is applied over."))
            }
            if (!is.null(smalln_msg)){
                submsgs <- c(submsgs, smalln_msg)
            }
            if (!is.null(ties_msg)){
                submsgs <- c(submsgs, ties_msg)
            }
            if (dropped_absent > 0){
                submsgs <- c(submsgs, paste0("NOTE: ", dropped_absent, " feature(s) absent from every sample were dropped before testing (drop_features_absent_in_all_samples = TRUE)."))
            }
            if (!is.null(na_msg)){
                submsgs <- c(submsgs, na_msg)
            }
            if (!is.null(cap_msg)){
                submsgs <- c(submsgs, cap_msg)
            }
            if ((n_up + n_dn) == 0){
                closest_p <- suppressWarnings(min(vdf$pplot, na.rm = TRUE))
                finite_eff <- abs(vdf$l2fc[is.finite(vdf$l2fc)])
                biggest_eff <- if (length(finite_eff) > 0) signif(max(finite_eff), 3) else "Inf"
                submsgs <- c(submsgs, paste0("NOTE: NO feature passed both thresholds. The closest was ", sigmeas_pretty, " = ", signif(closest_p, 3), " and the largest finite absolute effect was ", biggest_eff, "."))
            }
            if (!is.null(underpowered_msg)){
                submsgs <- c(submsgs, underpowered_msg)
            }
            if (!is.null(lab_msg)){
                submsgs <- c(submsgs, lab_msg)
            }
            if (colour_source != "fixed JAMS defaults"){
                submsgs <- c(submsgs, paste0("Point colours are ", colour_source, "."))
            }
            subtit <- paste0(submsgs, collapse = "\n")

            if (volcano_mode == "correlation"){
                xlab_txt <- paste0(effect_axis_name, " against ", compareby, "\n<--  falls as ", compareby, " rises          rises as ", compareby, " rises  -->")
            } else {
                xlab_txt <- effect_axis_name
                if (identical(effect_axis_name, "log2 fold change")){
                    xlab_txt <- paste0("log2 fold change (", fun_for_l2fc, ")")
                }
                xlab_txt <- paste0(xlab_txt, "\n<--  higher in ", down_class, "          higher in ", up_class, "  -->")
            }
            ylab_txt <- paste0("-log10( ", if (isTRUE(adj_here)) "FDR-adjusted p-value" else "p-value", " )")

            #### Draw ####
            featfontface <- if (analysis %in% taxonomic_spaces) "italic" else "plain"

            p <- ggplot(vdf, aes(x = l2fc_plot, y = neglog10p_plot, colour = Status))

            if (show_threshold_lines){
                p <- p + geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", colour = "grey45", linewidth = 0.3)
                if (effect_threshold > 0){
                    p <- p + geom_vline(xintercept = c((effect_threshold * -1), effect_threshold), linetype = "dashed", colour = "grey45", linewidth = 0.3)
                }
                p <- p + geom_vline(xintercept = 0, linetype = "solid", colour = "grey85", linewidth = 0.25)
            }

            notcapped <- vdf[!vdf$Capped, , drop = FALSE]
            iscapped <- vdf[vdf$Capped, , drop = FALSE]
            if (nrow(notcapped) > 0){
                p <- p + geom_point(data = notcapped, shape = 16, size = pt_sz, alpha = point_alpha)
            }
            if (nrow(iscapped) > 0){
                p <- p + geom_point(data = iscapped, shape = 17, size = (pt_sz * 1.35), alpha = min(1, (point_alpha + 0.2)))
            }

            if (length(highlight_here) > 0){
                hldf <- vdf[vdf$Feature %in% highlight_here, , drop = FALSE]
                if (nrow(hldf) > 0){
                    p <- p + geom_point(data = hldf, shape = 21, size = (pt_sz * 2.2), stroke = 0.45, colour = "black", fill = NA, inherit.aes = FALSE, mapping = aes(x = l2fc_plot, y = neglog10p_plot))
                }
            }

            p <- p + scale_colour_manual(values = statuscols, labels = statuslabels, breaks = c(lvl_up, lvl_ns, lvl_dn), drop = FALSE, name = NULL)
            p <- p + coord_cartesian(xlim = c((xcap * -1.03), (xcap * 1.03)), ylim = c(0, (ycap * 1.02)))

            if (!is.null(labdf) && (nrow(labdf) > 0)){
                p <- p + ggrepel::geom_text_repel(data = labdf, mapping = aes(x = l2fc_plot, y = neglog10p_plot, label = Label), size = lab_fs, lineheight = 0.82, fontface = featfontface, colour = "grey15", segment.color = "grey55", segment.size = 0.2, min.segment.length = 0, box.padding = 0.32, point.padding = 0.18, max.overlaps = Inf, seed = 2138, show.legend = FALSE, inherit.aes = FALSE)
            }

            p <- p + labs(title = maintit, subtitle = subtit, x = xlab_txt, y = ylab_txt)
            p <- p + theme_bw()
            p <- p + theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 6, colour = "grey25"), axis.title = element_text(size = 8), axis.text = element_text(size = 7), legend.text = element_text(size = 7), legend.position = "bottom", legend.key.size = unit(4, "mm"), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"))
            p <- p + guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

            print(p)
            pages_drawn <- (pages_drawn + 1)
            volcanoes_drawn <- (volcanoes_drawn + 1)
            if (volcano_mode == "correlation"){
                flog.info(paste0("Volcano plot drawn for ", where_label, ": ", n_up, " positively correlated with ", compareby, ", ", n_dn, " negatively correlated, ", n_ns, " not significant."))
            } else {
                flog.info(paste0("Volcano plot drawn for ", where_label, ": ", n_up, " up in ", up_class, ", ", n_dn, " up in ", down_class, ", ", n_ns, " not significant."))
            }

        } #End loop over pairs of classes

        if (pages_drawn >= max_total_pages){
            break
        }

    } #End subsetby loop

    if (volcanoes_drawn < 1){
        report_refusal(headline = "NO VOLCANO PLOT COULD BE DRAWN FOR ANY COMPARISON", reasons = c(paste0("Comparison variable: \"", compareby_original, "\"", if (dichotomized_compareby) " (dichotomized)" else "", "."), if (length(nothing_was_drawn_reasons) > 0) "The reasons recorded were:" else NULL, utils::head(unique(nothing_was_drawn_reasons), 12)), advice = "Each failed comparison is also reported on its own page above, and in the JAMS log.")
    } else {
        flog.info(paste0("Volcano plots generated: ", volcanoes_drawn, "."))
    }

    if (returnstats == TRUE){
        return(svec)
    }

    return(invisible(NULL))
}
