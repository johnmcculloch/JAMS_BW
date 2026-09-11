#' titrate_feature_counts(ExpObj = NULL, glomby = NULL, only_allow_CSBs = FALSE, samplesToKeep = NULL, PPM_thresholds = NULL, completeness_thresholds = NULL, n_ppm = 10, n_comp = 10, prevalence_pct = 5, PPM_normalize_to_bases_sequenced = FALSE, ignoreunclassified = TRUE, class_to_ignore = "N_A", return_plots = TRUE, verbose = TRUE)
#'
#' A FAST, correlation-free companion to titrate_correlation_stability. For a taxonomic
#' SummarizedExperiment object, it counts how many features survive each combination of a
#' relative-abundance (PPM) threshold and a genome-completeness threshold, on a grid, at a fixed
#' prevalence. No correlations, no permutations: this is just repeated feature counting, so it runs
#' in a fraction of a second on a laptop even for thousands of features.
#'
#' The goal is to let a user see, for THEIR dataset, where the surviving-feature count stops moving
#' as thresholds tighten. Because raw counts decline smoothly (there is rarely a sharp elbow), the
#' function also reports the LOCAL RATE of decline (the marginal feature loss per grid step), which
#' is where turning points such as the recurring ~10% genome-completeness floor actually show up.
#'
#' Crucially, it ALSO reports how much SEQUENCING DEPTH the surviving features represent, per cell.
#' Feature-count loss and depth loss can diverge enormously: axing 90% of features might cost only
#' 5% of depth (you removed low-abundance noise) or 50% of depth (you removed real community),
#' depending on each sample's evenness. Depth is summarised as the MEDIAN retained fraction across
#' samples with the WORST-hit sample (min) alongside, and a divergence panel shows depth% minus
#' feature% directly, so a threshold that looks aggressive by count but cheap by depth is visible.
#'
#' Returned plots: feature_count, decline_rate, depth_retained, count_vs_depth_divergence. The
#' returned grid data frame carries n_features, pct_retained, pct_depth_median, pct_depth_min.
#'
#' @param ExpObj A JAMS-style taxonomic SummarizedExperiment object (e.g. ConsolidatedGenomeBin).
#' @param glomby Taxonomic level to agglomerate to (e.g. "Species"). NULL for the native level.
#' @param only_allow_CSBs Restrict to Consolidated Species Bins (MetaBAT2 bins) only.
#' @param samplesToKeep Optional vector of samples to restrict to.
#' @param PPM_thresholds Optional explicit vector of PPM thresholds. If NULL, a sensible ladder of
#'   length n_ppm is built automatically from the data (see n_ppm).
#' @param completeness_thresholds Optional explicit vector of completeness thresholds. If NULL, an
#'   evenly spaced ladder of length n_comp from 0 to a data-driven maximum is built.
#' @param n_ppm Number of PPM thresholds to auto-generate when PPM_thresholds is NULL. The ladder
#'   starts at 0 and follows a roughly log spacing (0, then quantile-like steps up to a high value)
#'   so the informative low-PPM region is well sampled. Default 10.
#' @param n_comp Number of completeness thresholds to auto-generate when completeness_thresholds is
#'   NULL. Evenly spaced from 0. Default 10.
#' @param prevalence_pct Fixed prevalence (the "in at least X percent of samples" companion applied
#'   to BOTH the PPM and completeness rules). Held constant across the grid. Default 5. Set the
#'   partner argument prevalence_min_samples to override this with an absolute sample count.
#' @param prevalence_min_samples Optional absolute sample-count prevalence. If given (non-NULL), a
#'   feature must meet the threshold in at least this many samples, OVERRIDING prevalence_pct. This
#'   exists because "5% of samples" means very different things at n=10 vs n=300; a count is often
#'   the more honest criterion, and setting it to 1 recovers the "present in at least one sample"
#'   behaviour (useful for detecting a single-sample bug while still censoring empty/low-quality
#'   bins via the completeness axis, e.g. completeness >= 20 in >= 1 sample).
#' @param PPM_normalize_to_bases_sequenced Passed to filter_experiment for PPM computation.
#' @param ignoreunclassified Drop the unclassified/none feature before counting.
#' @param class_to_ignore Metadata classes to drop (signature parity; selection is by samplesToKeep).
#' @param return_plots If TRUE, returns list(grid = <df>, plots = <list>). If FALSE, returns the
#'   data frame invisibly.
#' @param plot_on_the_fly If TRUE, prints plots and returns NULL.
#' @param verbose Emit brief progress via flog.info.
#' @export

titrate_feature_counts <- function(ExpObj = NULL, glomby = NULL, only_allow_CSBs = FALSE, samplesToKeep = NULL, PPM_thresholds = NULL, completeness_thresholds = NULL, n_ppm = 10, n_comp = 10, prevalence_pct = 5, prevalence_min_samples = NULL, PPM_normalize_to_bases_sequenced = FALSE, ignoreunclassified = TRUE, class_to_ignore = "N_A", return_plots = TRUE, plot_on_the_fly = FALSE, verbose = TRUE){

    require(ggplot2)

    #Vet + agglomerate + CSB-restrict up front, once.
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = NULL, only_allow_CSBs = only_allow_CSBs, glomby = glomby, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    analysisname <- if (!is.null(glomby)) glomby else analysis
    facet_label <- paste0(analysisname, if (only_allow_CSBs) " | CSB-only" else " | all-contigs")

    #Build the PPM matrix (for the abundance rule and % feature retention) AND the raw BaseCounts
    #matrix (sequencing depth, for the % depth retention). Both come from the same vetted object.
    baseobj <- filter_experiment(SEobj = obj, featcutoff = c(0, 0), samplesToKeep = NULL, featuresToKeep = NULL, only_allow_CSBs = FALSE, normalization = "relabund", PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = c(0, 0), give_info = FALSE)
    PPMmat <- as.matrix(assays(baseobj)$PPM)
    BaseMat <- as.matrix(assays(baseobj)$BaseCounts)

    have_GC <- "GenomeCompleteness" %in% names(assays(baseobj))
    GCmat <- if (have_GC) as.matrix(assays(baseobj)$GenomeCompleteness) else NULL

    #Drop unclassified/none rows from ALL matrices consistently.
    if (ignoreunclassified){
        dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified", paste0(analysisname, "__Unclassified"))
        keep <- which(!(rownames(PPMmat) %in% dunno) & rownames(PPMmat) != "" & !is.na(rownames(PPMmat)))
        PPMmat <- PPMmat[keep, , drop = FALSE]
        BaseMat <- BaseMat[rownames(PPMmat), , drop = FALSE]
        if (have_GC) GCmat <- GCmat[rownames(PPMmat), , drop = FALSE]
    }

    nsamp <- ncol(PPMmat)
    nfeat_total <- nrow(PPMmat)

    #IMPORTANT for honest depth accounting: % depth retained is measured against the TOTAL depth of
    #the (vetted) object BEFORE the ignoreunclassified drop would understate the denominator. We use
    #the per-sample total of the kept feature universe as the denominator, i.e. the same feature set
    #the grid draws from, so retention fractions are internally consistent with the count panel.
    sample_total_depth <- colSums(BaseMat)
    #Guard against any all-zero sample column (would divide by zero).
    sample_total_depth[sample_total_depth == 0] <- NA_real_

    #Resolve the prevalence rule: absolute count takes precedence over percentage.
    if (!is.null(prevalence_min_samples)){
        min_samples <- max(1, min(prevalence_min_samples, nsamp))
        prev_desc <- paste0(">= ", min_samples, " sample(s)")
    } else {
        min_samples <- (min(prevalence_pct, 100) / 100) * nsamp
        prev_desc <- paste0(">= ", prevalence_pct, "% of samples")
    }

    #Auto-build PPM ladder if not supplied: 0 plus a roughly log-spaced ladder up to a high
    #percentile of the per-feature maximum PPM, so the informative low end is well sampled.
    if (is.null(PPM_thresholds)){
        feat_max_ppm <- apply(PPMmat, 1, max)
        top <- as.numeric(quantile(feat_max_ppm[feat_max_ppm > 0], 0.95, names = FALSE))
        if (!is.finite(top) || top <= 0) top <- max(1, max(PPMmat))
        #Log-spaced from ~1 to top, prepended with 0.
        ladder <- unique(round(c(0, exp(seq(log(1), log(top), length.out = max(1, n_ppm - 1))))))
        PPM_thresholds <- sort(ladder)
    }

    #Auto-build completeness ladder if not supplied: evenly spaced 0..max observed (capped at 100).
    if (is.null(completeness_thresholds)){
        if (have_GC){
            gcmax <- min(100, max(GCmat, na.rm = TRUE))
            if (!is.finite(gcmax) || gcmax <= 0) gcmax <- 100
            completeness_thresholds <- unique(round(seq(0, gcmax, length.out = n_comp)))
        } else {
            completeness_thresholds <- 0
        }
    }
    if (!have_GC && any(completeness_thresholds > 0)){
        if (verbose) flog.warn("No GenomeCompleteness assay present; completeness thresholds > 0 collapsed to 0.")
        completeness_thresholds <- 0
    }

    #For a single (ppm, comp) cell, return the surviving feature count AND the per-sample fraction
    #of sequencing depth those survivors represent. The depth part is what reveals whether a big
    #drop in feature COUNT is cheap (axing low-abundance noise) or expensive (axing real community);
    #these can diverge enormously depending on each sample's evenness.
    survivors_stats <- function(ppm_thr, comp_thr){
        pass_ppm <- if (ppm_thr <= 0) rep(TRUE, nfeat_total) else (rowSums(PPMmat >= ppm_thr) >= min_samples)
        if (have_GC && comp_thr > 0){
            pass_gc <- (rowSums(GCmat >= comp_thr) >= min_samples)
        } else {
            pass_gc <- rep(TRUE, nfeat_total)
        }
        surv <- pass_ppm & pass_gc
        n_surv <- sum(surv)
        #Per-sample retained depth fraction: depth of survivors / total depth, per column.
        if (n_surv > 0){
            retained_depth <- colSums(BaseMat[surv, , drop = FALSE])
        } else {
            retained_depth <- rep(0, nsamp)
        }
        depth_frac <- retained_depth / sample_total_depth  #vector length nsamp, NA where sample empty
        list(n = n_surv,
             depth_median = stats::median(depth_frac, na.rm = TRUE),
             depth_min = suppressWarnings(min(depth_frac, na.rm = TRUE)))
    }

    grid <- expand.grid(PPM = sort(unique(PPM_thresholds)), Completeness = sort(unique(completeness_thresholds)), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    stats_list <- mapply(survivors_stats, grid$PPM, grid$Completeness, SIMPLIFY = FALSE)
    grid$n_features <- sapply(stats_list, function(x) x$n)
    grid$pct_retained <- round(100 * grid$n_features / nfeat_total, 1)
    #Depth retention summaries across samples (as percentages).
    grid$pct_depth_median <- round(100 * sapply(stats_list, function(x) x$depth_median), 1)
    grid$pct_depth_min <- round(100 * sapply(stats_list, function(x) { v <- x$depth_min; if (is.finite(v)) v else NA_real_ }), 1)
    grid$facet <- facet_label

    if (verbose) flog.info(sprintf("%s: counted survivors over %d PPM x %d completeness cells (%d starting features, %d samples).", facet_label, length(unique(grid$PPM)), length(unique(grid$Completeness)), nfeat_total, nsamp))

    #Local rate of decline: fractional feature loss relative to the looser (one-step-down) neighbour
    #in each direction, averaged. Small magnitude = a stable region (tightening barely costs features);
    #large magnitude = a steep face (each notch removes many features). This is where turning points
    #such as the ~10% completeness floor become visible as a drop-off in the rate.
    ppm_levels <- sort(unique(grid$PPM))
    comp_levels <- sort(unique(grid$Completeness))
    key <- function(p, c) paste(p, c, sep = "§")
    lut <- setNames(grid$n_features, key(grid$PPM, grid$Completeness))
    grid$decline_rate <- NA_real_
    for (i in 1:nrow(grid)){
        p <- grid$PPM[i]; c <- grid$Completeness[i]
        here <- lut[[key(p, c)]]
        rates <- c()
        pi <- match(p, ppm_levels); ci <- match(c, comp_levels)
        if (pi > 1){ prev <- lut[[key(ppm_levels[pi - 1], c)]]; if (prev > 0) rates <- c(rates, (prev - here) / prev) }
        if (ci > 1){ prev <- lut[[key(p, comp_levels[ci - 1])]]; if (prev > 0) rates <- c(rates, (prev - here) / prev) }
        if (length(rates) > 0) grid$decline_rate[i] <- mean(rates)
    }

    if (!return_plots){
        return(invisible(grid))
    }

    gdf <- grid
    gdf$PPM_f <- factor(gdf$PPM, levels = ppm_levels)
    gdf$Comp_f <- factor(gdf$Completeness, levels = comp_levels)

    tile_theme <- theme_minimal() + theme(plot.title = element_text(size = 9), axis.text = element_text(size = 7), axis.text.x = element_text(angle = 45, hjust = 1))
    subt <- paste0(facet_label, " | prevalence ", prev_desc, " | ", nfeat_total, " starting features | n = ", nsamp, " samples")

    p_count <- ggplot(gdf, aes(x = PPM_f, y = Comp_f, fill = n_features)) +
        geom_tile(colour = "grey80") +
        geom_text(aes(label = sprintf("%d\n%.0f%%", n_features, pct_retained)), size = 2.3) +
        scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
        labs(x = "min PPM threshold", y = "min genome completeness (%)", fill = "features", title = paste0("Surviving feature count (count and % of starting retained)\n", subt)) +
        tile_theme

    #Decline-rate surface: dark = steep (many features lost per notch), pale = stable plateau.
    p_rate <- ggplot(gdf, aes(x = PPM_f, y = Comp_f, fill = decline_rate)) +
        geom_tile(colour = "grey80") +
        geom_text(aes(label = ifelse(is.na(decline_rate), "", sprintf("%.0f%%", 100 * decline_rate))), size = 2.3) +
        scale_fill_gradient(low = "#f7f7f7", high = "#b30000", na.value = "grey90", labels = scales::percent) +
        labs(x = "min PPM threshold", y = "min genome completeness (%)", fill = "feature loss\nper step", title = paste0("Local rate of feature loss (pale = stable plateau; dark = steep)\n", subt)) +
        tile_theme

    #Depth-retention surface: what % of sequencing depth the survivors represent. Median across
    #samples (fill and top number), with the WORST-hit sample's retention beneath (min). A cell can
    #retain few FEATURES yet high DEPTH (axed only noise) or, alarmingly, low depth (axed real
    #community). The min matters because a sample retaining, say, 15% depth is effectively gutted
    #even if the median looks comfortable.
    p_depth <- ggplot(gdf, aes(x = PPM_f, y = Comp_f, fill = pct_depth_median)) +
        geom_tile(colour = "grey80") +
        geom_text(aes(label = sprintf("%.0f%%\n(min %.0f%%)", pct_depth_median, pct_depth_min)), size = 2.2) +
        scale_fill_gradient(low = "#fff5eb", high = "#7f2704", limits = c(0, 100)) +
        labs(x = "min PPM threshold", y = "min genome completeness (%)", fill = "% depth\n(median)", title = paste0("Sequencing-depth retained by survivors (median across samples; worst sample in parens)\n", subt)) +
        tile_theme

    #Divergence surface: the gap between feature loss and depth loss. This is the direct read on the
    #evenness question. depth_retained - feature_retained, in percentage points:
    #  large POSITIVE (blue) = kept far more depth than features => you axed low-abundance noise (good);
    #  near ZERO (white)     = features and depth fell together => you axed a representative slice;
    #  NEGATIVE (red)        = lost more depth than features => you axed high-abundance community (bad).
    gdf$divergence <- gdf$pct_depth_median - gdf$pct_retained
    p_diverge <- ggplot(gdf, aes(x = PPM_f, y = Comp_f, fill = divergence)) +
        geom_tile(colour = "grey80") +
        geom_text(aes(label = sprintf("%+.0f", divergence)), size = 2.3) +
        scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0) +
        labs(x = "min PPM threshold", y = "min genome completeness (%)", fill = "depth% -\nfeature%", title = paste0("Depth retained minus features retained (blue = axed mostly noise; red = axed real depth)\n", subt)) +
        tile_theme

    titrate_feature_counts_list <- list(grid = grid, plots = list(feature_count = p_count, decline_rate = p_rate, depth_retained = p_depth, count_vs_depth_divergence = p_diverge))

    if (plot_on_the_fly){
        print(titrate_feature_counts_list$plots)
        return(NULL)
    } else {
        return(titrate_feature_counts_list)
    }

}