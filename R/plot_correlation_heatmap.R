#' plot_correlation_heatmap(ExpObj = NULL, glomby = NULL, stattype = "spearman", subsetby = NULL, maxnumfeatallowed = 10000, minabscorrcoeff = NULL, ntopvar = NULL, featuresToKeep = NULL, only_allow_CSBs = FALSE, samplesToKeep = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, show_GenomeCompleteness_boxplot = TRUE, PPM_normalize_to_bases_sequenced = FALSE, normalization = "relabund", showGram = FALSE, showPhylum = FALSE, addtit = NULL, cdict = NULL, ignoreunclassified = TRUE, class_to_ignore = "N_A", returnstats = FALSE)
#'
#' Plots correlation heatmaps annotated by the metadata or a correlelogram of features
#' @export

plot_correlation_heatmap <- function(ExpObj = NULL, glomby = NULL, stattype = "spearman", subsetby = NULL, maxnumfeatallowed = 10000, minabscorrcoeff = NULL, ntopvar = NULL, featuresToKeep = NULL, only_allow_CSBs = FALSE, samplesToKeep = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, show_GenomeCompleteness_boxplot = TRUE, PPM_normalize_to_bases_sequenced = FALSE, normalization = "relabund", showGram = FALSE, showPhylum = FALSE, addtit = NULL, cdict = NULL, ignoreunclassified = TRUE, class_to_ignore = "N_A", returnstats = FALSE) {

    #Account for JAMS2 spaces
    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")

    #Validate normalization choice up front so a typo fails loudly rather than silently
    #falling back to BaseCounts later.
    if (!(normalization %in% c("relabund", "clr"))){
        stop("normalization must be either \"relabund\" (relative abundance in PPM) or \"clr\" (centred log-ratio).")
    }

    #CRITICAL FIX for only_allow_CSBs + glomby:
    #ExpObjVetting applies only_allow_CSBs BEFORE glomby. But CSB restriction is only
    #valid on a "ConsolidatedGenomeBin" analysis (see filter_experiment), and CSB row
    #names carry the "CSB_" prefix, which is stripped/collapsed on agglomeration. So we
    #must (1) restrict to CSBs while the object is still ConsolidatedGenomeBin, and only
    #then (2) agglomerate. ExpObjVetting already orders these operations that way, so we
    #simply pass both through. The historical crash came later, from annotation blocks
    #assuming a Phylum column that no longer exists post-glom. Those are made defensive below.
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, only_allow_CSBs = only_allow_CSBs, glomby = glomby, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    #Phylum annotation only makes sense when the effective taxonomic level being plotted
    #is at Class or below (down to IS1/terminal bins). Abrogate with a warning otherwise.
    phylum_annot_levels <- c("Class", "Order", "Family", "Genus", "Species", "IS1", "LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")
    if (showPhylum && !(analysisname %in% phylum_annot_levels)){
        flog.warn(paste0("showPhylum = TRUE is not meaningful at the '", analysisname, "' level (Phylum annotation requires a sub-Class/parent-Phylum resolution). Turning showPhylum off for this plot."))
        showPhylum <- FALSE
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff)

    #Multi-tiered subsetting, harmonized with plot_relabund_heatmap / plot_Ordination.
    if (!(is.null(subsetby))){
        subset_list <- multiple_subsetting_sample_selector(SEobj = obj, phenotable = NULL, subsetby = subsetby, compareby = NULL, cats_to_ignore = class_to_ignore)
        #Guard against a NULL return (e.g. non-discrete subset variables)
        if (is.null(subset_list)){
            flog.warn("Subsetting could not be resolved (check that subsetby variables are discrete). Defaulting to no subsetting.")
            subset_points <- "none"
            subsetby <- NULL
        } else {
            subset_df <- subset_list$Subsets_stats
            subset_df <- subset_df[which(subset_df$Subset_Tier_Level != 0), , drop = FALSE]
            if (any(subset_df$Num_samples_in_subset < 2)){
                LowSampSubsets <- subset_df[which(subset_df$Num_samples_in_subset < 2), "Subset_Tier_Class_Name"]
                flog.warn(paste("Subsets", paste0(LowSampSubsets, collapse = ", "), "contain fewer than 2 samples, and will thus be omitted."))
                subset_df <- subset_df[which(subset_df$Num_samples_in_subset >= 2), , drop = FALSE]
            }
            if (nrow(subset_df) < 1){
                flog.warn("There are no surviving subset points. Defaulting to no subsetting.")
                subset_points <- "none"
                subsetby <- NULL
            } else {
                subset_points <- subset_df$Subset_Tier_Class_Name
            }
        }
    } else {
        subset_points <- "none"
    }

    #Initialize Stats and Graph Vector lists
    svec <- list()
    s <- 1
    n <- 1

    #subset by metadata column
    for (sp in 1:length(subset_points)) {

        if (!(is.null(subsetby))){
            samplesToKeep_sp <- subset_list[[subset_points[sp]]]
            flog.info(paste("Plotting within", subset_points[sp]))
            subsetname <- subset_points[sp]
        } else {
            samplesToKeep_sp <- rownames(colData(obj))
            subsetname <- "no_sub"
        }

        #Need at least 2 samples for correlation.
        if (length(samplesToKeep_sp) < 2){
            flog.warn(paste("Fewer than 2 samples within", subsetname, "- skipping."))
            next
        }

        currobj <- filter_experiment(SEobj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = samplesToKeep_sp, featuresToKeep = featuresToKeep, only_allow_CSBs = FALSE, normalization = normalization, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff)

        numfeats <- nrow(currobj)

        #There must be at least four features for a meaningful correlation heatmap
        if (nrow(currobj) > 3){
            #Compose an appropriate title for the plot
            if (length(unique(subset_points)) > 1){
                maintit <- paste("Feature Correlation Heatmap", analysisname, paste("within", subset_points[sp]), sep = " | ")
            } else {
                maintit <- paste("Feature Correlation Heatmap", analysisname, sep = " | ")
            }
            if (!is.null(addtit)) {
                maintit <- paste(addtit, maintit, sep = "\n")
            }

            #Get counts matrix. Pick the assay matching the requested normalization.
            #filter_experiment attaches $PPM when normalization includes "relabund" (or a
            #PPM/SD filter fires) and $CLR when it includes "clr". Guard against either being
            #absent and fall back to BaseCounts with a warning.
            if (normalization == "clr"){
                if ("CLR" %in% names(assays(currobj))){
                    countmat <- as.matrix(assays(currobj)$CLR)
                } else {
                    flog.warn("CLR assay not found on the filtered object; falling back to BaseCounts. Check that normalization = \"clr\" was honoured.")
                    countmat <- as.matrix(assays(currobj)$BaseCounts)
                }
            } else {
                if ("PPM" %in% names(assays(currobj))){
                    countmat <- as.matrix(assays(currobj)$PPM)
                } else {
                    flog.warn("PPM assay not found on the filtered object; falling back to BaseCounts. Set normalization = \"relabund\" for relative-abundance correlations.")
                    countmat <- as.matrix(assays(currobj)$BaseCounts)
                }
            }

            if (ignoreunclassified == TRUE){
                dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified", paste0(analysisname, "__Unclassified"))
                rowsToKeep <- which(!(rownames(countmat) %in% dunno))
                countmat <- countmat[rowsToKeep, , drop = FALSE]
            }

            #Protect against empty rows that would poison correlation. For relative-abundance
            #matrices an all-zero row is empty and safe to drop; for CLR, values can be
            #negative and a row need not sum to >0 to be informative, so only drop rows with
            #no variance (which cannot be correlated) plus unnamed/NA rows.
            if (normalization == "clr"){
                rowsToKeep <- which(rowSds(countmat) > 0 & rownames(countmat) != "" & !is.na(rownames(countmat)))
            } else {
                rowsToKeep <- which(rowSums(countmat) > 0 & rownames(countmat) != "" & !is.na(rownames(countmat)))
            }
            countmat <- countmat[rowsToKeep, , drop = FALSE]

            #Rename rows to include description if not taxonomic data
            if (!analysis %in% taxonomic_spaces){
                feattable <- rowData(currobj)
                feattable$Feature <- paste(feattable$Accession, feattable$Description, sep = "-")
                rownames(countmat) <- feattable$Feature[match(rownames(countmat), feattable$Accession)]
            }
            matrixSamples <- colnames(countmat)
            matrixRows <- rownames(countmat)

            if (!is.null(ntopvar)){
                ntop <- min(ntopvar, nrow(countmat))
                featsds <- rowSds(countmat)
                featIndices <- names(featsds[order(featsds, decreasing = TRUE)[1:ntop]])
                countmat <- countmat[featIndices, , drop = FALSE]
                ntopvarmsg <- paste("Top", ntop, "most variant features across samples")
            } else {
                ntopvarmsg <- NULL
            }

            #After all pruning, re-check there is enough to correlate.
            if (nrow(countmat) < 4){
                flog.warn(paste("Fewer than 4 features remain for", subsetname, "after filtration - skipping."))
                next
            }

            docorrelations <- TRUE
            if (!(is.null(maxnumfeatallowed))) {
                if (nrow(countmat) > maxnumfeatallowed){
                    flog.warn(paste("There are", nrow(countmat), "features to pairwise correlate, which is more than", maxnumfeatallowed, "allowed. This would entail", (nrow(countmat) ^ 2), "comparisons. If you are sure you want that many, set maxnumfeatallowed to a higher value."))
                    docorrelations <- FALSE
                }
            }

            if (docorrelations == TRUE){
                #Calculate matrix stats and get new matrix with correlations.
                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = FALSE, statsonlog = FALSE, stattype = stattype, classesdf = NULL)

                if (!is.null(minabscorrcoeff)){
                    flog.info(paste("Eliminating features which do not correlate with other features with a coefficient of at least", minabscorrcoeff))
                    matstats <- filter_correlations(corrmat = matstats, mincorrelcoeff = minabscorrcoeff)
                    minabscorrcoeffmsg <- paste("Largest correlation coefficient at least", minabscorrcoeff)
                } else {
                    minabscorrcoeffmsg <- NULL
                }

                #After correlation filtering there must still be a plottable matrix.
                if (is.null(matstats) || nrow(matstats) < 2 || ncol(matstats) < 2){
                    flog.warn(paste("Fewer than 2 features survive correlation filtering for", subsetname, "- skipping."))
                    next
                }

                #Plot heatmap - set colour scale
                CorrHmColours <- c("blue4", "lightgoldenrodyellow", "red1")
                heatmapCols <- colorRamp2(c(-1, 0, 1), CorrHmColours)

                fontcoefficient <- (-0.05 * nrow(matstats)) + 7.5
                fontsizey <- round((((-1 / 300) * (nrow(matstats))) + 0.85 * fontcoefficient), 2)
                fontsizey <- max(0.5, fontsizey)

                #Add annotations if taxonomic. Gram and Phylum are now resolved via the shared
                #resolve_tax_annotation_colours() helper (per-feature Gram, PhyCols-based Phylum
                #colours), requested independently via showGram / showPhylum. The genome-completeness
                #boxplot remains governed by show_GenomeCompleteness_boxplot. Everything stays
                #defensive: any component that cannot resolve degrades to no annotation.
                ha1 <- NULL
                ha2 <- NULL

                if (analysis %in% taxonomic_spaces){

                    #Genome completeness boxplot list, only if the assay is present.
                    gcl <- NULL
                    have_GC <- "GenomeCompleteness" %in% names(assays(currobj))
                    if (have_GC){
                        genomecompletenessdf <- assays(currobj)$GenomeCompleteness
                        gc_feats <- rownames(matstats)[rownames(matstats) %in% rownames(genomecompletenessdf)]
                        if (length(gc_feats) > 0){
                            genomecompletenessstats <- as.matrix(genomecompletenessdf[gc_feats, , drop = FALSE])
                            gcl <- lapply(1:nrow(genomecompletenessstats), function (x){
                                vals <- as.numeric(genomecompletenessstats[x, ])
                                vals <- vals[which(vals != 0)]
                                #anno_boxplot cannot cope with a zero-length vector; give it a 0.
                                if (length(vals) == 0) vals <- 0
                                return(vals)
                            })
                        } else {
                            have_GC <- FALSE
                        }
                    }

                    #Resolve Gram/Phylum colours via the shared helper, each requested independently.
                    annot <- resolve_tax_annotation_colours(feature_table = rowData(currobj), row_order = rownames(matstats), want_phylum = showPhylum, want_gram = showGram)
                    have_phylo_annot <- any(c(!is.null(annot$phylum), !is.null(annot$gram)))

                    #Helper to append whichever Gram/Phylum components resolved onto an argument list.
                    append_tax_annots <- function(arglist){
                        cols <- arglist$col
                        if (is.null(cols)) cols <- list()
                        if (!is.null(annot$gram)){
                            arglist$Gram <- annot$gram
                            cols$Gram <- annot$gram_cols
                        }
                        if (!is.null(annot$phylum)){
                            arglist$Phylum <- annot$phylum
                            cols$Phylum <- annot$phylum_cols
                        }
                        if (length(cols) > 0){
                            arglist$col <- cols
                        }
                        return(arglist)
                    }

                    #Build the LEFT annotation. Assemble conditionally so we never hand a NULL
                    #vector to anno_boxplot / HeatmapAnnotation.
                    if (show_GenomeCompleteness_boxplot && have_GC){
                        left_args <- list(
                            Pct_Genome_Compl = anno_boxplot(gcl, width = unit(4, "cm"), pch = 20, size = unit(1, "mm"), axis_param = list(labels_rot = 90))
                        )
                        if (have_phylo_annot){
                            left_args <- append_tax_annots(left_args)
                        }
                        left_args$annotation_name_gp <- gpar(fontsize = 6, col = "black")
                        ha1 <- do.call(rowAnnotation, left_args)
                    } else if (have_phylo_annot){
                        #No boxplot, but we can still show Gram/Phylum on the left.
                        left_args <- append_tax_annots(list(annotation_name_gp = gpar(fontsize = 6, col = "black")))
                        ha1 <- do.call(rowAnnotation, left_args)
                    }

                    #Bottom annotation mirrors Gram/Phylum, only when available.
                    if (have_phylo_annot){
                        bottom_args <- append_tax_annots(list(annotation_name_gp = gpar(fontsize = 6, col = "black"), show_legend = FALSE))
                        ha2 <- do.call(HeatmapAnnotation, bottom_args)
                    }
                }

                #Build plot title
                normspacemsg <- switch(normalization, "clr" = "CLR-transformed counts", "relabund" = "Relative abundance (PPM)", normalization)
                stattit <- paste0("Correlation measure = ", stattype, " | ", normspacemsg)
                plotit <- paste(maintit, stattit, presetlist$filtermsg, ntopvarmsg, minabscorrcoeffmsg, sep = "\n")

                svec[[s]] <- matstats
                stattitle <- paste(analysisname, stattype, subsetname, sep = "_")
                names(svec)[s] <- stattitle
                s <- s + 1

                ht1 <- Heatmap(matstats, name = paste(stattype, "correlation coefficient"), column_title = plotit, column_title_gp = gpar(fontsize = 10), col = heatmapCols, column_dend_height = unit(5, "mm"), cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = fontsizey), column_names_gp = gpar(fontsize = fontsizey), heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm"), title = paste(stattype, "correlation coefficient"), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)), left_annotation = ha1, bottom_annotation = ha2)

                draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "left", padding = unit(c(2, 20, 2, 2), "mm"))
                n <- n + 1
            } #End conditional of going ahead and doing correlations
        } else {
            flog.warn(paste("Fewer than 4 features left over after filtering for", subsetname, "- skipping."))
        } #End conditional if there are any features left over after filtering
    } #End for loop for plotting within each subset point

    if (returnstats == TRUE){
        return(svec)
    } else {
        return(print("Heatmaps generated."))
    }
}