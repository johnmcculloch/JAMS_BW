#' calculate_matrix_stats(countmatrix = NULL, uselog = NULL, statsonlog = TRUE, stattype = NULL, classesdf = NULL, invertbinaryorder = FALSE, fun_for_l2fc = "geom_mean", threshPA = 0, numthreads = 1, nperm = 10000, paired = FALSE)
#'
#' Returns a data frame with the statistics for a feature count matrix ordered by highest variance or lowest Mann-Whitney-Wilcoxon test between binary categories.
#' This function is used internally by the JAMS plotting functions, plot_Ordination, plot_relabund_heatmap, plot_correlation_heatmap and plot_relabund_features.
#'
#' @export

calculate_matrix_stats <- function(countmatrix = NULL, uselog = NULL, statsonlog = TRUE, stattype = NULL, classesdf = NULL, invertbinaryorder = FALSE, fun_for_l2fc = "geom_mean", threshPA = 0, numthreads = 1, nperm = 10000, paired = FALSE){

    #Test for silly stuff
    if ((stattype %in% c("binary", "permanova", "anova", "PA")) && (is.null(classesdf))){
        stop("If rows are to be selected by highest significant difference between classes in a discrete category, you must determine the category using the argument *classesdf*. See the documentation for make_classes_df")
    }

    #Decide if we are comparing anything
    if (!(is.null(classesdf))){
        #Decide if variable is continuous or discrete
        if (can_be_made_numeric(classesdf$cl)){
            numclass <- "continuous"
            classesdf$cl <- as.numeric(classesdf$cl)
        } else {
            #Determine if there are two or more classes
            numclass <- length(unique(classesdf$cl))
        }
        classesvector <- classesdf$cl
        #Ensure same order of samples on counts matrix and classesdf. Should be the same already, but better safe than sorry.
        if (nrow(countmatrix) > 1){
            countmatrix <- countmatrix[ , rownames(classesdf)]
        }

    } else {
        #No comparing, use variance
        numclass <- 0
        classesvector <- NULL
    }

    #If auto, decide what to do
    if (stattype %in% c("auto", "PA")){
        if (is.null(classesdf)){
            stattype <- "variance"
        } else {
            if(numclass > 2){
                stattype <- "anova"
            } else if (numclass == 2) {
                if (stattype != "PA"){
                    stattype <- "binary"
                }
            } else {
                stattype <- "variance"
            }
        }
    }

    if (all(c((uselog == TRUE), (statsonlog == FALSE)))){
        flog.info("Transforming log2 counts back to raw counts for calculating stats.")
        #log2 transform if applicable
        countmatrix <- convert_matrix_log2(mat = countmatrix, transformation = "from_log2")
    }

    #THE MULTIPLE-TESTING DENOMINATOR MUST BE THE NUMBER OF TESTS ACTUALLY PERFORMED.
    #
    #Every branch below ends in p.adjust(matstats$pval, method = x, n = ...). The n passed used to
    #be length(matstats$pval), i.e. one slot for every ROW, including rows whose statistic could
    #not be computed and came back NA. R's p.adjust excludes NA p-values from the correction but
    #honours whatever n it is given, so those rows were silently inflating every adjusted p-value.
    #A row that could not be tested is not a hypothesis, and must not be counted as one.
    #
    #THIS IS NOT A SMALL CORRECTION. The commonest untestable row is a feature that is absent from
    #every sample of the comparison, and such rows are abundant in two ordinary situations which
    #compound each other. First, an object accumulates features that are zero throughout whenever
    #samples are dropped from it - ExpObjVetting(samplesToKeep), filter_sample_by_class_to_ignore
    #and filter_experiment all take a pure column slice and leave behind every feature that
    #belonged only to the samples that went. A 178-sample multi-cohort human stool object built
    #this way was 68.8% such rows (32124 of 46691). Second, and independently, the subsetby loops
    #in plot_relabund_heatmap and plot_volcano, and every panel of compare_within_subsets, pin
    #featcutoff = c(0, 0) so that the feature set stays comparable across subsets - which means no
    #filtration removes the rows that are empty within the subset being tested.
    #
    #Measured on that object, the factor by which the Benjamini-Hochberg denominator was inflated:
    #  178 samples -> 3.3x     80 -> 5.1x     40 -> 7.1x     20 -> 10.9x     10 -> 16.2x
    #Every adjusted p-value in every subset panel was too large by that much, and the smaller the
    #subset the worse it got, which is the opposite of what anyone would assume.
    #
    #THE ROWS THEMSELVES ARE DELIBERATELY NOT DELETED. The old commented-out guard that used to
    #sit here, #rowsToKeep <- which(rowSums(countmatrix) > 0), was right to be disabled and must
    #not be reinstated. The matrix handed in is usually PPM (plot_relabund_heatmap:321) and getPPM
    #rounds to whole PPM, so a real feature carrying thousands of bases but sitting under 0.5 PPM
    #in every sample is an all-zero PPM row; deleting those would impose a silent half-PPM
    #abundance floor on every heatmap, volcano and correlation plot in the package. It is also
    #meaningless on a CLR matrix, where an absent feature is not zero at all and a present one can
    #sum below zero. And it would delete rows from matstats, which callers index by name:
    #plot_relabund_features would silently drop explicitly requested features, and plot_volcano's
    #drop_features_absent_in_all_samples switch would become a no-op in its FALSE direction while
    #its log still claimed parity with the heatmap. Correcting n fixes the defect and changes
    #nothing else - every row still comes back, carrying its NA where no test was possible.
    adjust_pvalues <- function(matstats = NULL){
        n_tested <- sum(!is.na(matstats$pval))
        n_untestable <- (nrow(matstats) - n_tested)
        if (n_untestable > 0){
            if (n_tested < 1){
                flog.warn(paste0("calculate_matrix_stats: not one of the ", nrow(matstats), " features yielded a computable p-value. No correction can be made and every adjusted p-value will be NA. This normally means the comparison is degenerate - check for constant or empty features."))

                n_tested <- nrow(matstats)
            } else {
                flog.warn(sprintf("calculate_matrix_stats: %d of %d features (%.1f%%) yielded no computable p-value, almost always because they are absent, or constant, across every sample in this comparison. They are reported with NA statistics but are EXCLUDED from the multiple-testing denominator, which is therefore %d and not %d. Counting them would have inflated every adjusted p-value by a factor of %.2f. Rows that are empty throughout accumulate whenever samples are dropped from an object, so they are commonest in subset panels.", n_untestable, nrow(matstats), ((100 * n_untestable) / nrow(matstats)), n_tested, nrow(matstats), (nrow(matstats) / n_tested)))
            }
        }
        adjlist <- lapply(p.adjust.methods, function(x){ p.adjust(matstats$pval, method = x, n = n_tested) })
        names(adjlist) <- paste("padj", p.adjust.methods, sep = "_")

        return(as.data.frame(adjlist))
    }

    #Calculate matrix stats and get new matrix.
    if (stattype == "variance"){

        flog.info("Calculating variance across samples.")
        featStatsSD <- apply(countmatrix, 1, sd)
        featStatsMAD <- apply(countmatrix, 1, mad)
        matstats <- data.frame(SD = as.vector(unlist(featStatsSD)), MAD = as.vector(unlist(featStatsMAD)))
        rownames(matstats) <- rownames(countmatrix)
        matstats <- matstats[order(matstats$SD, decreasing = TRUE), ]
        matstats$Method <- rep("variance", nrow(matstats))
        #Add rank variance to matrix names
        matstats$Rank <- sapply(1:nrow(matstats), function(x) { getOrdinalNumber1(x) })

    } else if (stattype == "binary"){
        if (numclass == 2){
            flog.info("Calculating p-values with Mann-Whitney-Wilcoxon test.")
            flog.info(paste("Log2 fold change will be calculated with the", fun_for_l2fc, "of values within each group."))
        } else {
            stop("To calculate p-values with Mann-Whitney-Wilcoxon test, the comparison variable must have exactly two classes.")
        }

        discretenames <- sort(unique(classesdf$cl))
        pair_wilcox_test <- ("wilcox_pairs" %in% colnames(classesdf))

        #When using paired give a global warning that it may be impossible to calculate exact p-values when there are zeroes, and suppress the warnings for each event.
        if (pair_wilcox_test){
            flog.warn("Using paired test")
            suppressWarnings(comparisons <- lapply(1:nrow(countmatrix), function(rw) { wilcox.test(x = countmatrix[rw, classesdf$Sample[which(classesdf$cl == discretenames[1])]], y = countmatrix[rw, classesdf$Sample[which(classesdf$cl == discretenames[2])]], paired = pair_wilcox_test, exact = TRUE) }))
        } else {
            flog.warn("Using unpaired test")
            #comparisons <- lapply(1:nrow(countmatrix), function(rw) { wilcox.test(x = countmatrix[rw, classesdf$Sample[which(classesdf$cl == discretenames[1])]], y = countmatrix[rw, classesdf$Sample[which(classesdf$cl == discretenames[2])]], paired = pair_wilcox_test) })
            #Revert to formula notation, as it is much, much, much faster, unless the df has less than 2 rows (like when agregating features in plot_relabund_features)
            if (nrow(countmatrix) > 1){
                comparisons <- lapply(1:nrow(countmatrix), function(rw) { wilcox.test(as.numeric(countmatrix[rw, classesdf$Sample]) ~ classesdf$cl) })
            } else {
                #Counts is a single row
                comparisons <- list()
                comparisons[[1]] <- wilcox.test(as.numeric(countmatrix[ , rownames(classesdf)]) ~ classesdf$cl)
            }
        }

        mwstat <- sapply(1:length(comparisons), function(x) { unname(comparisons[[x]]$statistic) } )
        mwpval <- sapply(1:length(comparisons), function(x) { unname(comparisons[[x]]$p.value) } )

        #Get Log2 Fold Change (l2fc)
        getl2fc <- function(countsvec = NULL, classesvector = NULL, discretenames= NULL, method = NULL, countsinlog = NULL){
            #If counts are in log space, transform them back
            if (countsinlog == TRUE){
                countsvecNL <- as.numeric(sapply(countsvec, function (x) { ((2 ^ x) - 1) }))
            } else {
                countsvecNL <- as.numeric(countsvec)
            }
            countsvecG1 <- countsvecNL[which(classesvector == discretenames[1])]
            countsvecG2 <- countsvecNL[which(classesvector == discretenames[2])]
            if(method == "geom_mean"){
                statscountsvecG1 <- exp(mean(log(as.numeric(countsvecG1 + 1))))
                statscountsvecG2 <- exp(mean(log(as.numeric(countsvecG2 + 1))))
            } else {
                statscountsvecG1 <- get(method)(countsvecG1)
                statscountsvecG2 <- get(method)(countsvecG2)
            }

            #Same arithmetic as gtools::foldchange2logratio(gtools::foldchange(num, denom), base = 2), inlined to drop gtools.
            fc <- ifelse(statscountsvecG1 < statscountsvecG2, -statscountsvecG2 / statscountsvecG1, statscountsvecG1 / statscountsvecG2)
            l2fc <- log(ifelse(fc < 0, 1 / -fc, fc), base = 2)
            if (is.na(l2fc)){
                l2fc <- 0
            }

            return(l2fc)
        }

        l2fc <- sapply(1:nrow(countmatrix), function(x) { getl2fc(countsvec = countmatrix[x, ], classesvector = classesvector, discretenames = discretenames, countsinlog = uselog, method = fun_for_l2fc) })

        if (invertbinaryorder == TRUE){
            l2fc <- (l2fc * -1)
        }
        matstats <- data.frame(pval = as.vector(unlist(mwpval)), stat = as.vector(unlist(mwstat)), l2fc = l2fc)
        rownames(matstats) <- rownames(countmatrix)
        matstats$absl2fc <- abs(matstats$l2fc)

        #Adjust p-values
        adjdf <- adjust_pvalues(matstats = matstats)
        matstats <- cbind(matstats, adjdf)
        matstats <- matstats[order(matstats$pval, decreasing = FALSE), ]
        matstats$Method <- rep("MannWhitneyWilcoxon", nrow(matstats))

        #Add binary directionality information
        matstats$Class_where_increased <- NA
        matstats$Class_where_increased[which(matstats$l2fc == 0)] <- "neither"
        matstats$Class_where_increased[which(matstats$l2fc > 0)] <- discretenames[1:2][c(!invertbinaryorder, invertbinaryorder)]
        matstats$Class_where_increased[which(matstats$l2fc < 0)] <- discretenames[1:2][c(invertbinaryorder, !invertbinaryorder)]

    } else if (stattype %in% c("anova", "permanova")){
        require(parallel)

        if (stattype == "permanova"){
            flog.info("Calculating p-values with permanova. Please be patient...")
            flog.info(paste("Calculating p-values with PERMANOVA using", nperm, "permutations on", numthreads, "threads. Please be patient..."))
            #vegan::adonis() is defunct since vegan 2.7-1 (2025); adonis2(by = "terms") gives the same per-term test.
            mwpval <- sapply(1:nrow(countmatrix), function(x) { vegan::adonis2(countmatrix[x,] ~ classesvector, parallel = numthreads, permutations = nperm, by = "terms")$`Pr(>F)`[1] } )
        } else if (stattype == "anova") {
            flog.info("Calculating p-values with ANOVA.")
            comparisons <- lapply(1:nrow(countmatrix), function(x) { summary(aov(countmatrix[x,] ~ classesvector))[[1]] })
            mwstat <- sapply(1:length(comparisons), function(x) { as.numeric(comparisons[[x]]["classesvector", "F value"]) } )
            mwpval <- sapply(1:length(comparisons), function(x) { as.numeric(comparisons[[x]]["classesvector", "Pr(>F)"]) } )
        }

        mwpval[is.na(mwpval)] <- 1

        matstats <- data.frame(Feature = rownames(countmatrix), pval = as.vector(unlist(mwpval)), stat = as.vector(unlist(mwstat)))

        rownames(matstats) <- rownames(countmatrix)
        #Adjust p-values
        adjdf <- adjust_pvalues(matstats = matstats)
        matstats <- cbind(matstats, adjdf)
        matstats$Feature <- NULL
        matstats <- matstats[order(matstats$pval, decreasing = FALSE), ]
        matstats$Method <- rep(stattype, nrow(matstats))

    } else if (stattype == "PA"){
        if (numclass == 2){
            flog.info(paste("Calculating p-values for presence/absence testing with Fisher's exact test. Threshold for presence is any value >", threshPA))
        } else {
            stop("To calculate p-values for presence/absence testing, the comparison variable must have exactly two classes.")
        }

        discretenames <- sort(unique(classesvector))
        nClass1 = sum(classesvector == discretenames[1])
        nClass2 = sum(classesvector == discretenames[2])

        getPAstats <- function(rownum = NULL, countmatrix = NULL, nClass1 = NULL, nClass2 = NULL, classesvector = NULL, threshPA = threshPA){
            PAvec <- countmatrix[rownum, ]
            if (threshPA != 0){
                PAvec[which(PAvec < threshPA)] <- 0
                PAvec[which(PAvec >= threshPA)] <- 1
            } else {
                PAvec[which(PAvec > 0)] <- 1
            }
            contingencytable <- matrix(ncol = 2, nrow = 2, data = 0)
            contingencytable[1, 1] = sum(PAvec[classesvector == discretenames[1]])
            contingencytable[1, 2] = sum(PAvec[classesvector == discretenames[2]])
            contingencytable[2, 1] = (nClass1 - contingencytable[1, 1])
            contingencytable[2, 2] = (nClass2 - contingencytable[1, 2])
            fisht <- fisher.test(contingencytable, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)
            fishmat <- cbind(OddsRatio = fisht$estimate, OR95lwr = fisht$conf.int[1], OR95upr = fisht$conf.int[2], pval = fisht$p.value)
            rownames(fishmat) <- rownames(countmatrix)[rownum]

            return(fishmat)
        }

        comparisons <- lapply(1:nrow(countmatrix), function(x) { getPAstats(rownum = x, countmatrix = countmatrix, nClass1 = nClass1, nClass2 = nClass2, classesvector = classesvector, threshPA = threshPA) })
        matstats <- do.call(rbind, comparisons)
        matstats <- as.data.frame(matstats)

        #Adjust p-values
        adjdf <- adjust_pvalues(matstats = matstats)
        matstats <- cbind(matstats, adjdf)
        matstats <- matstats[order(matstats$pval, decreasing = FALSE), ]
        if (invertbinaryorder == TRUE){
            matstats$OddsRatio <- (1 / matstats$OddsRatio)
            matstats$OR95lwr <- (1 / matstats$OR95lwr)
            matstats$OR95upr <- (1 / matstats$OR95upr)
        }
        matstats$Method <- rep("fisher", nrow(matstats))

    } else if (stattype %in% c("spearman", "pearson")){

        corstat <- stattype
        #Either correlate to variable or do pairwise feature comparison
        if (is.null(classesvector)){
            use_method <- "pairwise.complete.obs"
            flog.info(paste("Calculating pairwise correlation coefficients of all features using", use_method))
            tcountmatrix <- t(countmatrix)
            matstats <- stats::cor(tcountmatrix, method = corstat, use = use_method)
        } else {
            comparisons <- lapply(1:nrow(countmatrix), function(x) { cor.test(countmatrix[x, ], classesvector, method = corstat) })
            mwstat <- sapply(1:length(comparisons), function(x) { unname(comparisons[[x]]$estimate) } )
            mwpval <- sapply(1:length(comparisons), function(x) { unname(comparisons[[x]]$p.value) } )

            matstats <- data.frame(correl = mwstat, pval = mwpval)
            matstats$abscorrel <- abs(matstats$correl)
            rownames(matstats) <- rownames(countmatrix)
            #Adjust p-values
            adjdf <- adjust_pvalues(matstats = matstats)
            matstats <- cbind(matstats, adjdf)
            matstats <- matstats[order(matstats$pval, decreasing = FALSE), ]
            matstats$Method <- rep(corstat, nrow(matstats))
        }
    }

    return(matstats)
}
