#' Creates ordination plots based on PCA, tSNE or tUMAP

#'plot_Ordination(ExpObj = NULL, glomby = NULL, algorithm = "PCA", PCA_Components = c(1, 2), distmethod = "bray", samplesToKeep = NULL, samplesToHighlight = NULL, featuresToKeep = NULL, subsetby = NULL, compareby = NULL, colourby = NULL, colorby = NULL, shapeby = NULL, use_letters_as_shapes = FALSE, sizeby = NULL, dotsize = 2, connectby = NULL, connection_orderby = NULL, textby = NULL, ellipseby = NULL, tsne_perplx = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, permanova = TRUE, permanova_permutations = 10000, plotcentroids = FALSE, highlight_centroids = TRUE, max_neighbors = 15, show_centroid_distances = FALSE, calculate_centroid_distances_in_all_dimensions = FALSE, addtit = NULL, assay_for_matrix = "BaseCounts", asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, cdict = NULL, grid = TRUE, forceaspectratio = NULL, numthreads = 8, log2tran = TRUE, ignoreunclassified = TRUE, return_coordinates_matrix = FALSE, class_to_ignore = "N_A", ...)

#' @param ExpObj JAMS-style SummarizedExperiment object

#' @param glomby String giving the taxonomic level at which to agglomerate counts. This argument should only be used with taxonomic SummarizedExperiment objects. When NULL (the default), there is no agglomeration

#' @param algorithm String giving the algorithm to be used for dimensionality reduction. Permissible values are "tUMAP", "PCoA", "PCA" or "tSNE". For "tUMAP", the sample-by-feature matrix is processed using Uniform Manifold Approximation and Projection as implemented by the uwot package. For "PCA", the compositional dissimilarity (bray-curtis by default - see distmethod) of the sample-by-feature matrix is processed by Principal Component Analysis as implemented by the stats package. For "PCoA", the compositional dissimilarity (bray-curtis by default - see distmethod) of the sample-by-feature matrix is processed by Principal Component Analysis as implemented by the stats package. For "tSNE", the sample-by-feature matrix is processed using An R t-Distributed Stochastic Neighbor Embedding as implemented by the rtsne package.

#' @param PCA_Components Numerical vector of length 2 specifying which two components to plot on the 2-D ordination plot when using "PCoA" or "PCA" (see algotithm). Default is c(1, 2), meaning the first two components. The variance for each component is included on the axis.

#' @param distmethod String giving the dissimilarity index method for calculating the compositional dissimilarity of the sample-by-feature matrix. Default is "bray" (Bray-Curtis dissimilarity). For permissible values see the vegdist function of the vegan package ("euclidean", "bray", "jaccard", etc...).

#' @param samplesToKeep Vector with sample names to keep. If NULL, all samples within the SummarizedExperiment object are kept. Default is NULL.

#' @param samplesToHighlight Vector with sample names to highlight in respect to all other sample names in the plot. If NULL, all samples within the SummarizedExperiment object are shown at the same luminocity. Default is NULL.

#' @param featuresToKeep Vector with feature names to keep. If NULL, all features within the SummarizedExperiment object are kept. Default is NULL. Please note that when agglomerating features with the glomby argument (see above), feature names passed to featuresToKeep must be post-agglomeration feature names. For example, if glomby="Family", featuresToKeep must be family names, such as "f__Enterobacteriaceae", etc.

#' @param subsetby String specifying the metadata variable name for subsetting samples. If passed, multiple plots will be drawn, one plot for samples within each different class contained within the variable.  If NULL, data is not subset. Default is NULL.

#' @param compareby String specifying the metadata variable name for grouping samples when the hmtype argument is set to "comparative". This will calculate p-values for each feature using the Mann-Whitney-Wilcoxon U-test when there are exactly two classes within the variable, and the log2 foldchange between the two groups will be calculated. When there are three or more classes within the variable, the p-value will be calculated using ANOVA. If there is only a single class within the variable, hmtype will default to "exploratory" and features will be ranked by variance across samples.

#' @param permanova Requires a logical value. If set to TRUE, will include in the title plot the PERMANOVA stats for groups set with compareby. Default is TRUE.

#' @param permanova_permutations Numerical value specifying the number of permutations for PERMANOVA. Default is 10000.

#' @param colourby String specifying the metadata variable name for colouring in samples. If NULL, all samples will be black. Default is NULL.

#' @param colorby Alternative US spelling for the colourby argument. Use either.

#' @param shapeby String specifying the metadata variable name for attributing shapes to samples. If NULL, all samples will be a round dot (pch = 19). Default is NULL. If there are more than 27 classes within the variable, samples will be attributed letters (A-Z, then a-z) automatically. See also use_letters_as_shapes.

#' @param use_letters_as_shapes Requires a logical value. If set to TRUE, then force sample point shapes as being letters (A-Z, then a-z) independent of how many classes there are within the variable passed to shapeby. Default is FALSE.

#' @param sizeby String specifying the metadata variable name for attributing point size to samples. If NULL, all samples are plot with the same size, specified by dotsize. Default is NULL.

#' @param dotsize Numeric value for attributing point size to samples. Default is 2.

#' @param connectby String specifying the metadata variable name for drawing a line connecting samples belonging to the same class. If NULL, samples are not connected. Default is NULL.

#' @param connection_orderby String specifying the metadata variable name for determining the order in which samples connected by the variable specified in connectby should be drawn. This is only applicable, of course, if connectby is not NULL. The classes in the metadata variable specified in connection_orderby will be sorted either numerically from low to high, if the variable contains numeric classes, or sorted alphabetically if the variable contains discrete classes.

#' @param textby String specifying the metadata variable containing classes for annotating samples with text next to each sample point. Default is NULL.

#' @param ellipseby String specifying the metadata variable containing classes for encircling with an ellipse samples belonging to each class. Default is NULL.

#' @param tsne_perplx Numerical value with perplexity for tSNE. Default is NULL.

#' @param applyfilters Optional string specifying filtration setting "combos", used as a shorthand for setting the featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff arguments in JAMS plotting functions. If NULL, none of these arguments are set if not specified. Permissible values for applyfilters are "light", "moderate" or "stringent". The actual values vary whether the SummarizedExperiment object is taxonomical (LKT) or not. For a taxonomical SummarizedExperiment object, using "light" will set featcutoff=c(50, 5), GenomeCompletenessCutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(250, 15), GenomeCompletenessCutoff=c(10, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(2000, 15), GenomeCompletenessCutoff=c(30, 10), minl2fc=2, minabscorrcoeff=0.8. For non-taxonomical (i.e. functional) SummarizedExperiment objects, using "light" will set featcutoff=c(0, 0), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(50, 15), minl2fc=2.5, minabscorrcoeff=0.8. When using applyfilters, one can still set one or more of featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff, which will then take the user set value in lieu of those set by the applyfilters shorthand. Default is light.

#' @param featcutoff Requires a numeric vector of length 2 for specifying how to filter out features by relative abundance. The first value of the vector specifies the minimum relative abundance in Parts per Million (PPM) and the second value is the percentage of samples which must have at least that relative abundance. Thus, passing c(250, 10) to featcutoff would filter out any feature which does not have at least 250 PPM (= 0.025 percent) of relative abundance in at least 10 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the featcutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param GenomeCompletenessCutoff Requires a numeric vector of length 2 for specifying how to filter out features by genome completeness. This is, of course, only applicble for taxonomic shotgun SummarizedExperiment objects. When passed to non-taxonomic shotgun SummarizedExperiment objects, GenomeCompletenessCutoff will be ignored. The first value of the vector specifies the minimum genome completeness in percentage  and the second value is the percentage of samples which must have at least that genome completeness. Thus, passing c(50, 5) to GenomeCompletenessCutoff would filter out any taxonomic feature which does not have at least 50 percent of genome completeness in at least 5 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the GenomeCompletenessCutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param plotcentroids Requires a logical value. If set to TRUE, centroids of samples within each sample group belonging to classes within the variable specified with compareby will be plot and lines will be drawn from each sample to the group centroid. Default is FALSE.

#' @param highlight_centroids Requires a logical value. If set to TRUE, when using plotcentroids, the centroids will be highlighted and marked with a slightly larger relevant group shape. Default is TRUE.

#' @param max_neighbors Numerical value specifying the maximum number of neighbors when using the "tUMAP" algorithm for dimensionality reduction (see algorithm). Default is 15.

#' @param show_centroid_distances Requires a logical value. If set to TRUE, if centroids are to be plot (see plotcentroids), will include at the bottom of the plot a matrix showing the euclidean distance between the centroids of each group. Default is FALSE.

#' @param calculate_centroid_distances_in_all_dimensions Requires a logical value. If set to TRUE, when plotcentroids and how_centroid_distances are both also set to TRUE, the euclidean

#' @param addtit Optional string with text to append to heatmap main title. Default is NULL.

#' @param assay_for_matrix String specifying the SummarizedExperiment assay to be used for the heatmap. Permissible values are "BaseCounts" or "GeneCounts". "BaseCounts" (the default) will use the basepair counts for each feature (either taxonomical or functional). These values will be normalized into relative abundance in PPM unless specified by the normalization argument (see normalization and PPM_normalize_to_bases_sequenced). When using "GeneCounts" (only available in non-taxonomical SummarizedExperiment objects) the *number of genes* annotated as each feature will be used. The heatmap will be plot with a scale of 0 to the maximum number of genes for a single feature on the heatmap. For instance, using "GeneCounts" for, let's say, an ECNumber SummarizedExperiment will plot the number of genes bearing each Enzyme Commission Number annotation within each sample. Default is "BaseCounts".

#' @param PPM_normalize_to_bases_sequenced Requires a logical value. Non-filtered JAMS feature counts tables (the BaseCounts assay within SummarizedExperiment objects) always includes unclassified taxonomical features (for taxonomical SummarizedExperiment objects) or unknown/unattributed functional features (for non-taxonomical SummarizedExperiment objects), so the relative abundance for each feature (see normalization) will be calculated in Parts per Million (PPM) by dividing the number of bases covering each feature by the sum of each sample column **previous to any filtration**. Relative abundances are thus representative of the entirety of the genomic content for taxonomical objects, whereas for non-taxonomical objects, strictly speaking, it is the abundance of each feature relative to only the coding regions present in the metagenome, even if these are annotationally unatributed. In other words, intergenic regions are not taken into account. In order to relative-abundance-normalize a **non-taxonomical** SummarizedExperiment object with the total genomic sequencing content, including non-coding regions, set PPM_normalize_to_bases_sequenced = TRUE. Default is FALSE.

#' @param grid Requires a logical value. If set to FALSE, background will be one solid color within the plot, rather than include a grid behind the plot. Default is TRUE, meaning that the background will display a grid.

#' @param forceaspectratio .

#' @param numthreads Numeric value setting the number of threads to use for any multi-threaded process within this function. The default is 1.

#' @param log2tran .

#' @param ignoreunclassified Requires a logical value. If set to TRUE, for taxonomical SummarizedExperiment objects, the feature "LKT__Unclassified" will be omitted from being shown. In the case of non-taxonomical SummarizedExperiment objects, the completely unannotated features will be omitted. For example, for an ECNumber SummarizedExperiment object, genes *without* an Enzyme Commission Number annotation (feature "EC_none") will not be shown. Statistics are, however, computed taking the completely unclassifed feature into account, so p-values will not change.

#' @param return_coordinates_matrix .

#' @param class_to_ignore String or vector specifying any classes which should lead to samples being excluded from the comparison within the variable passed to compareby. Default is N_A. This means that within any metadata variable passed to compareby containing the "N_A" string within that specific variable, the sample will be dropped from that comparison.

#' @export


plot_Ordination <- function(ExpObj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, samplesToHighlight = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, assay_for_matrix = "BaseCounts", algorithm = "tUMAP", PCA_Components = c(1, 2), distmethod = "bray", compareby = NULL, colourby = NULL, colorby = NULL, shapeby = NULL, use_letters_as_shapes = FALSE, sizeby = NULL, connectby = NULL, connection_orderby = NULL, textby = NULL, ellipseby = NULL, dotsize = 2, log2tran = TRUE, tsne_perplx = NULL, max_neighbors = 15, permanova = TRUE, plotcentroids = FALSE, highlight_centroids = TRUE, show_centroid_distances = FALSE, calculate_centroid_distances_in_all_dimensions = FALSE, addtit = NULL, cdict = NULL, grid = TRUE, forceaspectratio = NULL, numthreads = 8, return_coordinates_matrix = FALSE, permanova_permutations = 10000, class_to_ignore = "N_A", ...){

    set.seed(2138)

    #Consider orthography of the word "colour"
    if (is.null(colourby)){
        colourby <- colorby
    }

    #Hardwire PctFromCtgscutoff, as this should never be used without filtering because of huge amounts of false positives when evaluating taxonomic information from unassembled reads. The use of classifying unassembled reads is deprecated in JAMS and the default is to NOT classify unassembled reads, so this is usually not an issue.
    PctFromCtgscutoff <- c(70, 50)
    transp <- TRUE

    #Check for silly stuff
    if (!(algorithm %in% c("PCA", "PCoA", "tSNE", "tUMAP"))){
        stop(paste0("It is not possible to use ordination algorithm ", algorithm, ". Please select between PCA, PCoA, tSNE or tUMAP"))
    }

    #Define what is being compared for permanova
    if (is.null(compareby)){
        compareby <- c(colourby, shapeby, ellipseby, textby, sizeby, connectby)[1]
    }
    #Remove samples bearing categories within class_to_ignore
    valid_vars <- c(compareby, subsetby)[which(!is.na(c(compareby, subsetby)))]

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, glomby = glomby, variables_to_fix = valid_vars, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }
    pcatitbase <- paste(algorithm, "of", analysisname)

    flog.info(pcatitbase)

    if (assay_for_matrix == "GeneCounts"){
        flog.warn("Counts matrix used for ordination will represent the *number of genes* for each feature, rather than its relative abundance. For using relative abundance (default), set assay_for_matrix = \"BaseCounts\"")
        asPPM <- FALSE
        log2tran <- FALSE
        pcatitbase <- paste(pcatitbase, "Counts matrix used for ordination represents the number of genes", sep = "\n")
    }

    if (asPPM == FALSE){
        applyfilters <- FALSE
        featcutoff <- c(0, 0)
        GenomeCompletenessCutoff <- NULL
        PctFromCtgscutoff <- NULL
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff)

    if (!(is.null(subsetby))){
        subset_points <- sort(unique(colData(obj)[, which(colnames(colData(obj)) == subsetby)]))
    } else {
        subset_points <- "none"
    }

    #Create list vector to hold plots
    gvec <- list()
    plotnum <- 1

    for (sp in 1:length(subset_points)){

        if (!(is.null(subsetby))){
            sp_samplesToKeep <- rownames(colData(obj))[which(colData(obj)[ , subsetby] == subset_points[sp])]
            flog.info(paste("Plotting within", subset_points[sp]))
            subsetname <- subset_points[sp]
            pcatit <- paste(pcatitbase, "within", subset_points[sp])
        } else {
            sp_samplesToKeep <- rownames(colData(obj))
            subsetname <- "no_sub"
            pcatit <- pcatitbase
        }
        pcatit <- paste(c(pcatit, presetlist$filtermsg), collapse = "\n")
        flog.info("filtering")
        currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = sp_samplesToKeep, featuresToKeep = featuresToKeep, asPPM = asPPM, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

        currpt <- as.data.frame(colData(currobj))

        if (PPM_normalize_to_bases_sequenced == TRUE){
            pcatit <- paste(c(pcatit, "Normalized to total number of bases sequenced in sample"), collapse = "\n")
        } else {
            pcatit <- paste(c(pcatit, "Normalized to number of bases for analysis in sample"), collapse = "\n")
        }

        #Get counts matrix
        countmat <- as.matrix(assays(currobj)[[assay_for_matrix]])

        #Protect against rows with empty data
        rowsToKeep <- which(rowSums(countmat) > 0 & rownames(countmat) != "")
        countmat <- countmat[rowsToKeep, ]

        if (ignoreunclassified == TRUE){
            dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified")
            rowsToKeep <- which(!(rownames(countmat) %in% dunno))
            countmat <- countmat[rowsToKeep, ]
        }

        #log2 transform if applicable
        if (log2tran == TRUE){
            #Transform to log2 space
            countmat <- convert_matrix_log2(mat = countmat, transformation = "to_log2")
            pcatit <- paste(c(pcatit, "Matrix log2 transformed"), collapse = "\n")
        }

        n <- nrow(countmat)
        comp <- 1:2
        rowVars <- rowSds(countmat)
        countmat <- countmat[order(rowVars, decreasing = TRUE), ]
        if (transp == TRUE) {
            countmat <- t(countmat)
        }

        sampsizemsg <- paste0("Number of samples = ", ncol(countmat), "; Number of features = ", nrow(countmat))
        flog.info("getting distance")

        if (!is.null(samplesToHighlight)){
            samplesToHighlight <- samplesToHighlight[samplesToHighlight %in% rownames(countmat)]
            currpt_stat <- currpt[(rownames(currpt) %in% samplesToHighlight), ]
            d <- vegdist(countmat[(rownames(countmat) %in% samplesToHighlight), ], method = distmethod, na.rm = TRUE)
            cats <- currpt_stat[ , compareby]
        } else {
            currpt_stat <- currpt
            d <- vegdist(countmat, method = distmethod, na.rm = TRUE)
            cats <- currpt_stat[, compareby]
        }

        #switch off permanova if there is only a single comparison
        if (length(cats) < 2){
            permanova <- FALSE
        }

        flog.info("getting permanova")
        if (permanova == TRUE){
            if (!(is.numeric(cats))){
                permanovaout <- vegan::adonis2(formula = as.formula(paste("d ~ ", compareby)), data = currpt_stat, permutations = permanova_permutations, parallel = numthreads)
                permanovap <- permanovaout$`Pr(>F)`[1]
                permanovaF <- permanovaout$`F`[1]
            } else {
                flog.info("Impossible to get permanova because compareby is continuous.")
                permanovap <- NULL
                permanovaF <- NULL
            }
        } else {
            permanovap <- NULL
            permanovaF <- NULL
        }

        flog.info("ordinating")
        if (algorithm == "tSNE"){
            #tSNE algorithm
            if (is.null(tsne_perplx)){
                tsne_perplx <- round(nrow(currpt) * 0.3, 0)
            }

            tsne_out <- Rtsne(countmat, dims = 3, initial_dims = 500, perplexity = tsne_perplx, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000)
            dford <- as.data.frame(tsne_out$Y)
            rownames(dford) <- rownames(currpt)
            colnames(dford)[1:2] <- c("PC1", "PC2")
            axisprefix <- "PC"
            dford_full <- dford
            xl <- "tSNE 1"
            yl <- "tSNE 2"

        } else if (algorithm == "tUMAP"){
            set.seed(2138)
            n_neighbors <- min((nrow(countmat) - 1), max_neighbors)

            tumap_out <- tumap(countmat, n_components = 2, n_neighbors = n_neighbors, verbose = FALSE, n_threads = numthreads, init = "spca")
            dford <- as.data.frame(tumap_out)
            rownames(dford) <- rownames(currpt)
            colnames(dford) <- c("PC1", "PC2")
            axisprefix <- "PC"
            dford_full <- dford
            xl <- "tUMAP 1"
            yl <- "tUMAP 2"

        } else {
            if (algorithm == "PCA"){
                pcaRes <- prcomp(d)
                ord <- pcaRes$x
                vars <- pcaRes$sdev^2
                vars <- round(vars/sum(vars), 5) * 100
                axisprefix <- "PC"
            } else {
                #Default to PCoA
                pcoaRes <- ape::pcoa(d, correction = "none")
                ord <- pcoaRes$vectors
                vars <- pcoaRes$values$Broken_stick
                vars <- round(vars, 5) * 100
                axisprefix <- "Axis."
            }

            #Make a data frame with how variance is explained by which components
            vardf <- data.frame(Component = colnames(ord), Variance = vars, Cumulative_variance = cumsum(vars))
            vardf$Component <- factor(vardf$Component, levels = vardf$Component)

            varplot <- ggplot(vardf, aes(x = Component)) + geom_bar(aes(y = Variance), fill = 'blue', stat = "identity") + geom_point(aes(y = Cumulative_variance), colour = "black", pch = 16, size = 1) + geom_path(aes(y = Cumulative_variance, group = 1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.6, size = rel(0.5))) + labs(title = "Variance explained by each PCoA component", subtitle = pcatit, x = 'Component', y = 'Variance')

            #Take into account different PCs if applicable
            if (any(PCA_Components != c(1, 2))){
                comp <- PCA_Components[1:2]
                flog.warn(paste("WARNING: plotting PCA with components", paste0(paste0("PC", comp), collapse = " and ")))
            }

            xl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[1]], vars[comp[1]])
            yl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[2]], vars[comp[2]])
            #zl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[3]], vars[comp[3]])
            dford_full <- as.data.frame(ord)
            dford <- as.data.frame(ord[, comp])
        }

        #Add colour, size, shape, text, ellipse and connect
        dford$Comparison <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == compareby)]
        if (!is.null(colourby)){
            dford$Colours <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == colourby)]
        }
        if (!is.null(shapeby)){
            dford$Shape <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == shapeby)]
        }
        if (!is.null(ellipseby)){
            dford$Ellipse <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == ellipseby)]
        }
        if (!is.null(textby)){
            dford$Text <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == textby)]
        }
        if (!is.null(sizeby)){
            dford$Size <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == sizeby)]
        }
        if (!is.null(connectby)){
            dford$Connect <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == connectby)]
            #If reorder dataframe if required, so geom_path will connect in the right order
            if (!is.null(connection_orderby)){
                dford$ConnectOrder <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == connection_orderby)]
                sampleorder <- NULL
                for (grp in unique(dford$Connect)){
                    currdford <- subset(dford, Connect == grp)[ , c("Connect", "ConnectOrder")]
                    if (can_be_made_numeric(currdford$ConnectOrder)){
                        sampleorder <-c(sampleorder, rownames(currdford[order(as.numeric(currdford$ConnectOrder)), ]))
                    } else {
                        sampleorder <-c(sampleorder, rownames(currdford[order(currdford$ConnectOrder), ]))
                    }
                }
                dford <- dford[sampleorder, ]
            }
        }

        if (!is.null(samplesToHighlight)){
            dford$Alpha <- 0.05
            dford[samplesToHighlight, "Alpha"] <- 1
        } else {
            dford$Alpha <- 1
        }

        #flog.info("getting centroids")
        if (!(is.numeric(cats))){
            if (!is.null(samplesToHighlight)){
                centroids <- aggregate(.~Comparison, data = dford[samplesToHighlight, c(paste0(axisprefix, comp), "Comparison")], FUN = mean)
                colnames(centroids)[c(2, 3)] <- paste0("mean", colnames(dford)[1:2])
                rownames(centroids) <- centroids[ , "Comparison"]
                centroiddf <- left_join(dford[samplesToHighlight, ], centroids, by = "Comparison")
            } else {
                centroids <- aggregate(.~Comparison, data = dford[ , c(paste0(axisprefix, comp), "Comparison")], FUN = mean)
                colnames(centroids)[c(2, 3)] <- paste0("mean", colnames(dford)[1:2])
                centroiddf <- left_join(dford, centroids, by = "Comparison")
                rownames(centroids) <- centroids[ , "Comparison"]
            }
            centroidmeandf <- centroiddf

            centroidmeandf <- centroidmeandf[ , 3:ncol(centroidmeandf)]
            centroidmeandf <- centroidmeandf[!duplicated(centroidmeandf[ , paste0("mean", colnames(dford)[1])]), ]
            rownames(centroidmeandf) <- centroidmeandf[ , "Comparison"]

        }

        flog.info("plotting")
        if (!is.null(samplesToHighlight)){
            aesthetic <- aes(x = get(colnames(dford)[1]), y = get(colnames(dford)[2]), alpha = Alpha)
        } else {
            aesthetic <- aes(x = get(colnames(dford)[1]), y = get(colnames(dford)[2]))
        }
        p <- ggplot(dford, aesthetic)

        if (return_coordinates_matrix){
            gvec[[plotnum]] <- dford
            plotnum <- plotnum + 1
        }

        if (!is.null(colourby)) {
            p <- p + aes(col = Colours)
        }

        if (!(is.null(shapeby))){
            p <- p + aes(shape = Shape)
            p <- add_shape_to_plot_safely(p = p, shapevec = dford$Shape, shapeby = shapeby, use_letters_as_shapes = use_letters_as_shapes, cdict = cdict)
        }

        if (!(is.null(sizeby))){
            p <- p + aes(size = Size)
            numsizes <- length(unique(dford$Size))
            p <- p + scale_shape_manual(values = dotsize:(dotsize + numsizes))
        }

        if (is.numeric(dford$Colours)){
            #Check if there is enough variance in the continuous data to plot a gradient
            if ((max(dford$Colours) - min(dford$Colours)) > 0){
                p <- p + scale_color_gradient(low = "blue", high = "red")
            }
        } else {
            #if there is a colour dictionary, then use that
            if (!(is.null(cdict))){
                ct <- cdict[[colourby]]
                groupcols <- setNames(as.character(ct$Colour), as.character(ct$Name))
                p <- p + scale_color_manual(values = groupcols)
            } else {
                #Use colour table if available
                if ("ctable" %in% names(metadata(currobj))){
                    discretenames <- sort(unique(dford$Colours))
                    colourshave <- discretenames[discretenames %in% rownames(metadata(currobj)$ctable)]
                    cores <- as.vector(rainbow(length(discretenames)))
                    names(cores) <- discretenames
                    cores[colourshave] <- metadata(currobj)$ctable[colourshave, "Hex"]
                    p <- p + scale_color_manual(values = cores)
                }
            }
        }

        if (!is.null(connectby)){
            p <- p + geom_path(color = "black", size = (dotsize * 0.25))
        }

        if (!is.null(textby)){
            require(ggrepel)
            p <- p + geom_text_repel(aes(label = Text), size = (dotsize * 0.75), min.segment.length = 0, segment.size = 0.05, segment.color = "black", segment.curvature = -0.05, segment.ncp = 3, segment.inflect = FALSE, max.overlaps = 50)
        }

        p <- p + geom_point(size = dotsize) + labs(x = xl, y = yl)

        if (any(c((!is.null(ellipseby)), (ellipseby != FALSE)))){
            p <- p + stat_ellipse(aes(group = Ellipse), type = "norm")
        }

        if (!(is.null(addtit))){
            pcatit <- paste(c(pcatit, addtit), collapse = "\n")
        }

        if (!is.null(permanovap)) {
            pcatit <- paste(c(pcatit, paste0("PERMANOVA p <", round(permanovap, 4), "; F-Model = ", round(permanovaF, 4))), collapse = "\n")
        }

        if(any(c((!is.null(permanovap)), (algorithm == "PCA")))){
            distmessage <- paste0("Dissimilarity index = ", distmethod)
            pcatit <- paste(c(pcatit, distmessage), collapse = "\n")
        }

        if (!(is.null(shapeby))){
            p <- p + labs(shape = shapeby)
        }

        if (!(is.null(sizeby))){
            p <- p + labs(size = sizeby)
        }

        p <- p + theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))

        if (!(is.null(forceaspectratio))){
            p <- p + theme(aspect.ratio = (1 / forceaspectratio))
        }

        if (all(c((!(is.numeric(cats))), plotcentroids))) {

            #Account for not using colours. I am sorry this is messy and unelegant, will tidy up later. Using safest method.
            if (!is.null(colourby)){
                p <- p + geom_segment(aes(x = get(colnames(centroiddf)[ncol(centroiddf) - 1]), y = get(colnames(centroiddf)[ncol(centroiddf)]), xend = get(colnames(centroiddf)[1]), yend = get(colnames(centroiddf)[2]), colour = Colours), data = centroiddf)
            } else {
                p <- p + geom_segment(aes(x = get(colnames(centroiddf)[ncol(centroiddf) - 1]), y = get(colnames(centroiddf)[ncol(centroiddf)]), xend = get(colnames(centroiddf)[1]), yend = get(colnames(centroiddf)[2])), data = centroiddf)
            }

            if (highlight_centroids){
                p <- p + geom_point(aes(x = get(colnames(centroidmeandf)[ncol(centroidmeandf) - 1]), y = get(colnames(centroidmeandf)[ncol(centroidmeandf)])), colour = "black", data = centroiddf, size = (dotsize * 3))
                if (!is.null(colourby)){
                    p <- p + geom_point(aes(x = get(colnames(centroidmeandf)[ncol(centroidmeandf) - 1]), y = get(colnames(centroidmeandf)[ncol(centroidmeandf)]), colour = Colours), data = centroiddf, size = (dotsize * 1.5))
                } else {
                    p <- p + geom_point(aes(x = get(colnames(centroidmeandf)[ncol(centroidmeandf) - 1]), y = get(colnames(centroidmeandf)[ncol(centroidmeandf)])), data = centroiddf, size = (dotsize * 1.5))
                }
            }
        }

        #pcatit <- paste(pcatit, sampsizemsg, sep = "\n")
        p <- p + ggtitle(pcatit)
        p <- p + labs(colour = colourby)
        p <- p + guides(alpha = "none")

        if (grid == FALSE ){
            p <- p + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
        }

        if (all(c((!(is.numeric(cats))), plotcentroids, show_centroid_distances))){

            if (!calculate_centroid_distances_in_all_dimensions){
                centroiddist <- as.matrix(dist(centroidmeandf[ , paste0("mean", colnames(dford)[1:2])], method = "euclidean"))

            } else {
                flog.info(paste("Calculating the euclidean distances between centroids in", ncol(dford_full), "dimesnions"))
                dford_full$Comparison <- currpt[match(rownames(dford_full), rownames(currpt)), which(colnames(currpt) == compareby)]

                if (!is.null(samplesToHighlight)){
                    centroids_full <- aggregate(.~Comparison, data = dford_full[samplesToHighlight, ], FUN = mean)
                } else {
                    centroids_full <- aggregate(.~Comparison, data = dford_full, FUN = mean)
                }
                colnames(centroids_full)[c(2:ncol(centroids_full))] <- paste0("mean", colnames(centroids_full)[c(2:ncol(centroids_full))])
                rownames(centroids_full) <- centroids_full[ , "Comparison"]
                centroiddist <- as.matrix(dist(centroids_full[,2:ncol(centroids_full)], method = "euclidean"))
            }

            centroiddist <- round(centroiddist, 3)

            centroiddistshow <- ggtexttable(centroiddist, theme = ttheme(base_style = "classic", base_size = 8))

            p <- ggarrange(p, centroiddistshow, ncol = 1, nrow = 2, heights = c(3, 1), labels = list("", "Euclidean distance between centroids"), font.label = list(size = 10, face = "italic"), vjust = 1, hjust = -1)
        }

        gvec[[plotnum]] <- p
        plotnum <- plotnum + 1

    }

    return(gvec)
}
