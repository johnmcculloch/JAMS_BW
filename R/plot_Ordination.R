#' plot_Ordination(ExpObj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, samplesToHighlight = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, algorithm = "PCA", distmethod = "jaccard", compareby = NULL, colourby = NULL, colorby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, textby = NULL, ellipseby = NULL, dotsize = 2, dotborder = NULL, log2tran = TRUE,  transp = TRUE, perplx = NULL, max_neighbors = 15, permanova = TRUE, plotcentroids = FALSE, highlight_centroids = TRUE, show_centroid_distances = FALSE, addtit = NULL, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 8, return_coordinates_matrix = FALSE, permanova_permutations = 10000, class_to_ignore = "N_A", ...)
#'
#' Creates ordination plots based on PCA, tSNE or tUMAP
#' @export

plot_Ordination <- function(ExpObj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, samplesToHighlight = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, algorithm = "PCA", distmethod = "jaccard", compareby = NULL, colourby = NULL, colorby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, textby = NULL, ellipseby = NULL, dotsize = 2, dotborder = NULL, log2tran = TRUE,  transp = TRUE, perplx = NULL, max_neighbors = 15, permanova = TRUE, plotcentroids = FALSE, highlight_centroids = TRUE, show_centroid_distances = FALSE, addtit = NULL, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 8, return_coordinates_matrix = FALSE, permanova_permutations = 10000, class_to_ignore = "N_A", ...){

    set.seed(2138)
    #Consider orthography of the word "colour"
    if (is.null(colourby)){
        colourby <- colorby
    }

    #Define what is being compared for permanova
    if (is.null(compareby)){
        compareby <- c(colourby, shapeby, ellipseby, textby, sizeby, pairby)[1]
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

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff)

    if (!(is.null(subsetby))){
        subset_points <- sort(unique(colData(obj)[, which(colnames(colData(obj)) == subsetby)]))
    } else {
        subset_points <- "none"
    }

    #Create list vector to hold plots
    gvec <- list()
    plotnum <- 1
    pcatitbase <- paste(algorithm, "of", analysisname)

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

        currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = sp_samplesToKeep, featuresToKeep = featuresToKeep, asPPM = asPPM, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

        currpt <- as.data.frame(colData(currobj))

        if (PPM_normalize_to_bases_sequenced == TRUE){
            pcatit <- paste(c(pcatit, "Normalized to total number of bases sequenced in sample"), collapse = "\n")
        } else {
            pcatit <- paste(c(pcatit, "Normalized to number of bases for analysis in sample"), collapse = "\n")
        }

        #Get counts matrix
        countmat <- as.matrix(assays(currobj)$BaseCounts)

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
        comp <- 1:3
        rowVars <- rowSds(countmat)
        countmat <- countmat[order(rowVars, decreasing = TRUE), ]
        if (transp == TRUE) {
            countmat <- t(countmat)
        }

        sampsizemsg <- paste0("Number of samples = ", ncol(countmat), "; Number of features = ", nrow(countmat))

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

        if (permanova == TRUE){
            if (!(is.numeric(cats))){
                permanovap <- vegan::adonis(as.formula(paste("d ~ ", compareby)), data = currpt_stat, permutations = permanova_permutations)$aov.tab$`Pr(>F)`[1]
            } else {
                flog.info("Impossible to get permanova because colourby is continuous")
                permanovap <- NULL
            }
        } else {
            permanovap <- NULL
        }

        if (algorithm == "tSNE"){
            #tSNE algorithm
            if (is.null(perplx)){
                perplx <- round(nrow(currpt) * 0.3, 0)
            }

            tsne_out <- Rtsne(countmat, dims = 3, initial_dims = 500, perplexity = perplx, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000)
            dford <- as.data.frame(tsne_out$Y)
            rownames(dford) <- rownames(currpt)
            colnames(dford)[1:3] <- c("PC1", "PC2", "PC3")
            xl <- "tSNE 1"
            yl <- "tSNE 2"
            zl <- "tSNE 3"

        } else if (algorithm == "tUMAP"){
            set.seed(2138)
            n_neighbors <- min((nrow(countmat) - 1), max_neighbors)

            tumap_out <- tumap(countmat, n_components = 2, n_neighbors = n_neighbors, verbose = FALSE, n_threads = threads, init = "spca")
            dford <- as.data.frame(tumap_out)
            rownames(dford) <- rownames(currpt)
            colnames(dford) <- c("PC1", "PC2")
            xl <- "tUMAP 1"
            yl <- "tUMAP 2"

        } else {
            #Not tSNE or tUMAP, so use PCA
            #distfun <- stats::dist
            pcaRes <- prcomp(d)
            ord <- pcaRes$x
            vars <- pcaRes$sdev^2
            vars <- round(vars/sum(vars), 5) * 100

            xl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[1]], vars[comp[1]])
            yl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[2]], vars[comp[2]])
            zl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[3]], vars[comp[3]])
            dford <- as.data.frame(ord[, comp])
        }

        #Add colour, size, shape, text, ellipse and pair
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
        if (!is.null(pairby)){
            dford$Pair <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == pairby)]
        }

        if (!is.null(samplesToHighlight)){
            dford$Alpha <- 0.05
            dford[samplesToHighlight, "Alpha"] <- 1
        } else {
            dford$Alpha <- 1
        }

        if (!(is.numeric(cats))){
            if (!is.null(samplesToHighlight)){
                centroids <- aggregate(cbind(PC1, PC2) ~ Comparison, dford[samplesToHighlight, ], mean)
                colnames(centroids)[c(2, 3)] <- c("meanPC1", "meanPC2")
                rownames(centroids) <- centroids[ , "Comparison"]
                centroiddf <- left_join(dford[samplesToHighlight, ], centroids, by = "Comparison")
            } else {
                centroids <- aggregate(cbind(PC1, PC2) ~ Comparison, dford, mean)
                colnames(centroids)[c(2, 3)] <- c("meanPC1", "meanPC2")
                centroiddf <- left_join(dford, centroids, by = "Comparison")
                rownames(centroids) <- centroids[ , "Comparison"]
            }
            centroidmeandf <- centroiddf
            centroidmeandf$PC1 <- NULL
            centroidmeandf$PC2 <- NULL
            centroidmeandf <- centroidmeandf[!duplicated(centroidmeandf$meanPC1), ]
        }

        if (!is.null(samplesToHighlight)){
            aesthetic <- aes(x = PC1, y = PC2, alpha = Alpha)
        } else {
            aesthetic <- aes(x = PC1, y = PC2)
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
            p <- add_shape_to_plot_safely(p = p, shapevec = dford$Shape, shapeby = shapeby, cdict = cdict)
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
            }
        }

        if (!is.null(pairby)){
            p <- p + aes(group = Pair) + geom_line()
        }

        if (!is.null(textby)){
            require(ggrepel)
            #p <- p + geom_text_repel(aes(label = Text), size = dotsize, segment.color = 'transparent')
            p <- p + geom_text_repel(aes(label = Text), size = (dotsize * 0.9), segment.size = 0.2)
        }

        p <- p + geom_point(size = dotsize) + labs(x = xl, y = yl)

        if (all(c((!is.null(ellipseby)), (ellipseby == FALSE)))){
            #p <- p + stat_ellipse(show.legend = TRUE, type = "norm")
            #if (!is.null(colourby)){
                #p <- p + stat_ellipse(aes(color = Colours, group = Ellipse), type = "norm")
            #} else {
                p <- p + stat_ellipse(aes(group = Ellipse), type = "norm")
            #}
        }

        if (!(is.null(addtit))){
            pcatit <- paste(c(pcatit, addtit), collapse = "\n")
        }

        if (!is.null(permanovap)) {
            pcatit <- paste(c(pcatit, (paste("PERMANOVA p <", round(permanovap, 4)))), collapse = "\n")
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
            p <- p + geom_segment(aes(x = meanPC1, y = meanPC2, xend = PC1, yend = PC2, colour = Colours), data = centroiddf)
            if (highlight_centroids){
                p <- p + geom_point(aes(x = meanPC1, y = meanPC2), colour = "black", data = centroidmeandf, size = (dotsize * 4)) + geom_point(aes(x = meanPC1, y = meanPC2, colour = Colours), data = centroidmeandf, size = (dotsize * 2))
            }
        }

        #pcatit <- paste(pcatit, sampsizemsg, sep = "\n")

        p <- p + ggtitle(pcatit)
        p <- p + labs(colour = colourby)
        p <- p + guides(alpha = FALSE)

        if (grid == FALSE ){
            p <- p + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
        }

        if (all(c((!(is.numeric(cats))), plotcentroids, show_centroid_distances))){
            centroiddist <- as.matrix(dist(centroidmeandf[ , c("meanPC1", "meanPC2")], method = "euclidean"))
            centroiddist <- round(centroiddist, 3)
            centroiddistshow <- ggtexttable(centroiddist, theme = ttheme(base_style = "classic", base_size = 8))
            p <- ggarrange(p, centroiddistshow, ncol = 1, nrow = 2, heights = c(3, 1), labels = list("", "Euclidean distance between centroids"), font.label = list(size = 10, face = "italic"), vjust = 1, hjust = -1)
        }

        gvec[[plotnum]] <- p
        plotnum <- plotnum + 1
    }

    return(gvec)
}
