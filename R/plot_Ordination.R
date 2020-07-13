#' plot_Ordination(ExpObj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, algorithm = "PCA", distmethod = "jaccard", colourby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, dotsize = 2, dotborder = NULL, log2tran = TRUE,  transp = TRUE, perplx = NULL, permanova = TRUE, ellipse = FALSE, plotcentroids = FALSE, show_centroid_distances = FALSE, addtit = NULL, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 1, class_to_ignore = "N_A", ...)
#'
#' Creates ordination plots based on PCA, tSNE or tUMAP
#' @export

plot_Ordination <- function(ExpObj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, algorithm = "PCA", distmethod = "jaccard", colourby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, dotsize = 2, dotborder = NULL, log2tran = TRUE,  transp = TRUE, perplx = NULL, max_neighbors = 15, permanova = TRUE, ellipse = FALSE, plotcentroids = FALSE, highlight_centroids = TRUE, show_centroid_distances = FALSE, addtit = NULL, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 1, class_to_ignore = "N_A", ...){

    set.seed(4140)

    #Remove samples bearing categories within class_to_ignore
    valid_vars <- c(colourby, shapeby, sizeby, subsetby)[which(!is.na(c(colourby, shapeby, sizeby, subsetby)))]

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
    gvec <- NULL
    gvec <- vector("list", length = length(subset_points))
    pcatitbase <- paste(algorithm, "of", analysisname)

    for (sp in 1:length(subset_points)){

        if (!(is.null(subsetby))){
            samplesToKeep <- rownames(colData(obj))[which(colData(obj)[ , subsetby] == subset_points[sp])]
            flog.info(paste("Plotting within", subset_points[sp]))
            subsetname <- subset_points[sp]
            pcatit <- paste(pcatitbase, "within", subset_points[sp])
        } else {
            samplesToKeep <- rownames(colData(obj))
            subsetname <- "no_sub"
            pcatit <- pcatitbase
        }
        pcatit <- paste(c(pcatit, presetlist$filtermsg), collapse = "\n")

        currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, asPPM = TRUE, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

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

        if (permanova == TRUE){
            d <- vegdist(countmat, method = distmethod, na.rm = TRUE)
            cats <- currpt[, colourby]
            if (!(is.numeric(cats))){
                permanovap <- vegan::adonis(as.formula(paste("d ~ ", colourby)), data = currpt)$aov.tab$`Pr(>F)`[1]
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
            if (permanova == FALSE){
                #get distance if missing
                #d <- distfun(mat, method = "euclidian")
                d <- vegdist(countmat, method = distmethod, na.rm = TRUE)
            }
            pcaRes <- prcomp(d)
            ord <- pcaRes$x
            vars <- pcaRes$sdev^2
            vars <- round(vars/sum(vars), 5) * 100

            xl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[1]], vars[comp[1]])
            yl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[2]], vars[comp[2]])
            zl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[3]], vars[comp[3]])
            dford <- as.data.frame(ord[, comp])
        }
        #Add colour, size, shape
        dford$Colours <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == colourby)]
        dford$Size <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == sizeby)]
        dford$Shape <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == shapeby)]
        dford$Pair <- currpt[match(rownames(dford), rownames(currpt)), which(colnames(currpt) == pairby)]
        centroids <- aggregate(cbind(PC1, PC2) ~ Colours, dford, mean)
        colnames(centroids)[c(2, 3)] <- c("meanPC1", "meanPC2")
        centroiddf <- left_join(dford, centroids, by = "Colours")
        rownames(centroids) <- centroids$Colours

        aesthetic <- aes(x = PC1, y = PC2)
        p <- ggplot(dford, aesthetic)

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

        if (ellipse == "auto"){
            if (all(c(!is.null(permanovap), (permanovap < 0.05)))){
                p <- p + stat_ellipse(show.legend = TRUE, type = "norm")
            }
        } else if (ellipse == TRUE) {
            p <- p + stat_ellipse(show.legend = TRUE, type = "norm")
        }

        if (!(is.null(addtit))){
            pcatit <- paste(c(pcatit, addtit), collapse = "\n")
        }

        if (!is.null(permanovap)) {
            pcatit <- paste(c(pcatit, (paste("PERMANOVA p <", permanovap))), collapse = "\n")
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
        p <- p + geom_point(size = dotsize) + labs(x = xl, y = yl)
        if (!(is.null(forceaspectratio))){
            p <- p + theme(aspect.ratio = (1 / forceaspectratio))
        }

        if (plotcentroids){
            p <- p + geom_segment(aes(x = meanPC1, y = meanPC2, xend = PC1, yend = PC2, colour = Colours), data = centroiddf)
            if (highlight_centroids){
                p <- p + geom_point(aes(x = meanPC1, y = meanPC2), colour = "black", data = centroids, size = (dotsize * 4)) + geom_point(aes(x = meanPC1, y = meanPC2, colour = Colours), data = centroids, size = (dotsize * 2))
            }
        }

        p <- p + ggtitle(pcatit)
        p <- p + labs(colour = colourby)

        if (grid == FALSE ){
            p <- p + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
        }

        if (all(c(plotcentroids, show_centroid_distances))){
            centroiddist <- as.matrix(dist(centroids[ , c("meanPC1", "meanPC2")], method = "euclidean"))
            centroiddist <- round(centroiddist, 3)
            centroiddistshow <- ggtexttable(centroiddist, theme = ttheme(base_style = "classic", base_size = 8))
            p <- ggarrange(p, centroiddistshow, ncol = 1, nrow = 2, heights = c(3, 1), labels = list("", "Euclidean distance between centroids"), font.label = list(size = 10, face = "italic"), vjust = 1, hjust = -1)
        }

        gvec[[sp]] <- p
    }

    gvec <- gvec[sapply(gvec, function(x){ !(is.null(x)) } )]

    return(gvec)
}
