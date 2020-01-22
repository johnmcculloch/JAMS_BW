#' plot_Ordination(ExpObj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, mgSeqnorm = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, algorithm = "PCA", colourby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, dotsize = 2, dotborder = NULL, log2tran = FALSE, transp = TRUE, perplx = NULL, permanova = FALSE, ellipse = FALSE, plotcentroids = FALSE, addtit = NULL, plot3D = FALSE, theta = 130, phi = 60, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 1, class_to_ignore = "N_A", ...)
#'
#' Creates ordination plots based on PCA, tSNE or tUMAP
#' @export

plot_Ordination <- function(ExpObj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, mgSeqnorm = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, algorithm = "PCA", colourby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, dotsize = 2, dotborder = NULL, log2tran = TRUE, transp = TRUE, perplx = NULL, permanova = TRUE, ellipse = FALSE, plotcentroids = FALSE, addtit = NULL, plot3D = FALSE, theta = 130, phi = 60, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 1, class_to_ignore = "N_A", ...){

    #Remove samples bearing categories within class_to_ignore
    valid_vars <- c(colourby, shapeby, sizeby, subsetby)[which(!is.na(c(colourby, shapeby, sizeby, subsetby)))]

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, glomby = glomby, variables_to_fix = valid_vars, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff)

    if (!(is.null(subsetby))){
        subset_points <- sort(unique(colData(obj)[, which(colnames(colData(obj)) == subsetby)]))
    } else {
        subset_points <- "none"
    }

    #Create list vector to hold plots
    gvec <- NULL
    gvec <- vector("list", length = length(subset_points))
    pcatitbase <- paste(algorithm, "of", analysis)

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
            d <- vegdist(countmat, method = "jaccard", na.rm = TRUE)
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

            set.seed(4140)
            tsne_out <- Rtsne(countmat, dims = 3, initial_dims = 500, perplexity = perplx, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000)
            dford <- as.data.frame(tsne_out$Y)
            rownames(dford) <- rownames(currpt)
            colnames(dford)[1:3] <- c("PC1", "PC2", "PC3")
            xl <- "tSNE 1"
            yl <- "tSNE 2"
            zl <- "tSNE 3"

        } else if (algorithm == "tUMAP"){

            set.seed(4140)
            if (nrow(countmat) < 20){
                n_neighbors <- (nrow(countmat) - 1)
            } else {
                n_neighbors <- 20
            }
            tumap_out <- tumap(countmat, n_components = 2, n_neighbors = n_neighbors, verbose = FALSE, n_threads = threads)
            dford <- as.data.frame(tumap_out)
            rownames(dford) <- rownames(currpt)
            colnames(dford)[1:2] <- c("PC1", "PC2")
            xl <- "tUMAP 1"
            yl <- "tUMAP 2"
            #zl <- "tSNE 3"

        } else {
            #Not tSNE or tUMAP, so use PCA
            distfun <- stats::dist
            if (permanova == FALSE){
                #get distance if missing
                #d <- distfun(mat, method = "euclidian")
                d <- vegdist(countmat, method = "jaccard", na.rm = TRUE)
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

        #centroids <- aggregate(cbind(PC1, PC2) ~ Colours, dford, mean)
        #colnames(centroids)[c(2,3)] <- c("meanx", "meany")
        #dford <- left_join(dford, centroids)

        if (plot3D == TRUE) {
            aesthetic <- aes(x=PC1, y=PC2, z=PC3)
        } else {
            aesthetic <- aes(x=PC1, y=PC2)
        }
        p <- ggplot(dford, aesthetic)

        if (!is.null(colourby)) {
            p <- p + aes(col = Colours)
        }

        if (!(is.null(shapeby))){
            p <- p + aes(shape = Shape)
            numshapes <- length(unique(dford$Shape))
            p <- p + scale_shape_manual(values = 15:(numshapes + 15))
        }

        if (!(is.null(sizeby))){
            p <- p + aes(size = Size)
            numsizes <- length(unique(dford$Size))
            p <- p + scale_shape_manual(values = dotsize:(dotsize + numsizes))
        }

        if (is.numeric(dford$Colours)){
            #Check if there is enough variance in the continuous data to plot a gradient
            if ((max(dford$Colours) - min(dford$Colours)) > 0){
                p <- p + scale_color_gradient(low="blue", high="red")
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

        if (ellipse != FALSE) {
            if (algorithm == "PCA"){
                if (ellipse == "auto"){
                    if (all(c(!is.null(permanovap), (permanovap < 0.05)))){
                        p <- p + stat_ellipse(type = "norm")
                    }
                } else if (ellipse == TRUE) {
                    p <- p + stat_ellipse(type = "norm")
                }
            } else {
                if (ellipse == TRUE) {
                    p <- p + stat_ellipse(type = "norm")
                }
            }
        }

        #if (plotcentroids == TRUE){
        #    p <- p + geom_point(data = centroids, size = 5) + geom_segment(aes(x=dford$meanx, y=dford$meany, xend=dford$PC1, yend=dford$PC2))
        #}

        if (!(is.null(addtit))){
            pcatit <- paste(c(pcatit, addtit), collapse = "\n")
        }

        if (!is.null(permanovap)) {
            pcatit <- paste(c(pcatit, (paste("PERMANOVA p <", permanovap))), collapse = "\n")
        }

        p <- p + ggtitle(pcatit)
        p <- p + labs(colour = colourby)

        if (!(is.null(shapeby))){
            p <- p + labs(shape = shapeby)
        }

        if (!(is.null(sizeby))){
            p <- p + labs(size = sizeby)
        }

        p <- p + theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))
        if (plot3D == TRUE) {
            p <- p + axes_3D(theta = theta, phi = phi) + stat_3D()
            p <- p + labs_3D(labs = c(xl,yl,zl), hjust = c(-0.25, 1,0.25)) + labs(x = "", y = "")
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
        } else {
            p <- p + geom_point(size = dotsize) + labs(x = xl, y = yl)
            if (!(is.null(forceaspectratio))){
                p <- p + theme(aspect.ratio = (1 / forceaspectratio))
            }
        }

        if (grid == FALSE ){
            p <- p + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
        }
        gvec[[sp]] <- p
    }

    gvec <- gvec[sapply(gvec, function(x){ !(is.null(x)) } )]

    return(gvec)
}
