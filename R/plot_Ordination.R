#' plot_Ordination(mgseqobj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, mgSeqnorm = FALSE, featmaxatleastPPM = 0, featcutoff = c(0, 0), applyfilters = NULL,  algorithm = "PCA", colourby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, dotsize = 2, dotborder = NULL, log2tran = TRUE, transp = TRUE, perplx = NULL, permanova = FALSE, ellipse = FALSE, plotcentroids = FALSE, plottit = NULL, plot3D = FALSE, theta = 130, phi = 60, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 8, class_to_ignore = "N_A", ...)
#'
#' Creates ordination plots based on PCA, tSNE or tUMAP
#' @export

plot_Ordination <- function(mgseqobj = NULL, glomby = NULL, subsetby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, mgSeqnorm = FALSE, featmaxatleastPPM = 0, featcutoff = c(0, 0), applyfilters = NULL,  algorithm = "PCA", colourby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, dotsize = 2, dotborder = NULL, log2tran = TRUE, transp = TRUE, perplx = NULL, permanova = FALSE, ellipse = FALSE, plotcentroids = FALSE, plottit = NULL, plot3D = FALSE, theta = 130, phi = 60, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 8, class_to_ignore = "N_A", ...){

    #Get appropriate object to work with
    obj <- mgseqobj

    #Exclude samples and features if specified
    if (!(is.null(samplesToKeep))){
        obj <- obj[, samplesToKeep]
    }

    if(!(is.null(featuresToKeep))){
        obj <- obj[featuresToKeep, ]
    }

    #Aggregate if required
    if (!(is.null(glomby))){
        obj <- agglomerate_features(mgseqobj = obj, glomby = glomby)
    }

    #Remove samples bearing categories within class_to_ignore
    valid_vars <- c(colourby, shapeby, sizeby, subsetby)[which(!is.na(c(colourby, shapeby, sizeby, subsetby)))]
    obj <- filter_sample_by_class_to_ignore(mgseqobj = obj, variables = valid_vars, class_to_ignore = class_to_ignore)

    #Define analysis type
    analysis <- attr(obj, "analysis")

    if (!is.null(applyfilters)){
        if (applyfilters == "stringent"){
            if (analysis == "LKT"){
                featcutoff <- c(2000, 15)
                genomecompleteness <- 0.3
            } else {
                featcutoff <- c(50, 15)
                genomecompleteness <- NULL
            }
        } else if (applyfilters == "moderate"){
            if (analysis == "LKT"){
                featcutoff <- c(500, 10)
                genomecompleteness <- 0.1
            } else {
                featcutoff <- c(10, 5)
                genomecompleteness <- NULL
            }
        }
    } else {
        featcutoff <- c(0, 0)
        genomecompleteness <- NULL
    }

    if (!(is.null(subsetby))){
        subset_points <- sort(unique((pData(obj)[, which(colnames(pData(obj)) == subsetby)])))
    } else {
        subset_points <- "none"
    }

    #Create list vector to hold plots
    gvec <- NULL
    gvec <- vector("list", length = length(subset_points))
    pcatitbase <- paste(algorithm, "of", analysis)

    for (sp in 1:length(subset_points)){
        if (!(is.null(subsetby))){
            samplesToKeep <- rownames(pData(obj))[which((pData(obj)[, which(colnames(pData(obj)) == subsetby)]) == subset_points[sp])]
            pcatit <- paste(pcatitbase, "within", subset_points[sp])
            currobj <- obj[, samplesToKeep]
        } else {
            currobj <- obj
            samplesToKeep <- rownames(pData(currobj))
            pcatit <- pcatitbase
        }

        currobj <- filter_experiment(mgseqobj = currobj, featmaxatleastPPM = featmaxatleastPPM, featcutoff = featcutoff, samplesToKeep = samplesToKeep, asPA = FALSE, asPPM = TRUE, mgSeqnorm = mgSeqnorm)

        countmat <- MRcounts(currobj)

        if (ignoreunclassified == TRUE){
            dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified")
            rowsToKeep <- which(!(rownames(countmat) %in% dunno))
            countmat <- countmat[rowsToKeep, ]
        }

        #log2 transform if applicable
        if (log2tran == TRUE){
            countmat2 <- sapply(1:ncol(countmat), function(x){ countmat[, x] <- (log2(countmat[, x] + 1)) })
            colnames(countmat2) <- colnames(countmat)
            countmat <- countmat2
        }

        n <- nrow(countmat)
        comp <- 1:3
        rowsToKeep <- names(which(rowSums(countmat) > 0))
        countmat <- countmat[rowsToKeep, ]
        rowVars <- rowSds(countmat)
        countmat <- countmat[order(rowVars, decreasing = TRUE), ]
        if (transp == TRUE) {
            countmat <- t(countmat)
        }

        if (algorithm == "tSNE"){
            #tSNE algorithm
            permanova <- FALSE
            if (is.null(perplx)){
                perplx <- round(nrow(pData(currobj)) * 0.3, 0)
            }

            set.seed(4140)
            tsne_out <- Rtsne(countmat, dims = 3, initial_dims = 500, perplexity = perplx, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000)
            dford <- as.data.frame(tsne_out$Y)
            rownames(dford) <- rownames(pData(currobj))
            colnames(dford)[1:3] <- c("PC1", "PC2", "PC3")
            xl <- "tSNE 1"
            yl <- "tSNE 2"
            zl <- "tSNE 3"

        } else if (algorithm == "tUMAP"){
            permanova <- FALSE
            set.seed(4140)
            tumap_out <- tumap(countmat, n_components = 2, n_neighbors = 15, verbose = FALSE, n_threads = threads)
            dford <- as.data.frame(tumap_out)
            rownames(dford) <- rownames(pData(currobj))
            colnames(dford)[1:2] <- c("PC1", "PC2")
            xl <- "tUMAP 1"
            yl <- "tUMAP 2"
            #zl <- "tSNE 3"

        } else {
            #Not tSNE, so use PCA
            distfun <- stats::dist
            #d <- distfun(mat, method = "euclidian")
            d <- vegdist(countmat, method = "jaccard", na.rm = TRUE)
            pcaRes <- prcomp(d)
            ord <- pcaRes$x
            vars <- pcaRes$sdev^2
            vars <- round(vars/sum(vars), 5) * 100

            pnmetadata <- pData(currobj)
            cats <- pnmetadata[, colourby]

            if (!(is.numeric(cats))){
                permanovap <- vegan::adonis(as.formula(paste("d ~ ", colourby)), data = pData(currobj))$aov.tab$`Pr(>F)`[1]
            } else {
                #print("Impossible to get permanova because colourby is continuous")
                permanova <- FALSE
                permanovap <- 1
            }

            xl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[1]], vars[comp[1]])
            yl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[2]], vars[comp[2]])
            zl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[3]], vars[comp[3]])
            dford <- as.data.frame(ord[, comp])
        }
        #Add colour, size, shape
        dford$Colours <- pData(currobj)[match(rownames(dford), rownames(pData(currobj))), which(colnames(pData(currobj)) == colourby)]
        dford$Size <- pData(currobj)[match(rownames(dford), rownames(pData(currobj))), which(colnames(pData(currobj)) == sizeby)]
        dford$Shape <- pData(currobj)[match(rownames(dford), rownames(pData(currobj))), which(colnames(pData(currobj)) == shapeby)]
        dford$Pair <- pData(currobj)[match(rownames(dford), rownames(pData(currobj))), which(colnames(pData(currobj)) == pairby)]

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
                    if (permanovap < 0.05){
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

        if (!(is.null(plottit))){
            pcatit <- plottit
        }

        if ((permanova != FALSE) && (algorithm == "PCA")) {
            pcatit <- paste(pcatit, paste("p <", permanovap))
        }

        if (mgSeqnorm == TRUE){
            pcatit <- paste(pcatit, (paste0("MetagenomeSeq normalization = ", as.character(mgSeqnorm))), sep="\n")
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
