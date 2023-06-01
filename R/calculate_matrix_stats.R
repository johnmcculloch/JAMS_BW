#' calculate_matrix_stats(countmatrix = NULL, uselog = NULL, statsonlog = TRUE, stattype = NULL, classesvector = NULL, invertbinaryorder = FALSE, fun_for_l2fc = "geom_mean", threshPA = 0, numthreads = 1, nperm = 999)
#'
#' Returns a data frame with the statistics for a feature count matrix ordered by highest variance or lowest Mann-Whitney-Wilcoxon test between binary categories.
#'
#' @export

calculate_matrix_stats <- function(countmatrix = NULL, uselog = NULL, statsonlog = TRUE, stattype = NULL, classesvector = NULL, invertbinaryorder = FALSE, fun_for_l2fc = "geom_mean", threshPA = 0, numthreads = 1, nperm = 10000){

    #Test for silly stuff
    if ((stattype %in% c("binary", "permanova", "anova", "PA")) && (is.null(classesvector))){
        stop("If rows are to be selected by highest significant difference between classes in a discrete category, you must determine the category using the argument *classesvector*")
    }

    #Decide if we are comparing anything
    if (!(is.null(classesvector))){
        #Decide if variable is continuous or discrete
        if (is.numeric(classesvector)){
            numclass <- "continuous"
        } else {
            #Determine if there are two or more classes
            numclass <- length(unique(classesvector))
        }
    } else {
        #No comparing, use variance
        numclass <- 0
    }

    #If auto, decide what to do
    if (stattype %in% c("auto", "PA")){
        if (is.null(classesvector)){
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

    #Protect against rows with empty data
    #rowsToKeep <- which(rowSums(countmatrix) > 0)
    #countmatrix <- countmatrix[rowsToKeep, ]

    #Calculate matrix stats and get new matrix.
    if(stattype == "variance"){

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

        discretenames <- sort(unique(classesvector))
        comparisons <- lapply(1:nrow(countmatrix), function(x) { wilcox.test(countmatrix[x, ] ~ classesvector) })
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

            l2fc <- foldchange2logratio(foldchange(statscountsvecG1, statscountsvecG2), base = 2)
            if (is.na(l2fc)){
                l2fc <- 0
            }

            return(l2fc)
        }

        #l2fc <- sapply(1:nrow(countmatrix), function(x) { getl2fc(countsvec = countmatrix[x, ], classesvector = classesvector, discretenames = discretenames, countsinlog = uselog, method = "median") })
        l2fc <- sapply(1:nrow(countmatrix), function(x) { getl2fc(countsvec = countmatrix[x, ], classesvector = classesvector, discretenames = discretenames, countsinlog = uselog, method = fun_for_l2fc) })

        if (invertbinaryorder == TRUE){
            l2fc <- (l2fc * -1)
        }
        matstats <- data.frame(pval = as.vector(unlist(mwpval)), stat = as.vector(unlist(mwstat)), l2fc = l2fc)
        rownames(matstats) <- rownames(countmatrix)
        matstats$absl2fc <- abs(matstats$l2fc)

        #Adjust p-values
        adjlist <- lapply(p.adjust.methods, function(x){ p.adjust(matstats$pval, method = x, n = length(matstats$pval)) })
        names(adjlist) <- paste("padj", p.adjust.methods, sep = "_")
        adjdf <- as.data.frame(adjlist)
        matstats <- cbind(matstats, adjdf)
        matstats <- matstats[order(matstats$pval, decreasing = FALSE), ]
        matstats$Method <- rep("MannWhitneyWilcoxon", nrow(matstats))

    } else if (stattype %in% c("anova", "permanova")){
        require(parallel)

        if (stattype == "permanova"){
            flog.info("Calculating p-values with permanova. Please be patient...")
            flog.info(paste("Calculating p-values with PERMANOVA using", nperm, "permutations on", numthreads, "threads. Please be patient..."))
            mwpval <- sapply(1:nrow(countmatrix), function(x) { vegan::adonis(countmatrix[x,] ~ classesvector, parallel = numthreads, permutations = nperm)$aov.tab$`Pr(>F)`[1] } )
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
        adjlist <- lapply(p.adjust.methods, function(x){ p.adjust(matstats$pval, method = x, n = length(matstats$pval)) })
        names(adjlist) <- paste("padj", p.adjust.methods, sep="_")
        adjdf <- as.data.frame(adjlist)
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
        adjlist <- lapply(p.adjust.methods, function(x){ p.adjust(matstats$pval, method = x, n = length(matstats$pval)) })
        names(adjlist) <- paste("padj", p.adjust.methods, sep = "_")
        adjdf <- as.data.frame(adjlist)
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
            flog.info("Calculating pairwise correlation coefficients of all features")
            tcountmatrix <- t(countmatrix)
            matstats <- stats::cor(tcountmatrix, method = corstat)
        } else {
            comparisons <- lapply(1:nrow(countmatrix), function(x) { cor.test(countmatrix[x, ], classesvector, method = corstat) })
            mwstat <- sapply(1:length(comparisons), function(x) { unname(comparisons[[x]]$estimate) } )
            mwpval <- sapply(1:length(comparisons), function(x) { unname(comparisons[[x]]$p.value) } )

            matstats <- data.frame(correl = mwstat, pval = mwpval)
            matstats$abscorrel <- abs(matstats$correl)
            rownames(matstats) <- rownames(countmatrix)
            #Adjust p-values
            adjlist <- lapply(p.adjust.methods, function(x){ p.adjust(matstats$pval, method = x, n = length(matstats$pval)) })
            names(adjlist) <- paste("padj", p.adjust.methods, sep = "_")
            adjdf <- as.data.frame(adjlist)
            matstats <- cbind(matstats, adjdf)
            matstats <- matstats[order(matstats$pval, decreasing = FALSE), ]
            matstats$Method <- rep(corstat, nrow(matstats))
        }
    }

    return(matstats)
}
