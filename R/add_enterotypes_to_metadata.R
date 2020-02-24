#' add_enterotypes_to_metadata(ptable = NULL, ExpObj = NULL, PPM_normalize_to_bases_sequenced = TRUE, force_n_clusters = NULL)
#' Defines enterotypes and attributes samples to enterotypes
#'
#' @export

add_enterotypes_to_metadata <- function(ptable = NULL, ExpObj = NULL, PPM_normalize_to_bases_sequenced = TRUE, force_n_clusters = NULL){

    require(ade4)
    require(cluster)
    require(clusterSim)

    analysis <- metadata(ExpObj)$analysis

    currobj <- filter_experiment(ExpObj = ExpObj, asPPM = TRUE, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL)


    if (is.null(ptable)){
        ptable <- as.data.frame(colData(currobj))
    }

    #Get rid of any taxa with 0 counts across samples and transform to straight relabund in ratios.
    enterotype_table_relabund <- assays(currobj)$BaseCounts
    enterotype_table_relabund <- enterotype_table_relabund / 1000000

    #Code from https://enterotype.embl.de
    #Define functions
    JSD <- function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
    KLD <- function(x, y) sum(x * log(x / y))
    dist.JSD <- function(inMatrix, pseudocount = 0.000001, ...) {
        KLD <- function(x, y) sum(x *log(x / y))
        JSD<- function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
        matrixColSize <- length(colnames(inMatrix))
        matrixRowSize <- length(rownames(inMatrix))
        colnames <- colnames(inMatrix)
        resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        inMatrix = apply(inMatrix,1:2,function(x) ifelse (x == 0, pseudocount, x))
        for(i in 1:matrixColSize) {
            for(j in 1:matrixColSize) {
                resultsMatrix[i, j] = JSD(as.vector(inMatrix[, i]),
                as.vector(inMatrix[, j]))
            }
        }
        colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
        as.dist(resultsMatrix)->resultsMatrix
        attr(resultsMatrix, "method") <- "dist"
        return(resultsMatrix)
    }

    #Apply functions
    enterotype_table_relabund.dist <- dist.JSD(enterotype_table_relabund)
    pam.clustering <- function(x, k) {

        cluster <- as.vector(pam(as.dist(x), k, diss = TRUE)$clustering)

        return(cluster)
    }

    #Start with k=3 and see how it goes
    enterotype_table_relabund.cluster <- pam.clustering(enterotype_table_relabund.dist, k = 3)

    nclusters <- index.G1(t(enterotype_table_relabund), enterotype_table_relabund.cluster, d = enterotype_table_relabund.dist, centrotypes = "medoids")
    nclusters <- NULL
    for (k in 1:20) {
        if (k==1) {
            nclusters[k] <- NA
        } else {
            enterotype_table_relabund.cluster_temp <- pam.clustering(enterotype_table_relabund.dist, k)
            nclusters[k] <- index.G1(t(enterotype_table_relabund), enterotype_table_relabund.cluster_temp,  d = enterotype_table_relabund.dist, centrotypes = "medoids")
        }
    }

    #Define which of nclusters is highest
    nclusters[which(is.na(nclusters))] <- 0
    optk <- which(nclusters == max(nclusters))
    #plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

    if (!is.null(force_n_clusters)){
        flog.warn(paste("Forcing clusterization into", force_n_clusters, "enterotypes"))
        flog.warn(paste("Calculated optimal number of enterotypes =", optk))
        optk <- force_n_clusters
    }

    #Redefine now with the best k
    enterotype_table_relabund.cluster <- pam.clustering(enterotype_table_relabund.dist, k = optk)

    noise.removal <- function(dataframe, percent = 0.01, top = NULL){
        dataframe -> Matrix
        bigones <- rowSums(Matrix) * 100 / (sum(rowSums(Matrix))) > percent
        Matrix_1 <- Matrix[bigones, ]
        print(percent)

        return(Matrix_1)
    }

    enterotype_table_relabund.denoized <- noise.removal(enterotype_table_relabund, percent = 0.01)
    obs.pca <- dudi.pca(data.frame(t(enterotype_table_relabund)), scannf = F, nf = 10)
    obs.bet <- bca(obs.pca, fac = as.factor(enterotype_table_relabund.cluster), scannf = F, nf = optk - 1)
    #s.class(obs.bet$ls, fac = as.factor(enterotype_table_relabund.cluster), grid = F)

    obs.pcoa <- dudi.pco(enterotype_table_relabund.dist, scannf = F, nf = optk)
    #s.class(obs.pcoa$li, fac = as.factor(enterotype_table_relabund.cluster), grid = F)
    #enterotit <- paste("Enterotypes by", analysis, sep=" ")
    #title(main = enterotit)
    #enteroplot <- recordPlot()

    sample2enterotype <- data.frame(Sample = rownames(obs.pcoa$li), Enterotype = as.factor(enterotype_table_relabund.cluster))
    varname <- paste(analysis, "Etype", sep="_")
    sample2enterotype$Enterotype <- paste(varname, sample2enterotype$Enterotype, sep = "_")

    #Add generated enterotype information back into pheno table
    ptable$Newcolumn <- sample2enterotype$Enterotype[match(ptable$Sample, sample2enterotype$Sample)]
    colnames(ptable)[which(colnames(ptable) == "Newcolumn")] <- varname

    return(ptable)
}
