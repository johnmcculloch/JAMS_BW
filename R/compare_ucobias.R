#' compare_ucobias
#'
#' Compares taxa by codon usage bias.
#' @export

compare_ucobias<-function(featuredata=NULL, contigsdata=NULL, taxlevel="LKT", taxaofinterest=NULL, ucobias=NULL, method="tSNE", frac=0.25){

    taxlevels<-colnames(opt$LKTdose)[colnames(opt$LKTdose) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT")]

    if(is.null(taxaofinterest)){
        taxaofinterest<-unique(featuredata$LKT)
    }

    if(taxlevel != "LKT"){
        tt<-contigsdata[, taxlevels]
        featuredata<-left_join(featuredata, tt, by = "LKT")
    }

    featuresIwant<-subset(featuredata, (get(taxlevel) %in% taxaofinterest))

    cats<-featuresIwant[,c("Feature", taxlevel)]
    ucobiasinterest<-subset(ucobias, Feature %in% featuresIwant$Feature)
    ucobiasinterest<-left_join(ucobiasinterest, cats)
    ucobiasinterest$Feature<-NULL
    colnames(ucobiasinterest)[which(colnames(ucobiasinterest) == taxlevel)]<-"category"

    bias <- sample_by_category(ucobiasinterest, frac)
    ucodata <- bias[,1:(ncol(bias)-1)]
    ucodata <- mutate_all(ucodata[,2:ncol(ucodata)], function(y) as.numeric(as.character(y))) 
    if (method == "pca") {
        pca <- prcomp(ucodata)
        ucoresult <- data.frame(x=pca$x[,1],y=pca$x[,2], category=bias[,ncol(bias)])
    } else if (method == "tSNE") {
        tsne <- Rtsne(dist(ucodata)) 
        ucoresult <- data.frame(x=tsne$Y[,1],y=tsne$Y[,2], category=bias[,ncol(bias)])
    }
    title <- paste("Codon usage bias plot. Fraction=", frac, method)
    ucobiasplot <- ggplot(ucoresult, aes(x=x,y=y, color=category)) + geom_point() + ggtitle(title)

    return(ucobiasplot)
}
