#' plot_bias(bias, method = c("pca", tSNE"), frac = 0.25))
#'
#' Makes biased plots. Just kidding.
#' @export

plot_bias <- function(bias=NULL, method = "pca", frac = 0.25) {
    bias <- sample_by_category(bias, frac)
    data <- bias[,1:(ncol(bias)-1)]
    data <- mutate_all(data[,2:ncol(data)], function(y) as.numeric(as.character(y))) 
    result <- data_frame()
    if (method == "pca") {
        pca <- prcomp(data)
        result <- data.frame(x=pca$x[,1],y=pca$x[,2], category=bias[,ncol(bias)])
    } else if (method == "tSNE") {
        tsne <- Rtsne(dist(data)) 
        result <- data.frame(x=tsne$Y[,1],y=tsne$Y[,2], category=bias[,ncol(bias)])
    }
    title <- paste("Codon usage bias plot. Fraction=", frac, method)
    ucobiasplot <- ggplot(result, aes(x=x,y=y, color=category)) + geom_point() + ggtitle(title)

    return(ucobiasplot)
}
