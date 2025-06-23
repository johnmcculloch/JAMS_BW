#' filter_correlations(corrmat = NULL, mincorrelcoeff = NULL)
#' Given a pairwise correlation matrix, eliminate features which do not present an absolute correlation coefficient smaller than mincorrelcoeff with all other features other than itself.
#'
#' @export

filter_correlations <- function(corrmat = NULL, mincorrelcoeff = NULL){

     if(nrow(corrmat) != ncol(corrmat)){
         stop("Correlation matrix must have equal numbers of rows and columns.")
     }

     featsIwant <- NULL

     for (rw in 1:nrow(corrmat)){
         featint <- rownames(corrmat)[rw]
         #print(paste("Checking:", featint))
         correlations <- corrmat[which(rownames(corrmat) != featint), featint]

         if(max(abs(correlations)) >= mincorrelcoeff){
             feat <- featint
         } else {
             feat <- NULL
         }

         featsIwant <- append(featsIwant, feat)

     }

     corrmat <- corrmat[featsIwant, featsIwant]

     return(corrmat)
 }
