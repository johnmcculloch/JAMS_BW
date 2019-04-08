#' shrink_perbasecoverage(perbasecoverage=NULL, percentage=2)
#' 
#' #Subset down to n% of the dataset to make it calculatable. To be used only by JAMSalpha.
#' @export

shrink_perbasecoverage<-function(perbasecoverage=NULL, percentage=2){
    numbases <- 1:nrow(perbasecoverage)
    breaks <- numbases[seq(1, length(numbases), (100/percentage))]
    perbasecoverage_reduced<-perbasecoverage[breaks, ]

    return(perbasecoverage_reduced)
}
