#' name_samples(list.data=NULL)
#'
#' Given a list.data object, returns a vector of sample names present in the object.
#' @export

name_samples <- function(list.data = NULL){
    loadedsamples<-gsub("_contigsdata", "", (names(list.data)[grep("_contigsdata", names(list.data))]))

    return(loadedsamples)
}
