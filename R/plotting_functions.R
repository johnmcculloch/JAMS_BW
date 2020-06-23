#' make_featureplot()
#'
#' JAMSbeta wrapper for returning a ggplot boxplot, violin plot or scatterplot from a data frame. This function was concieved to work within other main plotting functions in JAMS. Do not feel frustrated if it does not work out of this context. My playground, my rules.
#' @export


make_featureplot <- function(dat = NULL, logtransbase = NULL, subsetby = NULL, compareby = NULL, colourby = NULL, shapeby = NULL, facetby = NULL, wrap_facet = FALSE, overlay_boxplot = FALSE){
    require(ggplot2)

    if (class(dat)[1] != "data.frame"){
        stop("Input data must be a data frame")
    }





    return(p)
}
