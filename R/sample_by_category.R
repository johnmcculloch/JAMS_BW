#' sample_by_category
#'
#' Sample by category
#' @export

sample_by_category <- function(rows, frac = 0.25) {
    #counts <- table(rows$category) %>% .[.>= n]
    counts <- table(rows$category)
    result <- data.frame()
    for (name in names(counts)) {
        result <- rbind(result, sample_frac(rows[rows$category==name,], frac)) 
    }

    return(result)
}
