#' Author: Phillip Burger
#' Date: Sept 3, 2014
#' Purpose: Get the ordinal rank of a number.
#' http://www.phillipburger.net/wordpress/ordinal-number-suffix-function-in-r
#' @export

getOrdinalNumber1 <- function(num) {
    result <- ""
    if (!(num %% 100 %in% c(11, 12, 13))) {
        result <- switch(as.character(num %% 10), 
            "1" = {paste0(num, "st")}, 
            "2" = {paste0(num, "nd")},
            "3" = {paste0(num, "rd")},
            paste0(num, "th"))
    } else {
        result <- paste0(num, "th")
    }
    result
}
