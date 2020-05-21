#' print_table(tb=NULL, tabletitle=NULL, fontsize=6, numrows=20)
#'
#' JAMSalpha function
#' @export

print_table <- function(tb = NULL, tabletitle = NULL, fontsize = 6, numrows = 20){

    if (class(tb) != "data.frame"){
        stop("Table must be a data frame")
    }

    #Define function to wrap text
    wrap_cell_text <- function(celltext=NULL, width=30){
        celltextsplit <- strwrap(celltext, width = width, simplify = FALSE)
        wrappedtext <- sapply(celltextsplit, paste, collapse = "\n")

        return(wrappedtext)
    }

    #Maximum number of characters in a line should be about 160 for font of size 6
    maxchar <- as.integer((fontsize * ((160-80)/(6-12))) + 240)
    maxcharincell <- as.integer(maxchar / ncol(tb))
    tb[] <- lapply(tb, function (x) { wrap_cell_text(celltext = x, width = maxcharincell) })

    if (nrow(tb) > numrows){
        #split table into n rows at a time
        from <- seq(from = 1, to = nrow(tb), by = numrows)
        to <- seq(from = 1, to = nrow(tb), by = numrows) + (numrows - 1)
        for (p in 1:(length(from) - 1)){
            tbt <- tb[(from[p]:to[p]), ]
            table <- tableGrob(tbt, rows = NULL, theme = ttheme_default(base_size = fontsize))
            title <- textGrob(tabletitle, gp=gpar(fontsize = 10))
            padding <- unit(0.5, "line")
            table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0)
            table <- gtable_add_grob(table, list(title), t=c(1), l=c(1), r=ncol(table))
            grid.newpage()
            grid.draw(table)
        }
        tbt <- tb[(from[(length(from))]:nrow(tb)), ]
        table <- tableGrob(tbt, rows = NULL, theme = ttheme_default(base_size = fontsize))
        title <- textGrob(tabletitle, gp=gpar(fontsize=10))
        padding <- unit(0.5, "line")
        table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0)
        table <- gtable_add_grob(table, list(title), t=c(1), l=c(1), r=ncol(table))
        grid.newpage()
        grid.draw(table)
    } else {
        table <- tableGrob(tb, rows = NULL, theme = ttheme_default(base_size = fontsize))
        title <- textGrob(tabletitle, gp=gpar(fontsize=10))
        padding <- unit(0.5,"line")
        table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0)
        table <- gtable_add_grob(table, list(title), t=c(1), l=c(1), r=ncol(table))
        grid.newpage()
        grid.draw(table)
    }
}
