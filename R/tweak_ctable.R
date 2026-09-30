#' tweak_ctable(ctable = NULL, pheno = NULL, outfile = NULL, ncol = 3, launch_browser = TRUE)
#'
#' Interactive Shiny colour editor for a JAMS colour table (ctable), e.g. the output of
#' make_colour_dictionary(..., object_to_return = "ctable"). Presents one colour-wheel
#' picker per class and returns the edited ctable to the R session when you click Done.
#'
#' If a `pheno` metadata data frame is supplied, the pickers are organised into a titled box
#' per metadata variable, showing exactly which classes will be plotted together. Because a
#' ctable holds a single colour per class name, a class that occurs in more than one variable
#' (e.g. "Pittsburgh" in both Cohort and Study) appears in each of those boxes and is kept in
#' sync as you edit - editing it anywhere updates it everywhere, which is precisely how the
#' flat ctable will behave on your plots. Any ctable class not found in the supplied metadata
#' is collected into a trailing "not found in metadata" box so nothing is hidden.
#'
#' This function depends on the 'shiny' and 'colourpicker' packages, which are deliberately
#' NOT JAMS dependencies (this tool is used rarely). If either package is missing, a polite
#' message is printed and the ctable is returned unchanged.
#'
#' @param ctable A ctable data frame (columns Name/Colour/Hex, or Class_label/Class_colour), OR a path to a tsv of one.
#' @param pheno Optional metadata data frame (or path to a tsv). If supplied, pickers are grouped by variable. Columns whose values are not ctable classes (IDs, dates, continuous variables) are skipped automatically.
#' @param outfile Optional path. If set, the edited ctable is also written to this tsv on Done.
#' @param ncol Number of colour pickers per row within each box. Default 3.
#' @param launch_browser Logical. If TRUE (default) opens in the web browser; if FALSE, uses the default Shiny viewer (e.g. RStudio pane).
#'
#' @return The edited ctable data frame (same column layout as the input). On Cancel, the original ctable is returned unchanged.
#' @export

tweak_ctable <- function(ctable = NULL, pheno = NULL, outfile = NULL, ncol = 3, launch_browser = TRUE){

    ############################################################
    ## 1. Dependency gate (polite, non-fatal)
    ############################################################
    needed <- c("shiny", "colourpicker")
    missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0){
        message(paste0(
            "tweak_ctable() needs the following package(s), which are NOT JAMS dependencies: ",
            paste0(missing_pkgs, collapse = ", "), ".\n",
            "Install with:  install.packages(c(",
            paste0("\"", missing_pkgs, "\"", collapse = ", "), "))\n",
            "then try again. Returning the ctable unchanged for now."))
        return(ctable)
    }

    ############################################################
    ## 2. Input handling + flexible column resolution
    ############################################################
    if (is.character(ctable) && length(ctable) == 1){
        ctable <- read.table(ctable, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
    }
    if (is.null(ctable) || nrow(ctable) < 1){
        message("tweak_ctable(): no ctable rows to edit. Returning input unchanged.")
        return(ctable)
    }
    ctable <- as.data.frame(ctable, stringsAsFactors = FALSE)

    name_col <- intersect(c("Name", "Class_label"), colnames(ctable))[1]
    col_col  <- intersect(c("Hex", "Colour", "Class_colour"), colnames(ctable))[1]
    if (is.na(name_col) || is.na(col_col)){
        message("tweak_ctable(): could not find a name column (Name/Class_label) and a colour column (Hex/Colour/Class_colour). Returning input unchanged.")
        return(ctable)
    }

    classes    <- as.character(ctable[[name_col]])
    start_cols <- as.character(ctable[[col_col]])
    start_of   <- setNames(start_cols, classes)   # class -> starting colour

    ############################################################
    ## 2b. Optional pheno for variable grouping
    ############################################################
    grouped <- !is.null(pheno)
    if (grouped){
        if (is.character(pheno) && length(pheno) == 1){
            pheno <- read.table(pheno, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
        }
        pheno <- as.data.frame(pheno, stringsAsFactors = FALSE)
    }

    ############################################################
    ## 3. Build the picker map
    ##    - flat mode: one picker per class (canonical layout, unchanged behaviour)
    ##    - grouped mode: one picker per (variable, class); shared classes recur and sync
    ############################################################
    # picker_map: one row per rendered picker.
    picker_map <- data.frame(id = character(0), class = character(0), variable = character(0),
                             colour = character(0), stringsAsFactors = FALSE)
    gid <- 0
    var_order <- character(0)   # display order of boxes

    if (!grouped){
        for (i in seq_along(classes)){
            gid <- gid + 1
            picker_map <- rbind(picker_map, data.frame(id = paste0("cid", gid), class = classes[i],
                                                       variable = "__all__", colour = start_cols[i],
                                                       stringsAsFactors = FALSE))
        }
        var_order <- "__all__"
    } else {
        shown_classes <- character(0)
        for (v in colnames(pheno)){
            if (v == "Sample") next   # sample names are never colour classes
            vals <- unique(as.character(pheno[[v]]))
            cls  <- sort(classes[classes %in% vals])
            if (length(cls) < 1) next
            var_order <- c(var_order, v)
            for (cn in cls){
                gid <- gid + 1
                picker_map <- rbind(picker_map, data.frame(id = paste0("cid", gid), class = cn,
                                                           variable = v, colour = unname(start_of[cn]),
                                                           stringsAsFactors = FALSE))
            }
            shown_classes <- union(shown_classes, cls)
        }
        leftover <- setdiff(classes, shown_classes)
        if (length(leftover) > 0){
            var_order <- c(var_order, "__leftover__")
            for (cn in sort(leftover)){
                gid <- gid + 1
                picker_map <- rbind(picker_map, data.frame(id = paste0("cid", gid), class = cn,
                                                           variable = "__leftover__", colour = unname(start_of[cn]),
                                                           stringsAsFactors = FALSE))
            }
        }
    }

    # class -> its picker ids, and a canonical (first) id per class.
    ids_by_class <- split(picker_map$id, picker_map$class)
    canonical_id <- vapply(names(ids_by_class), function(cn){ ids_by_class[[cn]][1] }, character(1))

    ############################################################
    ## 4. UI helpers
    ############################################################
    colwidth <- max(1, floor(12 / ncol))

    make_box <- function(v){
        rows_map <- picker_map[picker_map$variable == v, , drop = FALSE]
        pk <- lapply(seq_len(nrow(rows_map)), function(k){
            colourpicker::colourInput(inputId = rows_map$id[k], label = rows_map$class[k],
                                      value = rows_map$colour[k], showColour = "both", returnName = FALSE)
        })
        groups <- split(seq_along(pk), ceiling(seq_along(pk) / ncol))
        body <- lapply(groups, function(r){
            shiny::fluidRow(lapply(r, function(k){ shiny::column(width = colwidth, pk[[k]]) }))
        })
        title <- if (v == "__all__") NULL
                 else if (v == "__leftover__") "Classes not found in supplied metadata"
                 else v
        if (is.null(title)){
            do.call(shiny::tagList, body)
        } else {
            shiny::wellPanel(shiny::h4(title), do.call(shiny::tagList, body))
        }
    }

    boxes_ui <- lapply(var_order, make_box)

    intro <- if (grouped){
        "Each box is a metadata variable; the classes inside it are what get plotted together, so their colours should be mutually distinct. A class shared by several variables (e.g. Pittsburgh in Cohort and Study) appears in each of its boxes and is kept in sync \u2014 editing it anywhere updates it everywhere, because a ctable stores one colour per class name. The report below flags any within-variable duplicates."
    } else {
        "Adjust any class colour with its picker. The report below updates live and lists any classes that currently share a hex \u2014 fix those to avoid two categories looking identical on a plot."
    }

    ui <- shiny::fluidPage(
        shiny::titlePanel("Tweak JAMS colour table"),
        shiny::fluidRow(shiny::column(12,
            shiny::helpText(intro),
            shiny::verbatimTextOutput("dupreport"),
            shiny::hr()
        )),
        do.call(shiny::tagList, boxes_ui),
        shiny::hr(),
        shiny::fluidRow(shiny::column(12,
            shiny::actionButton("done", "Done \u2013 return edited ctable", class = "btn-primary"),
            shiny::actionButton("cancel", "Cancel \u2013 return original")
        ))
    )

    ############################################################
    ## 5. Server
    ############################################################
    server <- function(input, output, session){

        # Keep every picker of a shared class in sync. The identical() guard makes the
        # update cascade converge (an update to the same value is a no-op).
        for (cn in names(ids_by_class)){
            idv <- ids_by_class[[cn]]
            if (length(idv) < 2) next
            lapply(idv, function(this_id){
                shiny::observeEvent(input[[this_id]], {
                    v <- input[[this_id]]
                    for (other in setdiff(idv, this_id)){
                        if (!identical(input[[other]], v)){
                            colourpicker::updateColourInput(session, other, value = v)
                        }
                    }
                }, ignoreInit = TRUE)
            })
        }

        # Current colour per unique class, read from its canonical picker.
        current_colours <- function(){
            setNames(vapply(names(ids_by_class), function(cn){
                v <- input[[ canonical_id[[cn]] ]]
                if (is.null(v) || !nzchar(v)) unname(start_of[cn]) else v
            }, character(1)), names(ids_by_class))
        }

        output$dupreport <- shiny::renderText({
            cur <- current_colours()
            if (grouped){
                lines <- character(0)
                report_vars <- setdiff(var_order, "__leftover__")
                for (v in report_vars){
                    cls <- picker_map$class[picker_map$variable == v]
                    hexes <- toupper(cur[cls])
                    dh <- unique(hexes[duplicated(hexes)])
                    for (h in dh){
                        lines <- c(lines, paste0("  [", v, "]  ", h, "  <-  ",
                                                 paste0(cls[hexes == h], collapse = ", ")))
                    }
                }
                if (length(lines) == 0){
                    "\u2713 No within-variable duplicate colours: inside every metadata variable, each class has a distinct colour."
                } else {
                    paste0("WITHIN-VARIABLE DUPLICATES (these classes would look identical on the same plot):\n",
                           paste0(lines, collapse = "\n"))
                }
            } else {
                cur <- toupper(cur)
                dup_hexes <- unique(cur[duplicated(cur)])
                if (length(dup_hexes) == 0){
                    "\u2713 No duplicate colours: every class currently has a unique colour."
                } else {
                    msgs <- vapply(dup_hexes, function(h){
                        paste0("  ", h, "  <-  ", paste0(names(cur)[cur == h], collapse = ", "))
                    }, character(1))
                    paste0("DUPLICATE COLOURS (classes sharing a hex):\n", paste0(msgs, collapse = "\n"))
                }
            }
        })

        gather <- function(){
            cur <- current_colours()
            out <- ctable
            newcols <- unname(cur[as.character(out[[name_col]])])
            # Any ctable row not represented by a picker (shouldn't happen) keeps its old colour.
            miss <- is.na(newcols)
            if (any(miss)) newcols[miss] <- out[[col_col]][miss]
            out[[col_col]] <- newcols
            if ("Hex" %in% colnames(out))    out[["Hex"]]    <- newcols
            if ("Colour" %in% colnames(out)) out[["Colour"]] <- newcols
            rownames(out) <- as.character(out[[name_col]])
            out
        }

        shiny::observeEvent(input$done, {
            out <- gather()
            if (!is.null(outfile)){
                utils::write.table(out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
                message(paste("tweak_ctable(): wrote edited ctable to", outfile))
            }
            shiny::stopApp(returnValue = out)
        })
        shiny::observeEvent(input$cancel, {
            shiny::stopApp(returnValue = ctable)
        })
    }

    ############################################################
    ## 6. Launch (blocks until Done/Cancel, then returns the value)
    ############################################################
    app <- shiny::shinyApp(ui = ui, server = server)
    shiny::runApp(app, launch.browser = launch_browser)
}
