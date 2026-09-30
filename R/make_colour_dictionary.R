#' make_colour_dictionary(variable_list = NULL, pheno = NULL, phenolabels = NULL, columnsToUse = NULL, columnsToIgnore = NULL, class_to_ignore = "N_A", colour_of_class_to_ignore = "#bcc2c2", legacy_ctable = NULL, colour_table = NULL, within_variable_palette = TRUE, diverging_palette = "Spectral", binary_positive_colour = "#0004e9", binary_negative_colour = "#E41A1C", extra_positive_tokens = NULL, extra_negative_tokens = NULL, guess_binary_polarity = TRUE, avoid_collisions = TRUE, shuffle = FALSE, object_to_return = "ctable", warn_on_collisions = TRUE)
#'
#' Returns a colour dictionary (either a legacy-style ctable data frame or a per-variable cdict list)
#' for the classes within the discrete variables of a phenotable.
#'
#' Colour assignment precedence, applied per variable and per class:
#'   1. legacy_ctable override (if the class name is found there)
#'   2. within-variable generated palette (diverging, harmonious)
#'   3. class_to_ignore is always painted colour_of_class_to_ignore
#'
#' For two-class variables, an OPTIONAL semantic guess (see guess_binary_polarity) tries to
#' place an "affirmative/positive" class (Responder, Good, Yes, Y, CR, PR, Sensitive, ...) on
#' binary_positive_colour and a "negative" class (Non_Responder, Bad, No, N, PD, Resistant,
#' Refractory, ...) on binary_negative_colour. This is a heuristic and can be semantically
#' wrong for valence-neutral pairs (e.g. Death_event Y/N, where Y is the bad outcome). Tokens
#' that are genuinely context-dependent (High/Low, Present/Absent) are deliberately NOT guessed
#' and fall back to a stable alphabetical assignment. Legacy overrides always win over the guess.
#'
#' When object_to_return = "ctable" and avoid_collisions = TRUE (the default), a preventive pass
#' guarantees that WITHIN EACH coloured variable every class resolves to a distinct colour in the
#' flat ctable. This fixes the classic artefact where two classes of one variable inherit the same
#' colour after the global collapse (e.g. two 5-class variables both taking the Spectral midpoint
#' #FFFFBF, one class then overwriting another during de-duplication). Legacy- and ignore-locked
#' colours are never moved; when a clash must be resolved the more "local" class (present in the
#' fewest variables) is nudged, and only away from its own siblings, so intended cross-variable
#' colour reuse (e.g. every binary-positive class sharing one blue) is preserved.
#'
#' @param variable_list Optional pre-computed output of define_kinds_of_variables. If NULL, it is derived from pheno/phenolabels.
#' @param pheno The curated metadata data frame (rownames = sample names).
#' @param phenolabels Optional phenolabels data frame. If NULL and variable_list is NULL, column types are imputed.
#' @param columnsToUse Optional character vector of metadata columns to colour. If supplied, only these (plus "Sample", which is always protected) are considered. Overrides the discrete-variable auto-detection.
#' @param columnsToIgnore Optional character vector of metadata columns to NOT colour. "Sample" is always ignored regardless.
#' @param class_to_ignore String (or vector) of class values to treat as missing/ignored. Painted colour_of_class_to_ignore. Default "N_A".
#' @param colour_of_class_to_ignore Hex colour for class_to_ignore. Default grey "#bcc2c2".
#' @param legacy_ctable Optional data frame OR path to a tsv of legacy colours to override generated colours where class names match. Accepts either columns Name/Hex (or Name/Colour) or the older Class_label/Class_colour. See details.
#' @param colour_table Deprecated alias for legacy_ctable, retained for backwards compatibility.
#' @param within_variable_palette Logical. If TRUE (default), each variable gets its own diverging palette. If FALSE, the historical global "atom" recycling behaviour is used.
#' @param diverging_palette String naming an RColorBrewer diverging palette used for variables with > 2 classes. Default "Spectral".
#' @param binary_positive_colour Hex colour for the affirmative/positive class of a 2-class variable. Default "#0004e9" (blue), matching the common Responder = blue legacy convention.
#' @param binary_negative_colour Hex colour for the negative class of a 2-class variable. Default "#E41A1C" (red).
#' @param extra_positive_tokens Optional character vector of additional class labels to treat as "positive". They are normalised (lowercased, non-alphanumerics stripped) before matching, so "Complete Response" and "complete_response" are equivalent.
#' @param extra_negative_tokens Optional character vector of additional class labels to treat as "negative", normalised the same way.
#' @param guess_binary_polarity Logical. If TRUE (default), apply the semantic positive/negative guess to 2-class variables. If FALSE, all 2-class variables use the plain alphabetical fallback (first class = binary_positive_colour, second = binary_negative_colour).
#' @param avoid_collisions Logical. If TRUE (default) and object_to_return = "ctable", run the preventive per-variable de-collision pass described above. Has no effect on the "cdict" return (per-variable dictionaries are collision-free by construction).
#' @param shuffle Logical. If TRUE, shuffles palette order. Default FALSE.
#' @param object_to_return Either "ctable" (a data frame with Name, Colour, Hex) or "cdict" (a named list, one data frame per variable). Default "ctable".
#' @param warn_on_collisions Logical. If TRUE (default), warns when the same class name appears in more than one coloured variable with differing generated colours before collapsing to the global ctable.
#'
#' @details
#' Reading a legacy ctable: the file/data frame is column-name agnostic. It looks for a
#' class-name column among c("Name", "Class_label") and a colour column among
#' c("Hex", "Colour", "Class_colour"), in that order of preference. If a legacy_ctable is
#' supplied but none of its class names match any class in the (used) metadata, a warning is
#' emitted rather than silently doing nothing.
#'
#' @export

make_colour_dictionary <- function(variable_list = NULL, pheno = NULL, phenolabels = NULL,
                                   columnsToUse = NULL, columnsToIgnore = NULL,
                                   class_to_ignore = "N_A", colour_of_class_to_ignore = "#bcc2c2",
                                   legacy_ctable = NULL, colour_table = NULL,
                                   within_variable_palette = TRUE, diverging_palette = "Spectral",
                                   binary_positive_colour = "#0004e9", binary_negative_colour = "#E41A1C",
                                   extra_positive_tokens = NULL, extra_negative_tokens = NULL,
                                   guess_binary_polarity = TRUE, avoid_collisions = TRUE,
                                   shuffle = FALSE, object_to_return = "ctable",
                                   warn_on_collisions = TRUE){

    require(RColorBrewer)

    ################################################################
    ## 0. Backwards compatibility: fold deprecated colour_table into legacy_ctable
    ################################################################
    if (!is.null(colour_table) && is.null(legacy_ctable)){
        flog.info("colour_table is deprecated; treating it as legacy_ctable.")
        legacy_ctable <- colour_table
    }

    ################################################################
    ## 1. Build the JAMS "atom" palette (kept for the non-within-variable path)
    ################################################################
    atoms <- c("0200FC", "ff7f02", "468B00", "CD2626", "267ACD", "640032", "FC00FA", "000000", "478559", "161748", "f95D9B", "39A0CA")
    JAMSpalette <- paste0("#", atoms)
    for (bset in c("Set1", "Paired", "Set2", "Set3")){
        JAMSpalette <- c(JAMSpalette, colorRampPalette(brewer.pal(8, bset))(length(atoms)))
    }
    if (shuffle){
        JAMSpalette <- sample(JAMSpalette)
    }

    ################################################################
    ## 2. Resolve which variables to colour
    ################################################################
    if (is.null(variable_list)){
        if (is.null(phenolabels)){
            flog.info("Imputing types of columns on the metadata.")
            Var_label <- colnames(pheno)
            Var_type <- sapply(Var_label, function (x) { infer_column_type(phenotable = pheno, colm = x, class_to_ignore = class_to_ignore) } )
            phenolabels <- data.frame(Var_label = unname(Var_label), Var_type = unname(Var_type), stringsAsFactors = FALSE)
        }
        variable_list <- define_kinds_of_variables(phenolabels = phenolabels, phenotable = pheno, maxclass = 30, maxsubclass = 100, class_to_ignore = class_to_ignore, verbose = FALSE)
    }

    candidate_vars <- unique(unname(unlist(variable_list[which(!(names(variable_list) %in% c("sample", "continuous")))])))

    # Apply columnsToUse / columnsToIgnore, ALWAYS protecting Sample.
    if (!is.null(columnsToUse)){
        columnsToUse <- setdiff(unique(columnsToUse), "Sample")
        missing_use <- columnsToUse[!(columnsToUse %in% colnames(pheno))]
        if (length(missing_use) > 0){
            flog.warn(paste("columnsToUse entries not found in metadata and skipped:", paste0(missing_use, collapse = ", ")))
        }
        discretes <- columnsToUse[columnsToUse %in% colnames(pheno)]
    } else {
        discretes <- candidate_vars
    }
    force_ignore <- unique(c("Sample", columnsToIgnore))
    discretes <- setdiff(discretes, force_ignore)

    if (length(discretes) < 1){
        stop("After applying columnsToUse/columnsToIgnore (and always protecting Sample), there are no variables left to colour.")
    }

    ################################################################
    ## 3. Read and normalise the legacy ctable, if supplied
    ################################################################
    legacy_lookup <- NULL  # named character vector: names = class, values = hex
    if (!is.null(legacy_ctable)){
        if (is.character(legacy_ctable) && length(legacy_ctable) == 1){
            legacy_ctable <- read.table(legacy_ctable, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
        }
        legacy_ctable <- as.data.frame(legacy_ctable, stringsAsFactors = FALSE)
        name_col <- intersect(c("Name", "Class_label"), colnames(legacy_ctable))[1]
        col_col  <- intersect(c("Hex", "Colour", "Class_colour"), colnames(legacy_ctable))[1]
        if (is.na(name_col) || is.na(col_col)){
            flog.warn("legacy_ctable supplied but could not find a class-name column (Name/Class_label) and a colour column (Hex/Colour/Class_colour). Ignoring it.")
        } else {
            legacy_lookup <- as.character(legacy_ctable[[col_col]])
            names(legacy_lookup) <- as.character(legacy_ctable[[name_col]])
            legacy_lookup <- legacy_lookup[!duplicated(names(legacy_lookup))]
        }
    }

    ################################################################
    ## 4. Helpers: normalise labels + polarity guess + palette + colour nudge
    ################################################################
    norm_tok <- function(x){ gsub("[^a-z0-9]", "", tolower(as.character(x))) }

    positive_tokens <- c("responder", "response", "responsive", "good", "yes", "y",
        "positive", "pos", "present", "cr", "pr", "completeresponse", "partialresponse",
        "alive", "survivor", "surviving", "sensitive", "improved", "improvement",
        "remission", "effective", "benefit", "clinicalbenefit", "durable",
        "favourable", "favorable", "true")
    negative_tokens <- c("nonresponder", "nonresponse", "nonresponsive", "bad", "no", "n",
        "negative", "neg", "absent", "pd", "progressivedisease", "progression",
        "progressive", "progressed", "dead", "death", "deceased", "resistant",
        "refractory", "relapse", "relapsed", "worsened", "worsening", "ineffective",
        "nonbenefit", "nobenefit", "unfavourable", "unfavorable", "false")

    if (!is.null(extra_positive_tokens)) positive_tokens <- unique(c(positive_tokens, norm_tok(extra_positive_tokens)))
    if (!is.null(extra_negative_tokens)) negative_tokens <- unique(c(negative_tokens, norm_tok(extra_negative_tokens)))

    classify_polarity <- function(cn){
        t <- norm_tok(cn)
        if (t %in% positive_tokens) return("positive")
        if (t %in% negative_tokens) return("negative")
        return(NA_character_)
    }

    assign_binary <- function(classnames){
        cols <- setNames(rep(NA_character_, 2), classnames)
        decided <- FALSE
        if (guess_binary_polarity){
            pol <- vapply(classnames, classify_polarity, character(1))
            pos_idx <- which(pol == "positive")
            neg_idx <- which(pol == "negative")
            if (length(pos_idx) == 1 && length(neg_idx) == 1){
                cols[classnames[pos_idx]] <- binary_positive_colour
                cols[classnames[neg_idx]] <- binary_negative_colour
                decided <- TRUE
            } else if (length(pos_idx) == 1 && length(neg_idx) == 0){
                cols[classnames[pos_idx]] <- binary_positive_colour
                cols[classnames[-pos_idx]] <- binary_negative_colour
                decided <- TRUE
            } else if (length(neg_idx) == 1 && length(pos_idx) == 0){
                cols[classnames[neg_idx]] <- binary_negative_colour
                cols[classnames[-neg_idx]] <- binary_positive_colour
                decided <- TRUE
            }
        }
        if (!decided){
            s <- sort(classnames)
            cols[s[1]] <- binary_positive_colour
            cols[s[2]] <- binary_negative_colour
        }
        return(cols)
    }

    generate_variable_palette <- function(classnames){
        n <- length(classnames)
        if (n == 0) return(setNames(character(0), character(0)))
        if (n == 1){
            return(setNames(binary_positive_colour, classnames))
        } else if (n == 2){
            return(assign_binary(classnames))
        } else {
            maxbrew <- tryCatch(brewer.pal.info[diverging_palette, "maxcolors"], error = function(e) 11)
            if (is.na(maxbrew)) maxbrew <- 11
            base_n <- min(max(3, n), maxbrew)
            base_cols <- brewer.pal(base_n, diverging_palette)
            cols <- colorRampPalette(base_cols)(n)
            return(setNames(cols, classnames))
        }
    }

    # Perturb a hex colour, in HSV space, to the nearest hex not already in `used`.
    # Deterministic (structured search, no RNG). Returns the original if no free slot found
    # or the colour is unparseable. Base-R only (grDevices), so no new dependency.
    nudge_colour <- function(hex, used){
        used <- toupper(used)
        base <- tryCatch(grDevices::col2rgb(hex), error = function(e) NULL)
        if (is.null(base)) return(hex)
        hsv0 <- grDevices::rgb2hsv(base)[, 1]
        h0 <- hsv0["h"]; s0 <- hsv0["s"]; v0 <- hsv0["v"]
        hue_steps <- c(0.04, -0.04, 0.08, -0.08, 0.12, -0.12, 0.16, -0.16, 0.20, -0.20, 0.25, -0.25, 0.33, -0.33, 0.5)
        val_steps <- c(0, -0.12, 0.12, -0.24, 0.24, -0.36)
        sat_steps <- c(0, 0.18, -0.18, 0.36)
        for (ss in sat_steps){
            for (vv in val_steps){
                for (hh in hue_steps){
                    h <- (h0 + hh) %% 1
                    s <- min(1, max(0.15, s0 + ss))
                    v <- min(1, max(0.15, v0 + vv))
                    cand <- toupper(grDevices::hsv(h, s, v))
                    if (!(cand %in% used)) return(cand)
                }
            }
        }
        return(hex)
    }

    ################################################################
    ## 5. Build per-variable colour dictionaries
    ################################################################
    cdict <- list()
    for (d in discretes){
        classnames <- sort(unique(as.character(pheno[[d]])))

        if (within_variable_palette){
            paint_classes <- classnames[!(classnames %in% class_to_ignore)]
            colvec <- generate_variable_palette(paint_classes)
            if (any(classnames %in% class_to_ignore)){
                greys <- setNames(rep(colour_of_class_to_ignore, sum(classnames %in% class_to_ignore)),
                                  classnames[classnames %in% class_to_ignore])
                colvec <- c(colvec, greys)
            }
        } else {
            pal <- JAMSpalette
            if (length(classnames) > length(pal)){
                pal <- rep(pal, ceiling(length(classnames) / length(pal)))
            }
            colvec <- setNames(pal[seq_along(classnames)], classnames)
            colvec[classnames %in% class_to_ignore] <- colour_of_class_to_ignore
        }

        # Legacy override wins wherever a class name matches.
        if (!is.null(legacy_lookup)){
            hits <- intersect(names(colvec), names(legacy_lookup))
            if (length(hits) > 0){
                colvec[hits] <- legacy_lookup[hits]
            }
        }

        ctab <- data.frame(Name = names(colvec), Colour = unname(colvec), stringsAsFactors = FALSE)
        ctab$Hex <- ctab$Colour
        cdict[[d]] <- ctab
    }

    ################################################################
    ## 6. Warn about cross-variable class-name collisions (informational)
    ################################################################
    if (warn_on_collisions){
        long <- do.call(rbind, lapply(names(cdict), function(v){
            data.frame(Variable = v, Name = cdict[[v]]$Name, Colour = cdict[[v]]$Colour, stringsAsFactors = FALSE)
        }))
        colliding_names <- names(which(tapply(long$Colour, long$Name, function(x) length(unique(x)) > 1)))
        colliding_names <- setdiff(colliding_names, class_to_ignore)
        if (length(colliding_names) > 0){
            flog.warn(paste("The following class names appear in more than one variable with differing generated colours and will be collapsed to a single colour in the global ctable (legacy values win, else first-seen wins):", paste0(colliding_names, collapse = ", ")))
        }
    }

    ################################################################
    ## 7. Report on legacy match rate
    ################################################################
    if (!is.null(legacy_lookup)){
        all_used_classes <- unique(unlist(lapply(cdict, function(x) x$Name)))
        matched <- intersect(all_used_classes, names(legacy_lookup))
        if (length(matched) == 0){
            flog.warn("A legacy_ctable was supplied but NONE of its class names matched any class in the metadata being coloured. Check that the legacy file classes and the metadata classes are spelled identically.")
        } else {
            flog.info(paste("Legacy ctable overrode", length(matched), "class colour(s):", paste0(matched, collapse = ", ")))
        }
    }

    ################################################################
    ## 8a. Early return for cdict (per-variable dicts are collision-free by construction)
    ################################################################
    if (object_to_return == "cdict"){
        return(cdict)
    }

    ################################################################
    ## 8b. Collapse to a single global ctable.
    ## Precedence: legacy value if present, otherwise first-seen generated colour.
    ################################################################
    ctable_new <- plyr::rbind.fill(cdict)
    if (!is.null(legacy_lookup)){
        ctable_new$is_legacy <- ctable_new$Name %in% names(legacy_lookup)
        ctable_new <- ctable_new[order(!ctable_new$is_legacy), ]
        ctable_new$is_legacy <- NULL
    }
    ctable_new <- ctable_new[!duplicated(ctable_new$Name), ]
    rownames(ctable_new) <- ctable_new$Name

    ################################################################
    ## 8c. Preventive de-collision pass (per-variable distinctness in the flat ctable)
    ################################################################
    if (avoid_collisions){
        # Per-variable membership of class names (only the coloured discretes).
        var_classes <- lapply(discretes, function(v) unique(as.character(pheno[[v]])))
        names(var_classes) <- discretes

        # Global name -> colour vector.
        gcol <- setNames(ctable_new$Hex, ctable_new$Name)

        # Locked names never move: legacy-backed OR the ignore class(es).
        locked <- setNames(names(gcol) %in% c(names(legacy_lookup), class_to_ignore), names(gcol))

        # How many coloured variables each class appears in (keep the most-shared, move the local).
        var_count <- function(nm){ sum(vapply(discretes, function(v) nm %in% var_classes[[v]], logical(1))) }
        vc <- setNames(vapply(names(gcol), var_count, integer(1)), names(gcol))

        max_iter <- 100
        for (iter in seq_len(max_iter)){
            moved_any <- FALSE
            for (v in discretes){
                cls <- var_classes[[v]]
                cls <- cls[cls %in% names(gcol)]
                if (length(cls) < 2) next
                hexes <- toupper(gcol[cls])
                dup_hex <- unique(hexes[duplicated(hexes)])
                if (length(dup_hex) == 0) next
                for (dh in dup_hex){
                    members <- cls[hexes == dh]
                    lock_members <- members[locked[members]]
                    # Choose keeper: prefer a locked member; else the most-shared; ties -> alphabetical.
                    if (length(lock_members) >= 1){
                        keeper <- lock_members[order(-vc[lock_members], lock_members)][1]
                    } else {
                        keeper <- members[order(-vc[members], members)][1]
                    }
                    movers <- setdiff(members, keeper)
                    movable <- movers[!locked[movers]]
                    unmovable <- movers[locked[movers]]
                    if (length(unmovable) > 0){
                        flog.warn(paste0("avoid_collisions: within variable '", v, "', legacy/ignore-locked class(es) ",
                            paste0(unmovable, collapse = ", "), " share colour ", dh,
                            " with '", keeper, "' and cannot be moved without violating a locked colour; leaving as-is."))
                    }
                    for (m in movable){
                        # Forbidden = colours of m's siblings across every variable m belongs to,
                        # so we never create a new within-variable clash, yet allow intended
                        # cross-variable reuse elsewhere.
                        forbidden <- character(0)
                        for (v2 in discretes){
                            if (m %in% var_classes[[v2]]){
                                sibs <- setdiff(var_classes[[v2]][var_classes[[v2]] %in% names(gcol)], m)
                                forbidden <- c(forbidden, toupper(gcol[sibs]))
                            }
                        }
                        newcol <- nudge_colour(gcol[m], unique(forbidden))
                        if (toupper(newcol) != toupper(gcol[m])){
                            gcol[m] <- newcol
                            moved_any <- TRUE
                        }
                    }
                }
            }
            if (!moved_any) break
        }

        # Residual check (e.g. two locked classes sharing a colour within a variable).
        residual <- character(0)
        for (v in discretes){
            cls <- var_classes[[v]]; cls <- cls[cls %in% names(gcol)]
            if (length(cls) < 2) next
            hexes <- toupper(gcol[cls])
            if (any(duplicated(hexes))) residual <- c(residual, v)
        }
        if (length(residual) > 0){
            flog.warn(paste("avoid_collisions: could not fully resolve within-variable colour clashes for variable(s):",
                paste0(unique(residual), collapse = ", "), "- consider tweak_ctable() for manual adjustment."))
        }

        # Write the de-collided colours back.
        ctable_new$Colour <- unname(gcol[ctable_new$Name])
        ctable_new$Hex <- ctable_new$Colour
    }

    return(ctable_new)
}