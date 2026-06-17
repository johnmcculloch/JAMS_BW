#' plot_taxsplit_tree()
#'
#' Helper plotting function for use exclusively within plot_relabund_features (q.v.). Do not attempt to use this outside of this context, or you will find yourself in a frightful muddle.
#' @export

plot_taxsplit_tree <- function(taxsplit_df = NULL,
                               strat_tt = NULL,
                               layout           = c("rectangular", "fan", "circular"),
                               open.angle       = 5,
                               tip_offset       = NULL,
                               fruit_offset     = NULL,
                               fruit_pwidth     = NULL,
                               bar_width        = 0.85,
                               tip_size         = NULL,
                               show_prevalence  = TRUE,
                               bar_fill         = "phylum",
                               low_col          = "#fde0dd",
                               high_col         = "#7a0177",
                               hilight_rank     = "Phylum",
                               taxon_hilight_palette  = NULL,
                               hilight_alpha    = 0.4,
                               tree_size        = 0.4,
                               legend_position  = "right",
                               plotitstrat = NULL) {

  required_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales",
                     "ape", "ggtree", "ggtreeExtra", "ggnewscale")
  invisible(lapply(required_pkgs, function(p) {
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("Package '%s' is required.", p))
  }))

  library(ggplot2); library(dplyr); library(tidyr)
  library(ape); library(ggtree); library(ggtreeExtra)

  layout <- match.arg(layout)

  if (is.null(tip_offset))   tip_offset   <- if (layout == "rectangular") 0.15 else 8
  if (is.null(fruit_offset)) fruit_offset <- if (layout == "rectangular") 3    else 0.15
  if (is.null(tip_size))     tip_size     <- if (layout == "rectangular") 2.5  else 0.9
  if (is.null(fruit_pwidth)) fruit_pwidth <- if (layout == "rectangular") 1    else 0.5

  #For the time being, only accept Compareby
  taxsplit_df <- taxsplit_df[ , !(colnames(taxsplit_df) %in% c("Shape", "Fill", "Connect", "Colour"))]

  meta_cols    <- c("Sample", "Accession", "Compareby")
  taxon_cols   <- setdiff(colnames(taxsplit_df), meta_cols)
  terminal_col <- colnames(strat_tt)[ncol(strat_tt)]

  accession_label <- unique(taxsplit_df[["Accession"]])
  if (length(accession_label) > 1) accession_label <- accession_label[1]

  # ── tree construction ────────────────────────────────────────────────────────
  make_taxonomy_tree <- function(tt) {
    tt <- as.data.frame(tt, stringsAsFactors = FALSE)
    tt[] <- lapply(tt, function(col) {
      col <- as.character(col)
      col[is.na(col) | col == ""] <- "Missing"
      col
    })
    rank_cols <- colnames(tt); n_ranks <- ncol(tt); tips <- tt[[n_ranks]]

    make_node <- function(label = "") {
      n <- new.env(parent = emptyenv())
      n$label <- label; n$children <- list(); n
    }
    add_path <- function(node, path) {
      if (length(path) == 0) return(invisible(NULL))
      key <- path[[1]]
      if (is.null(node$children[[key]])) node$children[[key]] <- make_node(key)
      add_path(node$children[[key]], path[-1])
      invisible(NULL)
    }
    node_to_newick <- function(node, is_root = FALSE) {
      kids <- node$children
      if (length(kids) == 0) return(node$label)
      child_txt <- vapply(kids, node_to_newick, character(1), is_root = FALSE)
      lbl <- if (is_root) "" else node$label
      if (nzchar(lbl)) paste0("(", paste(child_txt, collapse = ","), ")", lbl)
      else             paste0("(", paste(child_txt, collapse = ","), ")")
    }

    root <- make_node("")
    for (i in seq_len(nrow(tt))) {
      path <- character(0)
      for (rk in rank_cols[-n_ranks]) {
        val <- tt[i, rk]
        if (nzchar(val)) path <- c(path, paste0(rk, "_", val))
      }
      path <- c(path, tips[i])
      add_path(root, path)
    }
    ape::read.tree(text = paste0(node_to_newick(root, is_root = TRUE), ";"))
  }

  strat_sub <- strat_tt[strat_tt[[terminal_col]] %in% taxon_cols, , drop = FALSE]
  if (nrow(strat_sub) == 0) stop("No taxa in strat_tt match taxsplit_df columns.")
  phy <- make_taxonomy_tree(strat_sub)

  # ── summarise ────────────────────────────────────────────────────────────────
  long_df <- taxsplit_df %>%
    dplyr::select(all_of(c(meta_cols, taxon_cols))) %>%
    pivot_longer(cols = all_of(taxon_cols), names_to = "Taxon", values_to = "PPM") %>%
    mutate(PPM = replace_na(PPM, 0))

  summary_df <- long_df %>%
    group_by(Compareby, Taxon) %>%
    summarise(
      n_samples    = n(),
      n_nonzero    = sum(PPM > 0),
      prevalence   = n_nonzero / n_samples,
      mean_nonzero = ifelse(n_nonzero > 0, mean(PPM[PPM > 0]), 0),
      .groups      = "drop"
    )

  groups <- sort(unique(summary_df$Compareby))

  prev_summary <- taxsplit_df %>%
    mutate(total_PPM = rowSums(dplyr::select(., all_of(taxon_cols)))) %>%
    group_by(Compareby) %>%
    summarise(n_total = n(), n_nz = sum(total_PPM > 0), .groups = "drop") %>%
    mutate(lbl = sprintf("%s: %d/%d samples", Compareby, n_nz, n_total))
  subtitle_str <- paste(prev_summary$lbl, collapse = "   |   ")

  # ── base tree ────────────────────────────────────────────────────────────────
  if (layout == "fan") {
    p_tree <- ggtree(phy, layout = "fan", open.angle = open.angle) +
      geom_tree(size = tree_size, colour = "grey25")
  } else {
    p_tree <- ggtree(phy, layout = layout) +
      geom_tree(size = tree_size, colour = "grey25")
  }

  # ── phylum highlights ────────────────────────────────────────────────────────
  do_hilight <- !is.null(hilight_rank) && hilight_rank %in% colnames(strat_tt)
  taxon_phylum_col <- NULL

  if (do_hilight) {
    rank_prefix <- paste0("^", hilight_rank, "_")
    hilight_nodes <- p_tree$data %>%
      dplyr::filter(grepl(rank_prefix, label)) %>%
      dplyr::mutate(
        bare_name    = sub(paste0("^", hilight_rank, "_[a-zA-Z]+__[0-9]*_?"), "", label),
        display_name = sub(paste0("^", hilight_rank, "_"), "", label)
      ) %>%
      dplyr::filter(!bare_name %in% c("Missing", "Unclassified", "")) %>%
      dplyr::distinct(bare_name, display_name, node)
    if (nrow(hilight_nodes) > 0) {
      hilight_nodes <- as.data.frame(hilight_nodes)
      rownames(hilight_nodes) <- hilight_nodes$display_name
      n_hl <- nrow(hilight_nodes)
      if (requireNamespace("RColorBrewer", quietly = TRUE) && n_hl <= 8) {
        base_cols <- RColorBrewer::brewer.pal(max(3, n_hl), "Set2")[seq_len(n_hl)]
      } else {
        base_cols <- scales::hue_pal()(n_hl)
      }
      hilight_nodes$hilight_palette <- base_cols
      #Replace palette if required
      if (!is.null(taxon_hilight_palette)) {
        hilight_nodes[intersect(names(taxon_hilight_palette), hilight_nodes$display_name), "hilight_palette"]<- as.data.frame(taxon_hilight_palette)[intersect(names(taxon_hilight_palette), hilight_nodes$display_name), 1]
      }
      hilight_palette <- setNames(hilight_nodes$hilight_palette, hilight_nodes$bare_name)

#      hilight_nodes <- hilight_nodes %>%
#        dplyr::mutate(
#          fill_col     = hilight_palette[bare_name],
#          legend_label = display_name
#        ) %>%
#        dplyr::filter(!is.na(fill_col))

      hilight_nodes$fill_col <- hilight_nodes$hilight_palette
      hilight_nodes$legend_label <- hilight_nodes$display_name

      p_tree <- p_tree +
        ggtree::geom_hilight(
          data        = hilight_nodes,
          mapping     = aes(node = node, fill = legend_label),
          type        = "auto",
          alpha       = hilight_alpha,
          show.legend = TRUE
        ) +
        scale_fill_manual(
          values = setNames(hilight_nodes$fill_col, hilight_nodes$legend_label),
          name   = hilight_rank,
          guide  = guide_legend(override.aes = list(alpha = 1), order = 1)
        )

      # Build per-taxon phylum colour lookup for use in bar fill below
      taxon_to_phylum <- setNames(
        as.character(strat_sub[[hilight_rank]]),
        as.character(strat_sub[[terminal_col]])
      )
      bare_phylum <- sub("^[a-zA-Z]+__[0-9]*_?", "", taxon_to_phylum)
      taxon_phylum_col <- setNames(
        hilight_palette[bare_phylum],
        names(taxon_to_phylum)
      )
      taxon_phylum_col[is.na(taxon_phylum_col)] <- "grey60"
    }
  }

  # ── tip labels ───────────────────────────────────────────────────────────────
  if (layout == "fan") {
    p_tree <- p_tree +
      ggtree::geom_tiplab2(size = tip_size, offset = tip_offset, colour = "grey15")
  } else {
    p_tree <- p_tree +
      geom_tiplab(size = tip_size, align = TRUE, linetype = "dotted",
                  linesize = 0.2, offset = tip_offset, colour = "grey15")
  }

  # ── fruit rings ──────────────────────────────────────────────────────────────
  if (do_hilight) p_tree <- p_tree + ggnewscale::new_scale_fill()

  # Helper: compute fill_hex for a vector of prevalence + taxa
  compute_fill_hex <- function(prev_vec, taxa_vec) {
    taxa_vec <- as.character(taxa_vec)
    n <- length(taxa_vec)
    if (n == 0) return(character(0))
    
    # Treat NULL bar_fill as "phylum" — sensible default for stratification context
    effective_bar_fill <- if (is.null(bar_fill)) "phylum" else bar_fill
    
    out <- if (show_prevalence) {
      ramp <- colorRamp(c(low_col, high_col))
      vapply(prev_vec, function(p) {
        if (is.na(p)) return("#FFFFFF")
        t <- sqrt(pmin(pmax(p, 0), 1))
        rv <- ramp(t)
        rgb(rv[1, 1], rv[1, 2], rv[1, 3], maxColorValue = 255)
      }, character(1))
    } else if (identical(effective_bar_fill, "phylum") && !is.null(taxon_phylum_col)) {
      vapply(taxa_vec, function(tx) {
        col <- taxon_phylum_col[[tx]]
        if (is.null(col) || is.na(col)) "grey60" else as.character(col)
      }, character(1), USE.NAMES = FALSE)
    } else {
      rep(as.character(effective_bar_fill), n)
    }
    
    if (length(out) != n) out <- rep("grey60", n)
    out
  }

  for (i in seq_along(groups)) {
    grp <- groups[i]
    ring_data <- summary_df %>%
      filter(Compareby == grp) %>%
      dplyr::rename(label = Taxon) %>%
      dplyr::mutate(fill_hex = compute_fill_hex(prevalence, label))

    ring_title <- paste0(grp, "\nmean PPM")

    p_tree <- p_tree +
      ggtreeExtra::geom_fruit(
        data        = ring_data,
        geom        = geom_col,
        mapping     = aes(y = label, x = mean_nonzero, fill = fill_hex),
        width       = bar_width,
        pwidth      = fruit_pwidth,
        offset      = if (i == 1) fruit_offset else 0.08,
        color       = NA,
        axis.params = list(
          axis = "x", text.size = 2.6, hjust = 0.5, vjust = 1,
          text.angle = 0, nbreak = 3, line.size = 0.3, line.color = "grey40",
          title = ring_title, title.size = 2, title.height = 0.05
        ),
        grid.params = list(linetype = "dashed", colour = "grey85", size = 0.2)
      )
  }

  # ── single manual fill scale: hex values map to themselves ───────────────────
  all_hex <- unique(unlist(lapply(groups, function(grp) {
    rd <- summary_df %>% filter(Compareby == grp)
    compute_fill_hex(rd$prevalence, rd$Taxon)
  })))
  all_hex <- all_hex[!is.na(all_hex)]
  manual_pal <- setNames(all_hex, all_hex)

  p_tree <- p_tree + scale_fill_manual(values = manual_pal, guide = "none")

  # ── prevalence colourbar legend (decorative, when show_prevalence = TRUE) ────
  if (show_prevalence) {
    p_tree <- p_tree +
      ggnewscale::new_scale_fill() +
      geom_blank(
        data    = data.frame(x = NA_real_, y = NA_real_, prev = c(0, 1)),
        mapping = aes(x = x, y = y, fill = prev),
        inherit.aes = FALSE
      ) +
      scale_fill_gradient(
        low      = low_col,
        high     = high_col,
        limits   = c(0, 1),
        trans    = "sqrt",
        breaks   = c(0, 0.025, 0.1, 0.25, 0.5, 1),
        labels   = scales::label_percent(),
        name     = "Prevalence",
        na.value = "grey95",
        guide    = guide_colourbar(order = 2)
      )
  }

  # ── final labels ─────────────────────────────────────────────────────────────
  fill_descr <- if (show_prevalence) {
    "Fill = prevalence"
  } else if (identical(bar_fill, "phylum")) {
    "Fill = phylum"
  } else {
    sprintf("Fill = %s", bar_fill)
  }

  #Use just the accession if not passing a ready made title
  if (is.null(plotitstrat)){
    plotitstrat <- accession_label
  }

  p_tree <- p_tree +
    labs(
      title    = plotitstrat,
      subtitle = paste0(subtitle_str, "\n",
                        "Rings (inner to outer): ",
                        paste(groups, collapse = " -> "),
                        "  |  Bar length = conditioned mean PPM  |  ",
                        fill_descr)
    ) +
    theme(
      legend.position  = legend_position,
      legend.box       = "vertical",
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(3, "mm"),
      legend.spacing.y = unit(1, "mm"),
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 7, colour = "grey40")
    )

  return(p_tree)
}