#' Constrains a functional (non-taxonomic) JAMS-style SummarizedExperiment object which has feature-by-taxon stratification information (see JAMSbeta -u option) into a SummarizedExperiment object with BaseCounts and GeneCounts only pertaining to the wanted LKTs.

#' @param SEobj a functional (non-taxonomic) JAMS2-style SummarizedExperiment object which has feature-by-taxon stratification information (see JAMSbeta -u option)

#' @param wanted_LKTs vector of LKTs to constrain the functional features of the original SEobj into

#' @param drop_empty_features If TRUE (the default), features which have no bases attributable to any of the wanted LKTs in any sample are dropped from the returned object. If FALSE, every feature of the input object is retained, with zeroes where nothing is attributable.

#' @export

constrain_JAMS_SEobj_to_wanted_taxa <- function(SEobj = NULL, wanted_LKTs = NULL, drop_empty_features = TRUE){

    #Row pointers live in the SparseIndex assay. The JAMS1 long-form allfeaturesbytaxa_index
    #path this function used to rely on was removed in JAMS 2.3.1; JAMS2 objects never had that
    #slot, so on a current object the old code failed at the first line of real work.
    assert_JAMS2_SEobj(SEobj = SEobj, objname = "SEobj", require_stratification = TRUE, calling_function = "constrain_JAMS_SEobj_to_wanted_taxa")

    if (is.null(wanted_LKTs) || (length(wanted_LKTs) < 1)){
        stop("constrain_JAMS_SEobj_to_wanted_taxa: wanted_LKTs must be a non-empty vector of taxon names.")
    }

    original_metadata <- metadata(SEobj)

    #Collapse the wanted taxon columns of a stratification matrix down to one value per stored
    #row, then scatter those values back into a feature-by-sample matrix using SparseIndex.
    #Cells with no stratification entry (NA in SparseIndex) are genuinely zero: the feature was
    #not seen in that sample, so no bases are attributable to any taxon, wanted or otherwise.
    constrain_matrix <- function(sparse_mat = NULL, index_mat = NULL, wanted_cols = NULL, what = NULL){
        rowtotals <- Matrix::rowSums(sparse_mat[ , wanted_cols, drop = FALSE])
        rowkeys <- rownames(sparse_mat)
        if (is.null(rowkeys)){
            #Fall back to positional keys if the matrix was built without rownames.
            rowkeys <- as.character(seq_len(nrow(sparse_mat)))
        }
        names(rowtotals) <- rowkeys

        out <- matrix(0, nrow = nrow(index_mat), ncol = ncol(index_mat), dimnames = dimnames(index_mat))
        populated <- which(!is.na(index_mat))
        if (length(populated) > 0){
            #as.integer() before as.character() is deliberate: if the SparseIndex has been
            #promoted to double anywhere upstream, as.character(100000) is "1e+05", which
            #matches no rowname and would send a perfectly good object into the stop() below.
            lookedup <- rowtotals[as.character(as.integer(index_mat[populated]))]
            if (any(is.na(lookedup))){
                stop(paste0("constrain_JAMS_SEobj_to_wanted_taxa: the SparseIndex assay points at ", length(which(is.na(lookedup))), " row(s) which do not exist in the ", what, " stratification matrix. This object's stratification data is internally inconsistent and cannot be constrained safely."))
            }
            out[populated] <- as.numeric(lookedup)
        }

        return(out)
    }

    index_mat <- as.matrix(assays(SEobj)[["SparseIndex"]])

    sparse_mat <- original_metadata[["allfeaturesbytaxa_matrix"]]
    wanted_present <- intersect(wanted_LKTs, colnames(sparse_mat))
    if (length(wanted_present) < 1){
        stop(paste0("constrain_JAMS_SEobj_to_wanted_taxa: none of the ", length(wanted_LKTs), " taxa in wanted_LKTs are present among the ", ncol(sparse_mat), " taxa carrying features in this object. Check the taxon names, remembering that they are the contextualized names used in the stratification matrix (e.g. \"LKT__s__562_Escherichia_coli\")."))
    }
    if (length(wanted_present) < length(wanted_LKTs)){
        missing_LKTs <- setdiff(wanted_LKTs, wanted_present)
        flog.warn(paste("The following", length(missing_LKTs), "taxa in wanted_LKTs carry no features in this object and will be ignored:", paste0(missing_LKTs, collapse = ", ")))
    }

    cts <- constrain_matrix(sparse_mat = sparse_mat, index_mat = index_mat, wanted_cols = wanted_present, what = "BaseCounts")

    #Second, the GeneCounts, if present. Both matrices share the same row numbering, so the same
    #SparseIndex drives both.
    featcts <- NULL
    if ("allfeaturesbytaxa_GeneCounts_matrix" %in% names(original_metadata)){
        sparse_mat_GC <- original_metadata[["allfeaturesbytaxa_GeneCounts_matrix"]]
        wanted_present_GC <- intersect(wanted_present, colnames(sparse_mat_GC))
        if (length(wanted_present_GC) > 0){
            featcts <- constrain_matrix(sparse_mat = sparse_mat_GC, index_mat = index_mat, wanted_cols = wanted_present_GC, what = "GeneCounts")
        }
    }

    #Drop features with nothing attributable to the wanted taxa anywhere.
    if (drop_empty_features){
        featuresToKeep <- rownames(cts)[rowSums(cts) > 0]
        if (length(featuresToKeep) < 1){
            stop("constrain_JAMS_SEobj_to_wanted_taxa: no feature has any bases attributable to the wanted taxa, so there is nothing to return. Either the taxa chosen carry no annotated genes in this analysis space, or the object has been filtered down too far.")
        }
        if (length(featuresToKeep) < nrow(cts)){
            flog.info(paste("Keeping", length(featuresToKeep), "of", nrow(cts), "features which have bases attributable to the wanted taxa."))
        }
        cts <- cts[featuresToKeep, , drop = FALSE]
        if (!is.null(featcts)){
            featcts <- featcts[featuresToKeep, , drop = FALSE]
        }
    }

    SEobj <- SEobj[rownames(cts), ]

    assays_constrained <- list()
    assays_constrained$BaseCounts <- cts
    if (!is.null(featcts)){
        assays_constrained$GeneCounts <- featcts
    }

    phenoanal <- colData(SEobj)
    ftt <- rowData(SEobj)
    constrained_SEobj <- SummarizedExperiment(assays = assays_constrained, rowData = as.matrix(ftt), colData = as.matrix(phenoanal))

    #Put the original bells and whistles back, except the stratification payload itself. The
    #returned object describes only the wanted taxa, so carrying the full feature-by-taxon
    #matrices (and a SparseIndex pointing into them) would invite exactly the kind of dangling
    #pointer this rewrite exists to eliminate.
    relevant_md <- names(original_metadata)[!(names(original_metadata) %in% c("allfeaturesbytaxa_GeneCounts_index", "allfeaturesbytaxa_GeneCounts_matrix", "allfeaturesbytaxa_index", "allfeaturesbytaxa_matrix"))]
    for (md in relevant_md){
        metadata(constrained_SEobj)[[md]] <- original_metadata[[md]]
    }

    #Record what was done, so a constrained object is never mistaken for a complete one.
    metadata(constrained_SEobj)$constrained_to_taxa <- wanted_present

    #TotalBasesSequenced is deliberately left as the original sample totals, so that PPM figures
    #computed downstream remain a fraction of everything sequenced in the sample rather than a
    #fraction of the chosen taxa, which is almost always the question being asked.
    return(constrained_SEobj)
}
