#' Safely merges two compatible JAMS-style SummarizedExperiment objects into a single JAMS-style SummarizedExperiment object. SE objects must belong to the same analysis type, and have no Sample names in common. Ideally they should have been generated with the same version of JAMS using the same JAMS database. If the latter condition is not met, the consequences are unforseen and you are ill-advised to proceed.

#' @param ExpObj1 a JAMS-style SummarizedExperiment object

#' @param ExpObj2 a JAMS-style SummarizedExperiment object

#' @export

merge_JAMS_SEobj <- function(ExpObj1 = NULL, ExpObj2 = NULL){

    #Both objects must be JAMS2. Merging used to drop metadata$version, which silently
    #reclassified the merged object as JAMS1 and made filter_experiment divide
    #GenomeCompletenessCutoff by 100 - a wrong answer rather than a crash. The stamp is now
    #required on input and deliberately carried through to the output below.
    assert_JAMS2_SEobj(SEobj = ExpObj1, objname = "ExpObj1", calling_function = "merge_JAMS_SEobj")
    assert_JAMS2_SEobj(SEobj = ExpObj2, objname = "ExpObj2", calling_function = "merge_JAMS_SEobj")

    #Determine the characteristics of ExpObj
    characterize_ExpObj <- function(ExpObj = NULL){
        ExpObj_facts <- list()
        ExpObj_facts$class <- class(ExpObj)[1]
        ExpObj_facts$analysis <- metadata(ExpObj)$analysis
        ExpObj_facts$assays <- sort(names(assays(ExpObj)))
        ExpObj_facts$metadata <- sort(names(metadata(ExpObj)))
        ExpObj_facts$dimensions <- dim(ExpObj)

        return(ExpObj_facts)
    }

    ExpObj1_facts <- characterize_ExpObj(ExpObj1)
    ExpObj2_facts <- characterize_ExpObj(ExpObj2)
    Samples_in_common <- intersect(colnames(ExpObj1), colnames(ExpObj2))

    #Check we're merging the same kind of ExpObj. This used to compare ExpObj1 with itself, so
    #the test could never fire and objects from different analysis spaces merged happily.
    if (ExpObj1_facts$analysis != ExpObj2_facts$analysis){
        stop(flog.warn(paste("ExpObj1 is a", ExpObj1_facts$analysis, "SummarizedExperiment object, whereas ExpObj2 is a", ExpObj2_facts$analysis, "SummarizedExperiment object, and therefore cannot be merged. Please use same analysis kind as input for merging SummarizedExperiment objects.")))
    }

    #Check we're not merging ExpObjs with the same samples
    if (length(Samples_in_common) > 0){
        stop(flog.warn(paste("ExpObj1 and ExpObj2 have samples in common. Please use SummarizedExperiment objects with different sample names for merging.\nSamples in common are:", paste0(Samples_in_common, collapse = ", "))))
    }

    if (!identical(as.character(metadata(ExpObj1)$version), as.character(metadata(ExpObj2)$version))){
        flog.warn(paste("ExpObj1 was built with", paste0(as.character(metadata(ExpObj1)$version), collapse = " "), "but ExpObj2 was built with", paste0(as.character(metadata(ExpObj2)$version), collapse = " "), "- merging anyway, but be aware that annotation may not be strictly comparable between them. The merged object will be stamped with the version of ExpObj1."))
    }

    #SparseIndex is handled separately below: its cells are row pointers into the stratification
    #matrices, not counts, so they must be renumbered rather than merged, and its empty cells
    #must stay NA rather than becoming zero.
    assays_to_merge <- setdiff(intersect(ExpObj1_facts$assays, ExpObj2_facts$assays), "SparseIndex")

    #Merge count matrices
    merged_assays <- list()
    for (ass in assays_to_merge){
        flog.info(paste("Merging", ass))
        #optional = TRUE stops as.data.frame() from running make.names() over the sample names.
        #Without it a sample called "P-01" becomes "P.01" and one called "001_Ctrl" becomes
        #"X001_Ctrl", after which sample_order below no longer matches the real colnames of
        #either input object - which silently produces all-NA phenotable rows and, since the
        #stratification merge indexes by sample name, a hard subscript error.
        cts1 <- as.data.frame(assays(ExpObj1)[[ass]], optional = TRUE)
        cts2 <- as.data.frame(assays(ExpObj2)[[ass]], optional = TRUE)
        cts1$FeAtUrEnAmE <- rownames(cts1)
        cts2$FeAtUrEnAmE <- rownames(cts2)
        merged_cts <- base::merge(cts1, cts2, by = "FeAtUrEnAmE", all = TRUE)
        rownames(merged_cts) <- merged_cts$FeAtUrEnAmE
        merged_cts$FeAtUrEnAmE <- NULL
        merged_cts[is.na(merged_cts)] <- 0
        merged_cts <- as.matrix(merged_cts)
        merged_assays[[ass]] <- merged_cts
    }

    if (!("BaseCounts" %in% names(merged_assays))){
        stop("merge_JAMS_SEobj: both objects must have a BaseCounts assay in common.")
    }

    feat_order <- rownames(merged_assays[["BaseCounts"]])[order(rowSums(merged_assays[["BaseCounts"]]), decreasing = TRUE)]
    for (ass in assays_to_merge){
        merged_assays[[ass]] <- merged_assays[[ass]][feat_order, ]
    }
    sample_order <- colnames(merged_assays[["BaseCounts"]])

    #Cheap invariant: everything downstream indexes by sample name, so catch any mangling here
    #rather than three screens later as an inscrutable subscript error.
    if (!setequal(sample_order, c(colnames(ExpObj1), colnames(ExpObj2)))){
        stop("merge_JAMS_SEobj: the sample names were altered while merging the assays, most likely by being coerced into syntactically valid R names. The merge cannot continue safely.")
    }

    #Merge feature table(s)
    flog.info("Merging feature tables")
    ftt1 <- as.data.frame(rowData(ExpObj1))
    ftt2 <- as.data.frame(rowData(ExpObj2))
    ftt1$FeAtUrEnAmE <- rownames(ftt1)
    ftt2$FeAtUrEnAmE <- rownames(ftt2)
    #Not optimal, because there could be annotation conflict, but the user is responsible for merging compatible SE objects, so will keep it like this for the time being.
    common_cols <- intersect(colnames(ftt1), colnames(ftt2))
    ftt1 <- ftt1[ , common_cols, drop = FALSE]
    ftt2 <- ftt2[ , common_cols, drop = FALSE]
    merged_ftt <- rbind(ftt1, ftt2)
    merged_ftt <- merged_ftt[!duplicated(merged_ftt), ]
    #This should be the non-redundant feature table, but check and warn if there are discrepancies
    if (is.redundant(merged_ftt$FeAtUrEnAmE)){
        numdupes <- length(which(duplicated(merged_ftt$FeAtUrEnAmE) == TRUE))
        flog.warn(paste("ATTENTION: a total of", numdupes, "out of", nrow(merged_ftt), "unique feature table rows are divergent. This can be due to different annotation of features between the merging SummarizedExperiment objects. Making the feature names unique."))
        merged_ftt <- merged_ftt[!duplicated(merged_ftt$FeAtUrEnAmE), ]
    }
    rownames(merged_ftt) <- merged_ftt$FeAtUrEnAmE
    merged_ftt$FeAtUrEnAmE <- NULL
    #Should be the same, but enforce anyway
    merged_ftt <- merged_ftt[(rownames(merged_ftt) %in% feat_order), ]
    merged_ftt <- merged_ftt[feat_order, ]

    #Merge the metadata
    flog.info("Merging phenotable")
    pheno1 <- as.data.frame(colData(ExpObj1))
    pheno2 <- as.data.frame(colData(ExpObj2))
    pheno_merged <- base::merge(pheno1, pheno2, all = TRUE)
    rownames(pheno_merged) <- pheno_merged$Sample
    #Fill in the gaps with JAMS-style NA
    for (colm in 1:ncol(pheno_merged)){
        pheno_merged[is.na(pheno_merged[ , colm]), colm] <- "N_A"
    }
    pheno_merged <- pheno_merged[sample_order, ]

    #Merge the remaning bells and whistles (SEobj "metadata")
    SEobjMD_merged <- list()
    metadata_in_common <- intersect(ExpObj1_facts$metadata, ExpObj2_facts$metadata)
    stratification_slots <- c("allfeaturesbytaxa_index", "allfeaturesbytaxa_GeneCounts_index", "allfeaturesbytaxa_matrix", "allfeaturesbytaxa_GeneCounts_matrix")

    for (SEobjMD in metadata_in_common[!(metadata_in_common %in% stratification_slots)]){
        flog.info(paste("Merging", SEobjMD))
        if (SEobjMD == "analysis"){
            SEobjMD_merged[[SEobjMD]] <- ExpObj1_facts$analysis
        } else if (SEobjMD %in% c("TotalBasesSequenced", "TotalBasesSequencedinAnalysis")) {
            SEobjMD_merged[[SEobjMD]] <- cbind(metadata(ExpObj1)[[SEobjMD]], metadata(ExpObj2)[[SEobjMD]])
            #Ensure order, because, again, I am neurotic
            SEobjMD_merged[[SEobjMD]] <- t(as.matrix(SEobjMD_merged[[SEobjMD]]["NumBases", sample_order]))
            rownames(SEobjMD_merged[[SEobjMD]]) <- "NumBases"
        } else if (SEobjMD == "ctable") {
            SEobjMD_merged[[SEobjMD]] <- base::merge(metadata(ExpObj1)[[SEobjMD]], metadata(ExpObj2)[[SEobjMD]], all = TRUE)
            SEobjMD_merged[[SEobjMD]] <- SEobjMD_merged[[SEobjMD]][!duplicated(SEobjMD_merged[[SEobjMD]]$Name), ]
            rownames(SEobjMD_merged[[SEobjMD]]) <- SEobjMD_merged[[SEobjMD]]$Name
        } else if (SEobjMD == "tt") {
            #Taxonomy table used by retrieve_features_by_taxa to agglomerate taxon columns above
            #LKT. Union the two, keeping the first lineage seen for any repeated LKT.
            #bind_rows rather than rbind, because two objects built against different JAMS
            #databases can legitimately have different taxonomy table columns, which is exactly
            #the situation the version warning above contemplates. rbind would error there.
            tt_merged <- dplyr::bind_rows(as.data.frame(metadata(ExpObj1)$tt), as.data.frame(metadata(ExpObj2)$tt))
            tt_merged <- tt_merged[!duplicated(tt_merged$LKT), ]
            rownames(tt_merged) <- tt_merged$LKT
            SEobjMD_merged[[SEobjMD]] <- tt_merged
        } else {
            #Everything else (version, phenolabels, and any slot added in future) is carried
            #through from ExpObj1 rather than dropped. The old code assigned NULL here, which in
            #R removes the element instead of creating it, so metadata$version in particular was
            #silently lost on every merge.
            SEobjMD_merged[[SEobjMD]] <- metadata(ExpObj1)[[SEobjMD]]
        }
    }

    #Guarantee the generation stamp survives even if it was somehow absent from the intersection.
    SEobjMD_merged[["version"]] <- metadata(ExpObj1)$version

    #If these are functional SEobjs containing taxonomy stratification data, merge that.
    strat1 <- ("SparseIndex" %in% ExpObj1_facts$assays) && ("allfeaturesbytaxa_matrix" %in% ExpObj1_facts$metadata)
    strat2 <- ("SparseIndex" %in% ExpObj2_facts$assays) && ("allfeaturesbytaxa_matrix" %in% ExpObj2_facts$metadata)

    merged_SparseIndex <- NULL

    if (strat1 && strat2){
        flog.info("Merging taxonomy stratification data")

        #Compact one object's stratification down to the rows its SparseIndex actually points at,
        #and return both the compacted matrices and the old-to-new row number mapping. Rows are
        #renumbered from scratch so that the two objects occupy contiguous, non-overlapping bands
        #in the merged matrices, which is exactly the layout rbind_sparse_blocks requires.
        compact_stratification <- function(ExpObj = NULL, row_offset = 0, objname = NULL){
            idx <- as.matrix(assays(ExpObj)[["SparseIndex"]])
            used <- sort(unique(as.integer(idx[!is.na(idx)])))
            if (length(used) < 1){
                return(NULL)
            }

            out <- list()
            out$newnumbers <- setNames(as.integer(seq_along(used) + row_offset), as.character(used))

            #Rebuild this object's index with the new global row numbers. as.integer() before
            #as.character() is deliberate: if the SparseIndex has been promoted to double
            #anywhere upstream, as.character(100000) is "1e+05", which matches no key at all.
            newidx <- matrix(NA_integer_, nrow = nrow(idx), ncol = ncol(idx), dimnames = dimnames(idx))
            populated <- which(!is.na(idx))
            newidx[populated] <- out$newnumbers[as.character(as.integer(idx[populated]))]
            out$index <- newidx

            for (matname in c("allfeaturesbytaxa_matrix", "allfeaturesbytaxa_GeneCounts_matrix")){
                thismat <- metadata(ExpObj)[[matname]]
                if (is.null(thismat)){
                    out[[matname]] <- NULL
                    next
                }
                matkeys <- rownames(thismat)
                if (is.null(matkeys)){
                    matkeys <- as.character(seq_len(nrow(thismat)))
                    rownames(thismat) <- matkeys
                }
                if (!all(as.character(used) %in% matkeys)){
                    nmissing <- length(which(!(as.character(used) %in% matkeys)))
                    stop(paste0("merge_JAMS_SEobj: the SparseIndex assay of ", objname, " points at ", nmissing, " row(s) absent from its ", matname, ". That object's stratification data is internally inconsistent and cannot be merged safely."))
                }
                compacted <- thismat[as.character(used), , drop = FALSE]
                rownames(compacted) <- as.character(seq_along(used) + row_offset)
                out[[matname]] <- compacted
            }

            out$nrows <- length(used)

            return(out)
        }

        strat_1 <- compact_stratification(ExpObj = ExpObj1, row_offset = 0, objname = "ExpObj1")
        n1 <- if (is.null(strat_1)) 0 else strat_1$nrows
        strat_2 <- compact_stratification(ExpObj = ExpObj2, row_offset = n1, objname = "ExpObj2")
        n2 <- if (is.null(strat_2)) 0 else strat_2$nrows
        total_rows <- n1 + n2

        if (total_rows > 0){
            #Stack the matrices. rbind_sparse_blocks takes the union of the taxon columns and
            #assembles the result in one pass, which also keeps the CSC row indices sorted
            #because the bands are supplied in increasing row-offset order. It tolerates a NULL
            #block, which is what an object whose SparseIndex is entirely NA contributes.
            for (matname in c("allfeaturesbytaxa_matrix", "allfeaturesbytaxa_GeneCounts_matrix")){
                #Key the decision on whether each object HAS the slot, not on whether it happened
                #to contribute rows. Conflating the two meant that an object with an all-NA
                #SparseIndex looked like a missing slot, so both matrices were skipped while the
                #merged SparseIndex was still built and attached - an index pointing at nothing.
                present1 <- matname %in% ExpObj1_facts$metadata
                present2 <- matname %in% ExpObj2_facts$metadata
                if (!present1 && !present2){
                    next
                }
                if (xor(present1, present2)){
                    flog.warn(paste("Only one of the two objects carries", matname, "- it cannot be merged meaningfully and will be omitted from the merged object."))
                    next
                }
                block1 <- if (is.null(strat_1)) NULL else strat_1[[matname]]
                block2 <- if (is.null(strat_2)) NULL else strat_2[[matname]]
                SEobjMD_merged[[matname]] <- rbind_sparse_blocks(blocks = list(block1, block2), row_offsets = c(0, n1), total_rows = total_rows, rownames_out = as.character(seq_len(total_rows)))
            }

            #Only attach a SparseIndex if the matrix it points into was actually assembled.
            #rbind_sparse_blocks also returns NULL when there are no taxon columns at all.
            if (is.null(SEobjMD_merged[["allfeaturesbytaxa_matrix"]])){
                flog.warn("The feature-by-taxon stratification matrices could not be assembled from these two objects, so the SparseIndex assay is being dropped rather than left pointing at nothing.")
                SEobjMD_merged[["allfeaturesbytaxa_GeneCounts_matrix"]] <- NULL
            } else {
                #Scatter both renumbered indices into one feature-by-sample integer matrix. NA
                #means "this feature was not seen in this sample", which is materially different
                #from zero and must not be filled in. The two objects have disjoint samples, so
                #the writes below never overlap.
                merged_SparseIndex <- matrix(NA_integer_, nrow = length(feat_order), ncol = length(sample_order), dimnames = list(feat_order, sample_order))
                for (st in list(strat_1, strat_2)){
                    if (is.null(st)){
                        next
                    }
                    shared_feats <- intersect(rownames(st$index), feat_order)
                    merged_SparseIndex[shared_feats, colnames(st$index)] <- st$index[shared_feats, , drop = FALSE]
                }
                #Integer end to end: a double here would render large row numbers in scientific
                #notation when looked up with as.character(), which would never match a rowname.
                storage.mode(merged_SparseIndex) <- "integer"
            }
        }

    } else if (strat1 || strat2){
        flog.warn("Only one of the two SummarizedExperiment objects carries feature-by-taxon stratification data (the JAMSbeta -u payload). A half-stratified object would report zero attributable bases for the samples that were never stratified, which is indistinguishable from a real absence, so the stratification is being dropped entirely from the merged object. Rebuild both objects with -u if you need stratification after merging.")
    } else {
        flog.info("No feature stratification by taxonomy to merge.")
    }

    #Build the merged SummarizedExperiment object
    merged_SEobj <- SummarizedExperiment(assays = merged_assays, rowData = as.matrix(merged_ftt), colData = as.matrix(pheno_merged), metadata = SEobjMD_merged)

    if (!is.null(merged_SparseIndex)){
        assays(merged_SEobj)[["SparseIndex"]] <- merged_SparseIndex
    }

    return(merged_SEobj)
}
