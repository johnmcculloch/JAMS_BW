#' merge_sparse_matrix(matlist = NULL)
#'
#' Safely merges sparse matrices by taking the union of their row names and the union of their
#' column names, placing every stored value at its (rowname, colname) address. Used by
#' merge_JAMS_SEobj to join the taxon-stratification matrices of two JAMS SummarizedExperiment
#' objects.
#'
#' A (row, column) address may be occupied by at most one input matrix. If two inputs both store
#' a non-zero at the same address the merge is ambiguous, so the function warns and returns NULL,
#' which is the behaviour callers have always relied on.
#'
#' HISTORY / WHAT CHANGED AND WHY. The original implementation widened each matrix with
#' Reduce(cbind, ...), which appends a single column per call and copies the entire matrix every
#' time, then permuted the columns and rbind()ed. For the per-sample accumulation in
#' make_SummarizedExperiments that was the dominant cost of the whole JAMSbeta run (~26 h
#' projected for 1409 samples). It also had two defects in the shared-row branch, which neither
#' caller ever exercised because both guarantee disjoint row names:
#'   1. conflicts were detected per COLUMN rather than per CELL, and
#'   2. when rows were shared, values held by the accumulated matrix in columns that both
#'      matrices had in common were silently dropped.
#' This version is a single triplet-union pass: linear in the number of stored values, with
#' exact per-cell conflict detection and no data loss. For the disjoint-row case that both JAMS
#' callers produce, the result is identical to the old one, including row and column order.
#'
#' n.b. make_SummarizedExperiments no longer calls this in a loop; it uses rbind_sparse_blocks
#' (see) to stack all per-sample blocks in one pass.
#'
#' @param matlist A list of sparse matrices, each with both row names and column names.
#' @return A dgCMatrix, or NULL if the matrices cannot be merged unambiguously.
#' @export

merge_sparse_matrix <- function(matlist = NULL) {

    if (is.null(matlist)){
        return(NULL)
    }

    matlist <- matlist[!vapply(matlist, is.null, logical(1))]
    if (length(matlist) < 1){
        return(NULL)
    }

    #Coerce anything that is not already a general, double, column-compressed sparse matrix.
    #This runs BEFORE the single-matrix shortcut, so a lone dgTMatrix or a lone unnamed matrix
    #is still validated and normalised rather than handed straight back.
    for (k in seq_along(matlist)){
        if (!inherits(matlist[[k]], "dgCMatrix")){
            matlist[[k]] <- as(as(as(matlist[[k]], "dMatrix"), "generalMatrix"), "CsparseMatrix")
        }
        if (is.null(rownames(matlist[[k]])) || is.null(colnames(matlist[[k]]))){
            flog.warn("merge_sparse_matrix: every matrix must carry both row names and column names. Returning NULL.")
            return(NULL)
        }
        #Duplicated names within one matrix would be collapsed onto a single row or column by
        #the match() below, silently summing two distinct entities together. Refuse instead.
        if (anyDuplicated(rownames(matlist[[k]])) > 0 || anyDuplicated(colnames(matlist[[k]])) > 0){
            flog.warn(paste0("merge_sparse_matrix: matrix ", k, " has duplicated row or column names, which cannot be merged unambiguously. Returning NULL."))
            return(NULL)
        }
    }

    if (length(matlist) == 1){
        return(matlist[[1]])
    }

    all_rows <- unique(unlist(lapply(matlist, rownames), use.names = FALSE))
    all_cols <- unique(unlist(lapply(matlist, colnames), use.names = FALSE))
    nr <- length(all_rows)
    nc <- length(all_cols)

    if (nr > .Machine$integer.max || nc > .Machine$integer.max){
        flog.warn("merge_sparse_matrix: the merged matrix would exceed the dimension limits of the Matrix package. Returning NULL.")
        return(NULL)
    }

    I_list <- vector("list", length(matlist))
    J_list <- vector("list", length(matlist))
    X_list <- vector("list", length(matlist))

    for (k in seq_along(matlist)){
        m <- matlist[[k]]
        nnz_per_col <- diff(m@p)
        I_list[[k]] <- match(rownames(m), all_rows)[m@i + 1L]
        J_list[[k]] <- match(colnames(m), all_cols)[rep.int(seq_len(ncol(m)), nnz_per_col)]
        X_list[[k]] <- m@x
        m <- NULL
    }

    merged_i <- unlist(I_list, use.names = FALSE)
    merged_j <- unlist(J_list, use.names = FALSE)
    merged_x <- unlist(X_list, use.names = FALSE)
    I_list <- NULL
    J_list <- NULL
    X_list <- NULL

    if (is.null(merged_i)){
        merged_i <- integer(0)
        merged_j <- integer(0)
        merged_x <- numeric(0)
    }

    if (length(merged_i) > .Machine$integer.max){
        flog.warn("merge_sparse_matrix: the merged matrix would hold more non-zero values than the Matrix package can store. Returning NULL.")
        return(NULL)
    }

    #Per-cell conflict detection. Explicit zeros cannot conflict with anything, so they are
    #excluded; this is the same spirit as the original test, but exact.
    occupied <- which(merged_x != 0)
    if (length(occupied) > 1){
        #Cell address as a double, which is exact well beyond 2^31 cells.
        cell_id <- (as.numeric(merged_j[occupied]) - 1) * as.numeric(nr) + as.numeric(merged_i[occupied])
        if (anyDuplicated(cell_id) > 0){
            flog.warn(paste("Impossible to merge matrices as there are non unique row-column pairs."))
            return(NULL)
        }
        cell_id <- NULL
    }

    merged <- NULL
    if (length(merged_i) > 0){
        merged <- Matrix::sparseMatrix(i = as.integer(merged_i), j = as.integer(merged_j),
            x = as.numeric(merged_x),
            dims = c(as.integer(nr), as.integer(nc)),
            dimnames = list(all_rows, all_cols),
            index1 = TRUE, use.last.ij = FALSE)
    } else {
        #Every input was all-zero. Hand back an explicitly empty matrix of the right shape.
        merged <- new("dgCMatrix", i = integer(0), p = integer(nc + 1), x = numeric(0),
            Dim = c(as.integer(nr), as.integer(nc)),
            Dimnames = list(all_rows, all_cols))
    }

    return(merged)
}
