#' rbind_sparse_blocks(blocks = NULL, row_offsets = NULL, total_rows = NULL, rownames_out = NULL)
#'
#' Stacks a list of per-sample sparse matrices into ONE tall sparse matrix in a single pass,
#' taking the union of their column names. Each block occupies a contiguous, non-overlapping
#' band of rows given by row_offsets, which is exactly the situation in
#' make_SummarizedExperiments: every sample contributes its own fresh rows and no row is ever
#' shared between samples.
#'
#' WHY THIS EXISTS. The functional phase used to call merge_sparse_matrix() once per sample, so
#' the master matrix was rebuilt from scratch on every iteration. Worse, merge_sparse_matrix
#' widened matrices with Reduce(cbind, ...), which appends ONE column at a time and copies the
#' whole matrix each time. Cost was therefore roughly O(nsamples^2 x nnz x new_columns), which is
#' why 1409 samples projected to ~26 h. This function is O(total number of stored values) and
#' runs once.
#'
#' HOW IT IS CORRECT WITHOUT SORTING. Compressed-column (CSC) storage requires the row indices
#' within each column to be in increasing order. Because block k owns rows
#' (row_offsets[k] + 1) .. (row_offsets[k] + nrow(block k)), and the blocks are supplied in
#' increasing offset order, simply concatenating each column's entries block by block already
#' yields increasing row indices. So the output is assembled directly, with no sort.
#'
#' COLUMN ORDER is the order of first appearance across the blocks, which is identical to the
#' order the old iterative merge produced (merge kept the accumulated columns and appended the
#' new ones). Nothing downstream depends on it - retrieve_features_by_taxa sorts the taxon
#' columns itself - but reproducing it keeps the objects comparable with previously built ones.
#'
#' @param blocks List of sparse matrices (or NULL entries, which are treated as all-zero bands).
#'   Every non-NULL block must have column names.
#' @param row_offsets Integer vector, same length as blocks: the number of rows that precede
#'   each block in the output (so the first block normally has offset 0).
#' @param total_rows Integer. Number of rows of the result. Allowed to exceed the sum of the
#'   block row counts if some bands are deliberately empty.
#' @param rownames_out Optional character vector of length total_rows to use as rownames.
#' @return A dgCMatrix, or NULL if there is nothing to assemble.
#' @export

rbind_sparse_blocks <- function(blocks = NULL, row_offsets = NULL, total_rows = NULL, rownames_out = NULL){

    if (is.null(blocks) || length(blocks) < 1){
        return(NULL)
    }

    usable <- which(!vapply(blocks, is.null, logical(1)))
    if (length(usable) < 1){
        return(NULL)
    }

    if (is.null(row_offsets) || length(row_offsets) != length(blocks)){
        stop("rbind_sparse_blocks: row_offsets must be supplied and must be the same length as blocks.")
    }
    #The no-sort argument below depends on the bands being laid out in increasing order, so say
    #so now rather than emit a silently malformed matrix.
    usable_offsets <- as.numeric(row_offsets[usable])
    if (length(usable_offsets) > 1){
        if (any(diff(usable_offsets) < 0)){
            stop("rbind_sparse_blocks: row_offsets must be non-decreasing across the non-NULL blocks, otherwise the stacked row indices would not be sorted.")
        }
    }

    #Coerce anything that is not already a general, double, column-compressed sparse matrix.
    for (k in usable){
        if (!inherits(blocks[[k]], "dgCMatrix")){
            blocks[[k]] <- as(as(as(blocks[[k]], "dMatrix"), "generalMatrix"), "CsparseMatrix")
        }
        if (is.null(colnames(blocks[[k]]))){
            stop("rbind_sparse_blocks: every block must carry column names.")
        }
        #Duplicate column names within a block cannot be stacked: the two columns would be
        #written through the same output cursor, overrunning their allocation and leaving the
        #row indices unsorted, which silently produces a malformed dgCMatrix. Callers must
        #collapse same-named columns first (build_sample_stratification_block sums them).
        if (anyDuplicated(colnames(blocks[[k]])) > 0){
            stop(paste0("rbind_sparse_blocks: block ", k, " has duplicated column names (", paste0(unique(colnames(blocks[[k]])[duplicated(colnames(blocks[[k]]))]), collapse = ", "), "). Same-named columns must be collapsed before stacking."))
        }
    }

    #Each block must fit entirely inside its own band.
    block_ends <- vapply(usable, function(k) { (as.numeric(row_offsets[k]) + nrow(blocks[[k]])) }, numeric(1))
    if (is.null(total_rows)){
        total_rows <- max(block_ends)
    }
    if (any(block_ends > as.numeric(total_rows))){
        stop("rbind_sparse_blocks: at least one block extends past total_rows. row_offsets and total_rows do not describe the blocks supplied.")
    }
    if (as.numeric(total_rows) > .Machine$integer.max){
        stop(paste0("rbind_sparse_blocks: the stratification matrix would need ", total_rows, " rows, which exceeds the ", .Machine$integer.max, " row limit of a sparse matrix in the Matrix package. Reduce the number of samples or the number of analyses built in one go."))
    }
    if (!is.null(rownames_out) && (length(rownames_out) != total_rows)){
        stop("rbind_sparse_blocks: rownames_out must have exactly total_rows entries.")
    }

    #Global column dictionary, in order of first appearance.
    all_cols <- unique(unlist(lapply(blocks[usable], colnames), use.names = FALSE))
    ncol_total <- length(all_cols)
    if (ncol_total < 1){
        return(NULL)
    }
    col_pos <- seq_len(ncol_total)
    names(col_pos) <- all_cols

    #Pass 1: how many stored values will each output column hold?
    counts <- numeric(ncol_total)
    for (k in usable){
        b <- blocks[[k]]
        if (ncol(b) < 1){
            next
        }
        gj <- col_pos[colnames(b)]
        counts[gj] <- counts[gj] + diff(b@p)
    }

    p_out <- c(0, cumsum(counts))
    total_nnz <- p_out[ncol_total + 1]
    if (total_nnz > .Machine$integer.max){
        stop(paste0("rbind_sparse_blocks: the stratification matrix would hold ", total_nnz, " non-zero values, which exceeds what the Matrix package can store in a single dgCMatrix (", .Machine$integer.max, "). Build fewer analyses or fewer samples per run."))
    }

    i_out <- integer(total_nnz)
    x_out <- numeric(total_nnz)
    #Write cursor per output column, 0 based.
    fill <- p_out[seq_len(ncol_total)]

    #Pass 2: copy each block's columns into their slot in the output.
    for (k in usable){
        b <- blocks[[k]]
        nb <- ncol(b)
        if (nb < 1){
            next
        }
        bp <- b@p
        bi <- b@i
        bx <- b@x
        off <- row_offsets[k]
        gj <- col_pos[colnames(b)]
        for (lj in seq_len(nb)){
            n <- bp[lj + 1] - bp[lj]
            if (n == 0){
                next
            }
            g <- gj[lj]
            src <- (bp[lj] + 1):(bp[lj + 1])
            dst <- (fill[g] + 1):(fill[g] + n)
            #Both bi and i_out are 0 based, which is what the dgCMatrix i slot stores.
            i_out[dst] <- bi[src] + off
            x_out[dst] <- bx[src]
            fill[g] <- fill[g] + n
        }
    }

    #Assembled directly rather than through Matrix::sparseMatrix() because every CSC invariant is
    #already guaranteed above: i is 0 based and ascending within each column (blocks are supplied
    #in increasing row-offset order), and p is the cumulative count of stored values per column.
    #This also makes the all-zero case (no stored values at all) unambiguous.
    out <- new("dgCMatrix", i = as.integer(i_out), p = as.integer(p_out), x = as.numeric(x_out),
        Dim = c(as.integer(total_rows), as.integer(ncol_total)),
        Dimnames = list(rownames_out, all_cols))

    return(out)
}
