#' Safely merges two sparse matrices. This function was written for exclusive use by the merge_JAMS_SEobj function (see). It has not been tested or benchmarked for a use other than that one.

#' @param matlist - A list of sparse matrices

#' @export

merge_sparse_matrix <- function(matlist = NULL) {

    root_mat <- matlist[[1]]
    for (matnum in 2:length(matlist)){
        curr_mat <- matlist[[matnum]]

        Not_in_root_cols <- colnames(curr_mat)[!colnames(curr_mat) %in% colnames(root_mat)]
        Not_in_root_cols_list <- as.vector(numeric(length(Not_in_root_cols)), "list")
        names(Not_in_root_cols_list) <- Not_in_root_cols
        Present_in_root_cols <- colnames(curr_mat)[colnames(curr_mat) %in% colnames(root_mat)]

        Not_in_curr_cols <- colnames(root_mat)[!colnames(root_mat) %in% colnames(curr_mat)]
        Not_in_curr_cols_list <- as.vector(numeric(length(Not_in_curr_cols)), "list")
        names(Not_in_curr_cols_list) <- Not_in_curr_cols
        Present_in_curr_cols <- colnames(root_mat)[colnames(root_mat) %in% colnames(curr_mat)]

        Not_in_root_rows <- rownames(curr_mat)[!rownames(curr_mat) %in% rownames(root_mat)]
        Not_in_root_rows_list <- as.vector(numeric(length(Not_in_root_rows)), "list")
        names(Not_in_root_rows_list) <- Not_in_root_rows
        Present_in_root_rows <- rownames(curr_mat)[rownames(curr_mat) %in% rownames(root_mat)]

        Not_in_curr_rows <- rownames(root_mat)[!rownames(root_mat) %in% rownames(curr_mat)]
        Not_in_curr_rows_list <- as.vector(numeric(length(Not_in_curr_rows)), "list")
        names(Not_in_curr_rows_list) <- Not_in_curr_rows
        Present_in_curr_rows <- rownames(root_mat)[rownames(root_mat) %in% rownames(curr_mat)]

        root_expanded <- Reduce(cbind, c(root_mat, Not_in_root_cols_list))
        if (length(Not_in_root_cols) > 0){
            Startpos_root <- (ncol(root_expanded) - length(Not_in_root_cols)) + 1
            colnames(root_expanded)[Startpos_root:ncol(root_expanded)] <- names(Not_in_root_cols_list)
        }

        curr_expanded <- Reduce(cbind, c(curr_mat, Not_in_curr_cols_list))
        if (length(Not_in_curr_cols) > 0){
            Startpos_curr <- (ncol(curr_expanded) - length(Not_in_curr_cols)) + 1
            colnames(curr_expanded)[Startpos_curr:ncol(curr_expanded)] <- names(Not_in_curr_cols_list)
        }

        #Add row values for shared rownames if necessary
        if (length(Present_in_curr_rows) > 0){
            #Test to see if there are no conflicting values
            Non_sparse_cols_curr <- colnames(curr_expanded)[diff(curr_expanded[Present_in_curr_rows,]@p) > 0]
            Non_sparse_cols_root <- colnames(root_mat)[diff(root_mat[Present_in_curr_rows,]@p) > 0]

            if (length(intersect(Non_sparse_cols_curr, Non_sparse_cols_root)) > 0){
                flog.warn(paste("Impossible to merge matrices as there are non unique row-column pairs."))
                return(NULL)
            } else {
                #Add values from root to curr_expanded
                curr_expanded[Present_in_curr_rows, Not_in_curr_cols] <- root_mat[Present_in_curr_rows, Not_in_curr_cols]
            }

            #Remove these rows from root_expanded
            root_expanded <- root_expanded[!(rownames(root_expanded) %in% Present_in_curr_rows), ]
        }

        curr_expanded <- curr_expanded[ , colnames(root_expanded)]

        #Merge and call it the root_mat to iteratively expand
        root_mat <- rbind(root_expanded, curr_expanded)
    }

    return(root_mat)

}
