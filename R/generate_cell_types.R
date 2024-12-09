#' Generate cell type specific matrices with unique column names
#'
#' @param cell_types Vector of cell types to process (e.g., c("T_cells", "B_cells")).
#' @param sample_list List of matrices where rows are genes and columns are cells. Each matrix should have column names indicating cell types. Genes should be the same (same number and order) for all matrices.
#' @param prefix Prefix to include in the column names of the output.
#'
#' @return A list of sparse matrices (one for each cell type in cell_types), with unique column names and shared row names.
#' @export
#'
#' @examples
#'
#'sample1 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3,
#'dimnames = list(c("Gene1", "Gene2", "Gene3"), c("T_cells", "B_cells")))
#'sample2 <- matrix(c(7, 8, 9, 10, 11, 12), nrow = 3,
#'dimnames = list(c("Gene1", "Gene2", "Gene3"), c("T_cells", "B_cells")))
#'sample_list <- list(sample1, sample2)
#'
#'cell_types <- c("T_cells", "B_cells")
#'
#'prefix <- "Sample_"
#'
#'generate_cell_types(cell_types, sample_list, prefix)
generate_cell_types <- function(cell_types, sample_list, prefix) {
  if (!base::is.character(cell_types) || base::length(cell_types) == 0) {
    stop("Error: 'cell_types' must be a non-empty character vector.")
  }
  if (!base::is.list(sample_list) || base::length(sample_list) == 0) {
    stop("Error: 'sample_list' must be a non-empty list.")
  }
  if (!base::all(base::sapply(sample_list, base::is.matrix))) {
    stop("Error: All elements of 'sample_list' must be matrices.")
  }
  if (!base::is.character(prefix) || base::nchar(prefix) == 0) {
    stop("Error: 'prefix' must be a non-empty string.")
  }

  generate_individual <- function(ct, sample_list, prefix){

    final_list <- base::list()

    common_genes <- base::rownames(sample_list[[1]])
    if (base::is.null(common_genes)) {
      stop("Error: Samples must have rownames representing genes.")
    }
    if (!base::all(base::sapply(sample_list, function(x) base::identical(base::rownames(x), common_genes)))) {
      stop("Error: All samples must have the same rownames (genes).")
    }

    cell_counter <- 1
    for (sample_idx in seq_along(sample_list)) {
      sample_data <- sample_list[[sample_idx]]
      matching_cols <- base::which(base::colnames(sample_data) == ct)

      for (col_idx in matching_cols) {
        new_column_name <- base::paste0(prefix, ct, "_Sample", sample_idx, "_", cell_counter)
        cell_counter <- cell_counter + 1
        final_list[[new_column_name]] <- sample_data[, col_idx]
      }
    }
    final_dataframe <- base::as.data.frame(final_list, check.names = FALSE)
    base::rownames(final_dataframe) <- common_genes
    final_matrix <- base::as.matrix(final_dataframe)
    final_dgCmatrix <- Matrix::Matrix(final_matrix, sparse = TRUE)

    return(final_dgCmatrix)
  }

  ct_list <- base::list()
  for (celltype in cell_types){
    dgcmatrix <- generate_individual(celltype, sample_list, prefix)
    ct_list[[celltype]] <- dgcmatrix
  }

  return(ct_list)
}



