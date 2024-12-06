#' Get matrices for a specific cell type with specific cell IDs
#'
#' @param cell_type A string for which cell type it should be filtering for.
#' @param sample_list A list of matrices. Ensure all the matrices have the same number of rows (genes). The genes should be the same for all matrices and in the same order. The names of the columns should be the cell types.
#' @param prefix The prefix that your cell IDs should have.
#'
#' @return A sparse filtered matrix.
#' @export
#'
#' @examples
#'
#' set.seed(2024)
#' sample1 <- matrix(rnorm(20), nrow = 5, ncol = 4,
#' dimnames = list(paste0("Gene", 1:5),
#'  c("Ast", "Neu", "Oligo", "Ast")))
#' sample2 <- matrix(rnorm(20), nrow = 5, ncol = 4,
#'  dimnames = list(paste0("Gene", 1:5),
#'   c("Ast", "Ast", "Oligo", "Microglia")))
#' sample_list <- list(sample1, sample2)
#'
#' generate_cell_types("Ast", sample_list, "Cohort_name_")
generate_cell_types <- function(cell_type, sample_list, prefix) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required but is not installed.")
  }

  final_list <- base::list()

  common_genes <- base::rownames(sample_list[[1]])

  for (sample_data in sample_list) {
    cell_counter <- 1

    for (col_idx in base::seq_along(base::colnames(sample_data))){
      column <- base::colnames(sample_data)[col_idx]

      if (column == cell_type) {
        new_column_name <- base::paste0(prefix, cell_type, cell_counter)
        cell_counter <- cell_counter + 1
        final_list[[new_column_name]] <- sample_data[, col_idx]
      }
    }
  }
  final_dataframe <- base::as.data.frame(final_list, check.names = FALSE)
  base::rownames(final_dataframe) <- common_genes
  final_matrix <- base::as.matrix(final_dataframe)
  final_dgCmatrix <- Matrix::Matrix(final_matrix, sparse = TRUE)

  return(final_dgCmatrix)
}
