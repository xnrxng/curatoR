#' Generate cell type specific matrices with barcodes
#'
#' @param cell_types Vector of cell types to process (e.g., c("T_cells", "B_cells")).
#' @param sample_list List of matrices where rows are genes and columns are cells. Each matrix should have column names indicating cell types. Genes should be the same (same number and order) for all matrices.
#' @param cell_IDs A dataframe that should have two columns: "cell_type" and "cell_ID".
#'
#' @return A list of sparse matrices (one for each cell type in cell_types), with cell IDs as column names and shared row names.
#' @export
#'
#' @examples
#'
#' cell_IDs <- data.frame(
#' cell_type = c("T_cells", "B_cells", "T_cells", "B_cells",
#' "T_cells", "T_cells", "B_cells", "B_cells"),
#' cell_ID = c("Cell1", "Cell2", "Cell3", "Cell4",
#' "Cell5", "Cell6", "Cell7", "Cell8"))
#'
#' sample1 <- matrix(
#' c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
#' nrow = 3,
#' dimnames = list(c("Gene1", "Gene2", "Gene3"),
#' c("Cell1", "Cell2", "Cell3", "Cell4")))
#'
#' sample2 <- matrix(
#' c(13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24),
#' nrow = 3,
#' dimnames = list(c("Gene1", "Gene2", "Gene3"),
#' c("Cell5", "Cell6", "Cell7", "Cell8")))
#' sample_list <- list(sample1, sample2)
#' cell_types <- c("T_cells", "B_cells")
#' generate_cell_types_w_barcodes(cell_types, sample_list, cell_IDs)
generate_cell_types_w_barcodes <- function(cell_types, sample_list, cell_IDs){
  if (!base::is.character(cell_types) || base::length(cell_types) == 0) {
    base::stop("Error: 'cell_types' must be a non-empty character vector.")
  }
  if (!base::is.list(sample_list) || base::length(sample_list) == 0) {
    base::stop("Error: 'sample_list' must be a non-empty list.")
  }
  if (!base::all(base::sapply(sample_list, base::is.matrix))) {
    base::stop("Error: All elements of 'sample_list' must be matrices.")
  }
  if (!base::is.data.frame(cell_IDs)) {
    base::stop("Error: 'cell_IDs' must be a dataframe.")
  }
  if (!all(c("cell_type", "cell_ID") %in% base::colnames(cell_IDs))) {
    base::stop("Error: 'cell_IDs' must contain 'cell_type' and 'cell_ID' columns.")
  }

  common_genes <- base::rownames(sample_list[[1]])
  if (base::is.null(common_genes)) {
    base::stop("Error: The sample matrix must have rownames representing genes.")
  }
  if (!base::all(base::sapply(sample_list, function(x) base::identical(base::rownames(x), common_genes)))) {
    base::stop("Error: All matrices in 'sample_list' must have identical rownames (genes).")
  }

  generate_individual <- function(celltype, sample_list, cell_IDs) {
    final_list <- list()

    valid_cell_IDs <- cell_IDs[cell_IDs$cell_type == celltype, "cell_ID"]

    for (i in base::seq_along(sample_list)) {
      sample <- sample_list[[i]]

      filtered_sample <- sample[, base::colnames(sample) %in% valid_cell_IDs, drop = FALSE]

      final_list[[paste0("sample_", i)]] <- filtered_sample
    }

    final_dataframe <- do.call(cbind, final_list)
    base::rownames(final_dataframe) <- common_genes
    final_matrix <- base::as.matrix(final_dataframe)
    final_dgCmatrix <- Matrix::Matrix(final_matrix, sparse = TRUE)

    return(final_dgCmatrix)
  }

  ct_list <- list()
  for (ct in cell_types){
    dgcmatrix <- generate_individual(ct, sample_list, cell_IDs)
    ct_list[[ct]] <- dgcmatrix
  }

  return(ct_list)
}

