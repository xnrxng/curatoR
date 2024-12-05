generate_cell_types <- function(cell_type, sample_list, prefix) {
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
  final_dgCmatrix <- Matrix::as(final_matrix, "dgCMatrix")

  return(final_dgCmatrix)
}
