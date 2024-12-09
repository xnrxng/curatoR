#' Get a list of matrices with specified barcodes and genes
#'
#' @param specimenIDs A character vector of specimen IDs. Each ID corresponds to a subdirectory or file prefix in the data directory.
#' @param initial_path A string specifying the base path to the data directory containing subdirectories or files for each specimen. For example, if your files are in "data/data_raw/", that is your initial path.
#'
#' @return A list where each element is a sparse matrix. Rows represent genes, and columns represent barcodes.
#' @export
#'

create_combined_matrix <- function(specimenIDs, initial_path) {
  final_list <- base::list()

  for (patient in specimenIDs){
    barcode_file <- base::paste0(initial_path, patient, "barcodes.tsv")
    barcodes <- data.table::fread(barcode_file, header = FALSE) |>
      dplyr::pull()

    gene_file <- base::paste0(initial_path, patient, "features.tsv")
    genes <- data.table::fread(gene_file, header = FALSE) |>
      dplyr::select(V2) |>
      dplyr::pull()

    matrix_file <- base::paste0(initial_path, patient, "matrix.mtx")
    matrix_expr <- Matrix::readMM(matrix_file)
    base::rownames(matrix_expr) <- genes
    base::colnames(matrix_expr) <- barcodes
    final_list[[patient]] <- matrix_expr
  }
  return(final_list)
  }
