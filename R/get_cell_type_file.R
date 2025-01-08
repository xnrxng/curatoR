#' Get cleaned cell type annotations
#'
#' @param file_path A string that denotes where the file is.
#' @param remove_prefix TRUE or FALSE. Whether to remove the barcodes' prefix.
#' @param remove_suffix TRUE or FALSE. Whether to remove the barcodes' suffix.
#' @param id_col The name of the barcodes column.
#' @param type_col The name of the cell type column.
#'
#' @return A cleaned data table with two columns, cell_id and cell_type.
#' @export
get_cell_type_file <- function(file_path, remove_prefix = FALSE, remove_suffix = FALSE, id_col = NULL, type_col = NULL){
  # read the file
  raw_file <- data.table::fread(file_path)

  if (!is.null(id_col)){
    if (id_col %in% names(raw_file)) {
      data.table::setnames(raw_file, id_col, "cell_id")
    } else {
      stop(paste("Specified id_col", id_col, "does not exist in the file."))
    }
  }
  else{
  # find column that looks like barcodes and rename it to cell_id
  barcode_col <- grep("barcode|cell_id|cellid|barcodes|orig.ident", names(raw_file), ignore.case = TRUE, value = TRUE)
  if (length(barcode_col) == 1) {
    data.table::setnames(raw_file, barcode_col, "cell_id")
  } else {
    stop("Could not identify a unique column for barcodes.")
  }}

  if (!is.null(type_col)){
    if (type_col %in% names(raw_file)) {
      data.table::setnames(raw_file, type_col, "cell_type")
    } else {
      stop(paste("Specified type_col", type_col, "does not exist in the file."))
    }
    }else{
  # find column that looks like cell types and rename it to cell_type
  cell_type_col <- grep("celltype|cell_type|type|cluster|clusters|seurat_clusters", names(raw_file), ignore.case = TRUE, value = TRUE)
  if (length(cell_type_col) == 1) {
    data.table::setnames(raw_file, cell_type_col, "cell_type")
  } else {
    stop("Could not identify a unique column for cell types.")
  }}

  raw_file <- raw_file |>
    dplyr::select(cell_id, cell_type)

  if (remove_prefix){
    raw_file[, cell_id := sub("^[^\\t\\s\\-_]*[\\t\\s\\-_]+", "", cell_id)]
  }

  if (remove_suffix){
    raw_file[, cell_id := sub("[-_\\t\\s][^-_\\t\\s]*$", "", cell_id)]
  }

  raw_file <- unique(raw_file)

  more_celltype <- raw_file |>
    dplyr::group_by(cell_id) |>
    dplyr::filter(dplyr::n_distinct(cell_type) > 1) |>
    dplyr::ungroup()

  raw_file <- raw_file[!(raw_file$cell_id %in% more_celltype$cell_id), ]

  if(nrow(more_celltype) >= 1){
    print("These cell IDs had more than one cell type assigned, and hence, removed.")
    print(more_celltype$cell_id)
  }

  readr::write_tsv(raw_file, "cleaned_annotations.tsv")

  return(raw_file)
}
