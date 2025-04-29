#' Get MEX files
#'
#' @param path Parent folder containing each matrix.
#' @param sample_list List of sample names.
#' @param meta Metadata file containing the metadata information of each sample.
#' @param final_path Final parent directory where each MEX file will be found.
#'
#' @return A folder where each sample has its MEX files. Also a metadata file of each cell.
#' @export
#'
get_mex_files <- function(path, sample_list, meta, final_path) {

  final_metadata <- list()

  for (sample in sample_list) {
    dir.create(file.path(final_path), recursive = TRUE)
    dir.create(file.path(paste0(final_path, "/", sample)))

    sample_file <- paste0(path, "/", sample, "-annotated_matrix.txt")
    sample_data <- fread(sample_file, header = TRUE, sep = "\t", data.table = FALSE)
    rownames(sample_data) <- sample_data[, 1]
    sample_data <- sample_data[, -1]
    colnames(sample_data) <- gsub(" ", "_", colnames(sample_data))
    colnames(sample_data) <- gsub("\\..*", "", colnames(sample_data))
    cell_types <- colnames(sample_data)

    colnames(sample_data) <- make.unique(colnames(sample_data))
    new_columns <- paste0(colnames(sample_data), "_", sample)
    colnames(sample_data) <- new_columns

    write.table(
      new_columns,
      file = gzfile(paste0(final_path, "/", sample, "/barcodes.tsv.gz")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      col.names = FALSE
    )

    genes <- rownames(sample_data)
    write.table(
      genes,
      file = gzfile(paste0(final_path, "/", sample, "/features.tsv.gz")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      col.names = FALSE
    )

    mat <- as.matrix(sample_data)
    sparse_mat <- as(mat, "dgCMatrix")
    writeMM(sparse_mat, file = paste0(final_path, "/", sample, "/matrix.mtx"))
    gzip(paste0(final_path, "/", sample, "/matrix.mtx"), destname = paste0(final_path, "/", sample, "/matrix.mtx.gz"))


    individual_metadata <- meta[meta$Individual_ID == sample, ]

    sample_meta <- data.frame(
      cell_id = new_columns,
      cell_type = cell_types)

    for (col in names(individual_metadata)) {
      sample_meta[[col]] <- rep(individual_metadata[[col]], nrow(sample_meta))
    }

    final_metadata[[sample]] <- sample_meta

      gc()
    }

  all_metadata <- do.call(rbind, final_metadata)

  write_tsv(all_metadata, paste0(final_path, "/cell_metadata.tsv"))
}
