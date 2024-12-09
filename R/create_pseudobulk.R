#' Create pseudobulks
#'
#' @param expr Expression matrix where columns are cell IDs and rows are gene names.
#' @param meta Dataframe where each row should represent a cell. Cell IDs should be the rownames. It needs to have the columns "patientID", "disorder", "sex", and "age". "disorder" should have the values "yes" or "no".
#'
#' @return A list containing the pseudobulk matrix and its metadata.
#' @export
#'
#' @examples
#'set.seed(123)
#'expr <- matrix(sample(0:10, 50, replace = TRUE),
#'nrow = 5, ncol = 10,
#' dimnames = list(paste0("Gene", 1:5), paste0("Cell", 1:10)))
#' meta <- data.frame(
#' patientID = c("P1", "P2", "P3", "P1",
#' "P1", "P2", "P2", "P3", "P3", "P1"),
#' disorder = c("yes", "no", "no", "yes",
#' "yes", "no", "no", "no", "no", "yes"),
#' sex = c("male", "female", "female", "male",
#' "male", "female", "female", "female", "female", "male"),
#' age = c(70, 68, 64, 70, 70, 68, 68, 64, 64, 70),
#' row.names = paste0("Cell", 1:10)
#' )
#'
#' create_pseudo_bulk(expr, meta)
#'
create_pseudo_bulk <- function(expr, meta) {

  mm <- stats::model.matrix(~ 0 + patientID:disorder, data = meta)
  mat_mm <- expr %*% mm

  mat_mm = base::as.matrix(mat_mm)
  keep_genes <- base::rowSums(mat_mm > 0) > 0
  mat_mm <- mat_mm[keep_genes, ] |> base::as.matrix() |> base::as.data.frame()
  mat_mm <- mat_mm |> base::as.matrix() |> base::as.data.frame()

  base::colnames(mat_mm) <- base::gsub("replicate|label", "", base::colnames(mat_mm))

  keep_samples <- base::colSums(mat_mm) > 0
  mat_mm <- mat_mm[, keep_samples]
  mat_mm <- mat_mm[, ]

  targets <- base::data.frame(group_sample = base::colnames(mat_mm))|>
    dplyr::mutate(group = base::gsub(".*\\:", "", group_sample),
           patientID = base::sub(".*patientID(.*)\\:.*", "\\1", group_sample))

  targets <- targets |>
    dplyr::left_join(meta |> dplyr::select(patientID, disorder, sex, age) |> dplyr::distinct(),
              by = "patientID")

  targets$group <- base::as.factor(targets$group)
  targets$group <- stats::relevel(targets$group, ref = "disorderno")

  PB <- base::list(meta = targets, expr = mat_mm)

  return(PB)
}
