#' Filter a matrix based on expression levels
#'
#' @param ctmat A numeric matrix where rows represent genes and columns represent samples (e.g., cells).
#' @param geneThr A numeric value between 0 and 1 specifying the minimum proportion of cells in which a gene must be expressed to be retained. The default is 0.02.
#' @param sampleThr A numeric value between 0 and 1 specifying the minimum rank percentile for the number of genes expressed in a sample to retain that sample. The default is 0.05.
#'
#' @return The filtered matrix.
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' ctmat <- matrix(sample(c(0, 1, 2, 3), 30, replace = TRUE),
#' nrow = 5, ncol = 6,
#' dimnames = list(
#' paste0("Gene", 1:5),
#' paste0("Cell", 1:6)))
#' cleanCtmat(ctmat)
#'
cleanCtmat <- function(ctmat, geneThr = 0.02, sampleThr = 0.05) {
  nCellsTot <- base::ncol(ctmat)
  nCellsThr <- geneThr * nCellsTot

  genesNN <- base::data.frame(gene = base::rownames(ctmat), nCells = base::rowSums(ctmat > 0), check.names = FALSE)  |>
    dplyr::filter(nCells > nCellsThr)
  ctmat <- ctmat[genesNN$gene, ]

  cellsNN <- base::data.frame(cell = base::colnames(ctmat), nGenes = base::colSums(ctmat > 0), check.names = FALSE) |>
    dplyr::mutate(perc = dplyr::percent_rank(nGenes)) |>
    dplyr::filter(perc > sampleThr)
  ctmat <- ctmat[, cellsNN$cell]

  return(ctmat)
}
