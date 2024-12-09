#' CPMLog normalize and perform DEA
#'
#' @param design To specify the groups of interest. E.g. design <- model.matrix(~ group, data = PB$meta)
#' @param PB he output of create_pseudobulk.
#'
#' @return The results of DEA with the limma package.
#' @export
limma_dge <- function(design, PB){
  x = limma::voom(base::as.matrix(PB$expr), design)
  fit = limma::lmFit(x, design) |> limma::eBayes()
  res = fit |> limma::topTable(number = Inf)
  return(res)
}
