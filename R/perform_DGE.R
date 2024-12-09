#' TMM Normalize and perform DEA
#'
#' @param PB The output of create_pseudobulk.
#' @param design To specify the groups of interest. E.g. design <- model.matrix(~ group, data = PB$meta)
#'
#' @return The results of DEA with the edgeR package.
#' @export
perform_DGE <- function(PB, design) {

  if (!base::all(base::colnames(PB$expr) %in% PB$meta$group_sample)) {
    stop("Sample names in expression matrix do not match group sample names in meta data")
  }

  y <- edgeR::DGEList(counts = PB$expr, group = PB$meta$group)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  x <- edgeR::estimateDisp(y, design, trend.method = "locfit", tagwise = TRUE, prior.df = NULL)
  fit <- edgeR::glmFit(x, design = design)
  test <- edgeR::glmLRT(fit)

  res <- edgeR::topTags(test, n = Inf) |>
    as.data.frame() |>
    tibble::rownames_to_column('gene') |>
    dplyr::mutate(test = 'pseudobulk_edgeR')

  return(res)
}
