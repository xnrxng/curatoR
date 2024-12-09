#' CPM normalizes a matrix and log-transforms it
#'
#' @param mtx A matrix.
#' @param log TRUE or FALSE. Whether to log-transform the matrix as well. Default is false.
#'
#' @return A CPM normalized matrix.
#' @export
#'
#' @examples
#' #' set.seed(123)
#' ctmat <- matrix(sample(c(0, 1, 2, 3), 30, replace = TRUE),
#' nrow = 5, ncol = 6,
#' dimnames = list(
#' paste0("Gene", 1:5),
#' paste0("Cell", 1:6)))
#' do_cpm_log(ctmat, log = TRUE)
do_cpm_log <- function(mtx, log = FALSE) {
  colsums <- base::colSums(mtx)
  cpm_result <- base::t(base::t(mtx) / colsums * 1e6)

  if (log) {
    cpm_result <- base::log1p(cpm_result)
  }

  return(cpm_result)
}
