#' Load dataset with new column names
#'
#' Given an URL of an online dataset, the datasets' delimiters,
#' and a vector of new column names, loads the dataset as a tibble
#' with the new column names.
#'
#' @param dataset_path An URL or path that leads to a dataset.
#' @param col_names Vector with the new column names as strings (quoted).
#' @param delimiter The dataset's delimiter (could be ".", ";", ",", etc), quoted.
#'
#' @return Loaded dataset with the new column names, as a tibble.
#'
#' @export
#'
#' @examples
#' load_w_new_col("https://raw.githubusercontent.com/plotly/datasets/master/mtcars.csv",
#' c("model", "milespergallon", "cylinder_number",
#' "displacement", "horsepower", "rear_axle_ratio",
#' "weight", "quartermiletime", "engine", "transmission",
#' "forwardgears", "carb"), ",")

load_w_new_col <- function(dataset_path, col_names, delimiter) {
  # returns loaded dataset with new column names
  if (!is.character(col_names) || any(!is.character(col_names))) {
    stop("`col_names` must be a vector of strings (quoted)")
  }

  loaded_dataset <- readr::read_delim(toString(dataset_path), delim = delimiter)

  if (ncol(loaded_dataset) != length(col_names)) {
    stop("Number of columns in the dataset does not match the length of vector of column names")
  }

  colnames(loaded_dataset) <- col_names
  return(loaded_dataset)
}
