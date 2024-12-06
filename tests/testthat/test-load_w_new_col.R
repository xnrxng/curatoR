test_that("function should return a dataset as a tibble with new column names", {
  dataset_path <- "https://raw.githubusercontent.com/DSCI-310-2024/dsci-310_group-7_wine-quality-prediction/function-load-data/tests/testthat/simple_dataset.csv"
  col_names <- c("A", "B", "C", "D")
  delimiter <- ","
  loaded_data <- load_w_new_col(dataset_path, col_names, delimiter)
  expect_s3_class(loaded_data, "data.frame")
  expect_equal(ncol(loaded_data), 4)
  expect_equal(colnames(loaded_data), col_names)
  expect_equal(loaded_data, simple_dataset_w_col)
})

test_that("function should return a dataset as a tibble with new column names
          even for datasets with only one column", {
            dataset_path <- "https://raw.githubusercontent.com/DSCI-310-2024/dsci-310_group-7_wine-quality-prediction/function-load-data/tests/testthat/one_dataset.csv"
            col_names <- c("A")
            delimiter <- ","
            loaded_data <- load_w_new_col(dataset_path, col_names, delimiter)
            expect_s3_class(loaded_data, "data.frame")
            expect_equal(ncol(loaded_data), length(col_names))
            expect_equal(colnames(loaded_data), col_names)
            expect_equal(loaded_data, one_dataset_w_col)
          })

test_that("function should not read empty datasets", {
  dataset_path <- "https://raw.githubusercontent.com/DSCI-310-2024/dsci-310_group-7_wine-quality-prediction/function-load-data/tests/testthat/empty_dataset.csv"
  col_names <- c("A", "B")
  delimiter <- ","
  expect_error(loaded_data <- load_w_new_col(dataset_path, col_names, delimiter))
})

test_that("function should raise an error when delimiter is unquoted", {
  dataset_path <- "https://raw.githubusercontent.com/DSCI-310-2024/dsci-310_group-7_wine-quality-prediction/function-load-data/tests/testthat/simple_dataset.csv"
  col_names <- c("A", "B", "C", "D")
  expect_error(load_w_new_col(dataset_path, col_names, ,))
})

test_that("function should raise an error when a vector with strings is not given", {
  dataset_path <- "https://raw.githubusercontent.com/DSCI-310-2024/dsci-310_group-7_wine-quality-prediction/function-load-data/tests/testthat/simple_dataset.csv"
  col_names <- c(A, B, C, D)
  delimiter <- ","
  expect_error(load_w_new_col(dataset_path, col_names, delimiter))
})

test_that("function should raise an error when no delimiter is given", {
  dataset_path <- "https://raw.githubusercontent.com/DSCI-310-2024/dsci-310_group-7_wine-quality-prediction/function-load-data/tests/testthat/simple_dataset.csv"
  col_names <- c("A", "B", "C", "D")
  expect_error(load_w_new_col(dataset_path, col_names))
})

test_that("function should raise an error when no dataset path is given", {
  col_names <- c("A", "B", "C", "D")
  delimiter <- ","
  expect_error(load_w_new_col(NULL, col_names, delimiter))
})

test_that("function should raise an error when no vector is given", {
  dataset_path <- "https://raw.githubusercontent.com/DSCI-310-2024/dsci-310_group-7_wine-quality-prediction/function-load-data/tests/testthat/simple_dataset.csv"
  delimiter <- ","
  expect_error(load_w_new_col(dataset_path, NULL, delimiter))
})

test_that("function should raise an error when number of columns does not match
          number of items in the vector", {
            dataset_path <- "https://raw.githubusercontent.com/DSCI-310-2024/dsci-310_group-7_wine-quality-prediction/function-load-data/tests/testthat/simple_dataset.csv"
            col_names <- c("A", "B", "C")
            delimiter <- ","
            expect_error(load_w_new_col(dataset_path, col_names, delimiter))
          })
