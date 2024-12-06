library(tibble)

# creating a simple dataset
A <- c(0, 1, 2, 3, 4)
B <- c(5, 6, 7, 8, 9)
C <- c(10, 11, 12, 13, 14)
D <- c(15, 16, 17, 18, 19)
simple_dataset_w_col <- tibble::tibble(A, B, C, D)
simple_dataset_w_col

# creating a dataset with only one column
one_dataset_w_col <- tibble::tibble(A = A)
one_dataset_w_col
