% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_combined_matrix.R
\name{create_combined_matrix}
\alias{create_combined_matrix}
\title{Get a list of matrices with specified barcodes and genes}
\usage{
create_combined_matrix(specimenIDs, initial_path)
}
\arguments{
\item{specimenIDs}{A character vector of specimen IDs. Each ID corresponds to a subdirectory or file prefix in the data directory.}

\item{initial_path}{A string specifying the base path to the data directory containing subdirectories or files for each specimen. For example, if your files are in "data/data_raw/", that is your initial path.}
}
\value{
A list where each element is a sparse matrix. Rows represent genes, and columns represent barcodes.
}
\description{
Get a list of matrices with specified barcodes and genes
}
