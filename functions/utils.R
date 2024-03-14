# convert back to what WGCNA needs
make_symmetric_matrix <- function(mat) {
  mat |> 
    Matrix::forceSymmetric(uplo = "L") |> 
    as.matrix()
}

make_lower_matrix <- function(mat) {
  mat |> 
    Matrix::tril()
}