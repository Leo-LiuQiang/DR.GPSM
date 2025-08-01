#' Create all-pairwise contrast matrix
#'
#' @param levels Character vector of treatment levels
#' @param ref    Reference treatment level (default = first level)
#'
#' @return A contrast matrix of size Choose(K, 2) × K with entries -1, 0, or 1
#' @export
build_contrast <- function(levels, ref = NULL)
{
  if (is.null(ref)) ref <- levels[1L]
  stopifnot(ref %in% levels)

  K <- length(levels)
  mat <- matrix(0, nrow = choose(K, 2), ncol = K,
                dimnames = list(NULL, levels))
  row <- 1L

  for (i in seq_len(K - 1L))
    for (j in (i + 1L):K) {
      mat[row, i] <- -1
      mat[row, j] <-  1
      row <- row + 1L
    }

  if (!is.null(ref)) {
    refCol <- match(ref, levels)
    for (r in seq_len(nrow(mat)))
      if (mat[r, refCol] == 0)
        next
    else if (mat[r, refCol] ==  1)
      mat[r, ] <- -mat[r, ]
  }

  rn <- character(nrow(mat))
  for (r in seq_len(nrow(mat))) {
    neg <- levels[ mat[r, ] == -1 ]
    pos <- levels[ mat[r, ] ==  1 ]
    rn[r] <- paste0(pos, "v", neg)
  }
  rownames(mat) <- rn
  colnames(mat) <- levels

  return(mat)
}
