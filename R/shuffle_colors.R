#' Reshuffle a list of colors.
#'
#' This function reshuffles a vector of colours so that in a list of \code{n}
#' colors, every color \code{i} in \code{\{1, n / 2\}} will be placed next
#' to color \code{i + n / 2}. This is useful for visualization of
#' automatically generated color palettes based upon for example gradients
#' or HSL color circles.
#'
#' @param clrs A vector of colors to reshuffle.
#' @return A reshuffled vector of \code{clrs}.
#' @examples
#' shuffle_colors(gen_colors(6))
#' @export
shuffle_colors <- function(clrs){
  n.col <- length(clrs)
  ordr <- rep(1:ceiling(n.col/2), each = 2, length.out = n.col)
  indx <- seq(2, length(ordr), 2)
  ordr[indx] <- ordr[indx] + ceiling(n.col/2)
  clrs <- clrs[ordr]
  return(clrs)
}
