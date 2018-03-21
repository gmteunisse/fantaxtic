#' Generate a color palette.
#'
#' This function generates \code{n} maximally separated colors with equal
#' saturation and brightness from a base color.
#'
#' @param n The number of colors to generate.
#' @param clr The base color from which to generate other colors.
#' @return A vector of \code{n} colors.
#' @examples
#' gen_colors(5)
#' gen_colors(5, "blue")
#' @export
gen_colors <- function(n, clr = "#6495ed"){
  if (n == 1){
    return(clr)
  }
  clr.rgb <- t(col2rgb(clr))
  clr.hsv <- t(rgb2hsv(clr.rgb[1], clr.rgb[2], clr.rgb[3], 255))
  clrs <- 0:(n-1)
  clrs <- sapply(clrs, function(x){
    offset <- x/(n)
    h <- (clr.hsv[1] + offset) %% 1
    s <- clr.hsv[2]
    v <- clr.hsv[3]
    clr <- hsv(h, s, v)
    return(clr)
  })
  clrs <- shuffle_colors(clrs)
  return(clrs)
}
