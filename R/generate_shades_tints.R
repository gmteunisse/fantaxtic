#' Generate shade and tint variants for a base color.
#'
#' \code{gen_shades_tints} returns \code{i/2} shades and \code{i/2} tints
#' of a base color.
#'
#' @param i The number of tints and/or shades to be generated.
#' @param clr The base color to use (hex or named color).
#' @param incl.base Should the base color be part of the generated palette?
#' @return \code{i} variants of a base color.
#' @rdname gen_shades_tints
#' @examples
#' gen_shades(3)
#' gen_tints(3)
#' gen_shades_tints(3)
#' gen_shades(3, "blue")
#' gen_tints(3, "blue")
#' gen_shades_tints(3, "blue")
#' @export
gen_shades_tints <- function(i, clr = "#6495ed"){
  if (i == 0){
    return(c())
  }
  if (i == 1){
    return(clr)
  }
  if (i == 2){
    n_tnts <- 2
    n_shds <- 1
    tnts <- gen_tints(n_tnts, clr, incl.base = F)[1]
  } else {
    n_tnts <- ceiling(i/2)
    n_shds <- i - n_tnts
    tnts <- gen_tints(n_tnts, clr, incl.base = F)
  }
  shds <- gen_shades(n_shds, clr, incl.base = T)
  clrs <- c(shds, tnts)
  return(clrs)
}

#' \code{gen_shades} returns \code{i} shades of a base color.
#' @rdname gen_shades_tints
gen_shades <- function(i, clr = "#6495ed", incl.base = FALSE){
  if (i == 0){
    return(c())
  }
  if (i == 1){
    return(clr)
  }
  clr_rgb <- t(col2rgb(clr))
  shds <- (1 + round(i)):(2 * round(i) + 1)
  shds <- sapply(shds, function(x){
    shd_rgb <- clr_rgb * x / (2 * round(i) + 1)
    shd <- rgb(shd_rgb[1], shd_rgb[2], shd_rgb[3], maxColorValue = 255)
    return(shd)
  })
  if (incl.base){
    return(shds[1:i+1])
  }
  return(shds[1:i])
}

#' \code{gen_tints} returns \code{i} tints of a base color.
#' @rdname gen_shades_tints
gen_tints <- function(i, clr = "#6495ed", incl.base = FALSE){
  if (i == 0){
    return(c())
  }
  if (i == 1){
    return(clr)
  }
  clr_rgb <- t(col2rgb(clr))
  tnts <- 0:i
  tnts <- sapply(tnts, function(x){
    tnt_rgb <- clr_rgb + (255 - clr_rgb) * x / (round(i) * 1.5)
    tnt <- rgb(tnt_rgb[1], tnt_rgb[2], tnt_rgb[3], maxColorValue = 255)
    return(tnt)
  })
  if (incl.base){
    return(tnts[1:i])
  }
  return(tnts[1:i+1])
}
