#' Generate a color palette.
#'
#' This function generates \eqn{i_1, i_2, ... i_n} shades and tints specified
#' in \code{clr.tbl} for \eqn{n} colors \code{clr.pal}. If no \code{clr.tbl}
#' is specified, this function generates \eqn{n} maximally separated colors
#' with equal saturation and brightness from a base color.
#'
#' @param clr_tbl A \code{data.frame} or \code{matrix} with two columns: the
#' first column gives the labels to which colors need to be assigned; the
#' second column gives the number of tints and shades that need to be generated
#' for that color.
#' @param clr_pal A color palette with with \eqn{n} colors, where
#' \eqn{n =} \code{nrow(clr_tbl)}.
#' @param base_clr The base color from which to generate a color palette.
#' @return A list of \eqn{n} items, each containing \eqn{i_n} shades.
#' @examples
#' clr_tbl <- data.frame(Var = c("Phylum A", "Phylum B"),
#'                       Freq = c(5, 8))
#' gen_palette(clr_tbl)
#' gen_palette(clr_tbl, clr_pal = c("red", "blue"))
#' gen_palette(clr_tbl, base_clr = "blue")
#' @export
gen_palette <- function(clr_tbl, clr_pal = NULL, base_clr = "#6495ed"){

  #Define a palette
  base_pal <- c("#6495ed", "#ff7256", "#edbc64", "#8470ff", "#8ee5ee", "#EE8DD6")

  #Define the number of variants required per color
  clr_vars <- clr_tbl[,2]
  n_clrs <- length(clr_vars)

  #Check arguments and generate colors if required
  if (is.null(clr_pal)){
    if (base_clr != "#6495ed"){
      clr_pal <- gen_colors(n_clrs, base_clr)
    } else {
      if (n_clrs <= length(base_pal)){
        clr_pal <- base_pal[1:n_clrs]
      } else {
        clr_pal <- gen_colors(n_clrs, base_clr)
      }
    }
  } else {
    if (length(clr_pal) < n_clrs){
      stop(sprintf("Error: %d values required in clr.pal, %d provided.", n_clrs, length(clr_pal)))
    }
    clr_pal <- clr_pal[1:n_clrs]
  }

  #Give warnings when many colors are present
  if (n_clrs > 10){
    warning("Warning: using > 10 base colors, results may be hard to read")
  }
  if (max(clr_vars) > 7){
    warning("Warning: generating > 7 shades and tints; results may be hard to read.")
  }

  #Print the translation table from color.by to their central colors
  clr_mtrx <- cbind(clr_tbl, clr_pal)
  names(clr_mtrx) <- c("Level", "N.color.shades", "Central.color")
  print(clr_mtrx)

  #Generate i vars for each clr in  clr.pal
  clr_mtrx <- cbind(clr_vars, clr_pal)
  clr_list <- split(clr_mtrx, seq(nrow(clr_mtrx)))
  palette <- lapply(clr_list, function(v){
    i <- as.numeric(v[1])
    clr <- v[2]
    gen_shades_tints(i, clr)
  })
  return(palette)
}
