#' Generate a colour for each taxon
#'
#' This function creates a named colour for each taxon in a phyloseq object
#'  at a specified taxonomic level.
#'
#' @details
#' If no palette is provided, colours are generated using the
#' \code{\link[ggnested]{nested_palette}} function from the
#' \code{ggnested} package. Only the taxonomic group with the \code{merged_label}
#' will get a custom colour. If a palette is provided, its colours will replace
#' the generated colours in the order in which they were specified. If a named
#' palette is provided, they will be used to replace specific colours.
#'
#' @param ps_obj A phyloseq object
#' @param merged_label The label of the merged taxa (usually "Other")
#' @param merged_clr The colour to assign to the merged taxa
#' @param palette A custom colour palette. Can be named.
#' @param base_clr the colour from which to generate the palette if no palette is
#' provided.
#' @return A named character vector
#' @import phyloseq dplyr tidyr ggnested
#' @importFrom magrittr %>%
#' @examples
#' data(GlobalPatterns)
#'
#' # Get top taxa and generate a palette
#' top <- top_taxa(GlobalPatterns, n_taxa = 10)
#' pal <- taxon_colours(top$ps_obj, tax_level = "Phylum")
#' scales::show_col(pal)
#'
#' # Generate a palette with a different base_clr
#' pal <- taxon_colours(top$ps_obj, tax_level = "Phylum", base_clr = "blue")
#' scales::show_col(pal)
#'
#' # Provide a custom incomplete palette
#' pal <- taxon_colours(top$ps_obj, tax_level = "Phylum", palette = c("blue", "pink"))
#' scales::show_col(pal)
#'
#' # Provide a custom incomplete named palette
#' pal <- taxon_colours(top$ps_obj, tax_level = "Phylum", palette = c(Cyanobacteria = "blue", Proteobacteria = "pink"))
#' scales::show_col(pal)
#'
#' @export
taxon_colours <- function(ps_obj, tax_level, merged_label = "Other", merged_clr = "grey90", palette = NULL, base_clr = "#008CF0"){

  # Generate a palette
  pal_df <- tax_table(ps_obj) %>%
    data.frame() %>%
    filter(!!sym(tax_level) != merged_label) %>%
    nested_palette(group = tax_level,
                   subgroup = tax_level,
                   base_clr = base_clr)

  # Extract the unique colours and add the merged_clr
  pal <- c(unique(pal_df$group_colour), merged_clr)

  # Update the names
  names(pal) <- c(unique(pal_df[[tax_level]]), merged_label)

  # Incorporate user colours
  if(!is.null(palette)){
    if(!is.null(names(palette))){
      pal[names(palette)] <- palette
    } else {
       n <- min(length(palette), length(pal) - 1)
      pal[1:n] <- palette[1:n]
    }
  }

  return(pal)
}
