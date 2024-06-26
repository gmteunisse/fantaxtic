#' Plot a nested bar plot
#'
#' This function takes a phyloseq object and plots the abundance of taxa in each
#' sample at two levels: a top level (e.g. Phylum), using colours, and a nested
#' level (e.g. Species), using shades and tints of each colour. It is intended
#' to be used together with \code{\link{top_taxa}} or \code{\link{nested_top_taxa}},
#' but functions with any phyloseq object. The output is a \code{ggplot} object,
#' and can be manipulated as such, although take note of the documentation for
#' \code{\link[ggnested]{ggnested}}.
#'
#' @details
#' This function is a wrapper around \code{\link[ggnested]{ggnested}}, with some
#' accessory functions to ensure the plot comes out neatly It runs through
#' the following steps:
#' \itemize{
#'    \item{Generate a palette}{ using \code{\link{taxon_colours}}}
#'    \item{Name NA taxa}{ using \code{\link{name_na_taxa}}}
#'    \item{Label identical taxa}{ using \code{\link{label_duplicate_taxa}}}
#'    \item{Convert phyloseq to a data frame}{ using \code{\link[phyloseq]{psmelt}}}
#'    \item{Refactor merged taxa}{ using \code{\link{move_label}} and \code{\link{move_nested_labels}}}
#'    \item{Generate a nested barplot}{ using \code{\link[ggnested]{ggnested}}}
#' }
#'
#' See the documentation for each function for details.
#'
#' @param ps_obj A phyloseq object, ideally generated by \code{\link{top_taxa}} or
#' \code{\link{nested_top_taxa}}
#' @param top_level The taxonomic level with which to fill the bars
#' @param nested_level The taxonomic level at which to generate shades and tints
#' @param top_merged_label The \code{(top_)merged_label} that was used in
#' \code{\link{top_taxa}} or \code{\link{nested_top_taxa}}.
#' @param nested_merged_label The \code{nested_merged_label} that was used in
#' \code{\link{nested_top_taxa}}.
#' @param palette A custom palette. See \code{\link{taxon_colours}}.
#' @param base_clr A base colour for palette generation. See \code{\link{taxon_colours}}.
#' @param merged_clr A colour for the merged taxon. See \code{\link{taxon_colours}}.
#' @param include_rank Include rank when naming NA taxa? See \code{\link{name_na_taxa}}.
#' @param na_taxon_label Label format for NA taxa. See \code{\link{name_na_taxa}}.
#' @param asv_as_id Use ASVs as IDs for duplicate taxa? See \code{\link{label_duplicate_taxa}}.
#' @param duplicate_taxon_label Label format for duplicate taxa. See \code{\link{label_duplicate_taxa}}.
#' @param relative_abundances Plot relative abundances?
#' @param sample_order Character vector with all sample names in order.
#' @param ... options to be passed to \code{\link[ggnested]{ggnested}}
#' @return A ggplot2 object
#' @import phyloseq dplyr tidyr ggnested
#' @importFrom magrittr %>%
#' @examples
#'
#' # Plot the Phylum and Species of the 20 most abundant ASVs
#' ps_obj <- GlobalPatterns
#' top <- top_taxa(ps_obj, n_taxa = 20,
#'                 FUN = median)
#' p <- plot_nested_bar(top$ps_obj,
#'                      "Phylum",
#'                      "Species")
#' p
#'
#' # Order the samples by the abundance of Other and plot with faceting
#' sample_order <- p$data %>%
#'     group_by(Sample) %>%
#'     mutate(Abundance = Abundance / sum(Abundance)) %>%
#'     filter(Phylum == "Other") %>%
#'     arrange(Abundance) %>%
#'     pull(Sample)
#' plot_nested_bar(top$ps_obj,
#'                 "Phylum",
#'                 "Species",
#'                 sample_order = sample_order) +
#'   facet_wrap(~SampleType, scales = "free_x")
#'
#' # Plot the nested_top_taxa
#' top_level <- "Phylum"
#' nested_level <- "Species"
#' top_merged_label <- "Other"
#' nested_merged_label <- "Other"
#' top <- nested_top_taxa(ps_obj,
#'                        top_tax_level = top_level,
#'                        nested_tax_level = nested_level,
#'                        n_top_taxa = 3,
#'                        n_nested_taxa = 3)
#' plot_nested_bar(top$ps_obj, top_level, nested_level)
#'
#' # Order by SampleType and alter the colours for some taxa
#' sample_order <- sample_data(top$ps_obj) %>%
#' data.frame() %>%
#'   arrange(SampleType, X.SampleID) %>%
#'   pull(X.SampleID) %>%
#'   as.character()
#' plot_nested_bar(top$ps_obj,
#'                 top_level = top_level,
#'                 nested_level = "Genus",
#'                 palette = c(Bacteroidetes = "red",
#'                            Proteobacteria = "blue"),
#'                 sample_order = sample_order)
#'
#' @export
plot_nested_bar <- function(ps_obj,
                            top_level,
                            nested_level,
                            top_merged_label = "Other",
                            nested_merged_label = "Other <tax>",
                            palette = NULL,
                            base_clr = "#008CF0",
                            merged_clr = "grey90",
                            include_rank = T,
                            na_taxon_label = "<tax> (<rank>)",
                            asv_as_id = F,
                            duplicate_taxon_label = "<tax> <id>",
                            relative_abundances = T,
                            sample_order = NULL,
                            ...){

  # Create labels
  ps_tmp <- ps_obj %>%
    name_na_taxa(include_rank = include_rank,
                 na_label = na_taxon_label)
  ps_tmp <- ps_tmp %>%
    label_duplicate_taxa(tax_level = nested_level,
                         asv_as_id = asv_as_id,
                         duplicate_label = duplicate_taxon_label )

  # Generate a palette
  pal <- taxon_colours(ps_tmp,
                       tax_level = top_level,
                       merged_label = top_merged_label,
                       merged_clr = merged_clr,
                       palette = palette,
                       base_clr = base_clr)

  # Convert physeq to df
  psdf <- psmelt(ps_tmp)

  # Move the merged labels to the appropriate positions
  psdf <- move_label(psdf = psdf,
                     col_name = top_level,
                     label = top_merged_label,
                     pos = 0)
  psdf <- move_nested_labels(psdf,
                             top_level = top_level,
                             nested_level = nested_level,
                             top_merged_label = top_merged_label,
                             nested_label = gsub("<tax>", "", nested_merged_label),
                             pos = Inf)

  # Reorder samples
  if(!is.null(sample_order)){
    if(all(sample_order %in% unique(psdf$Sample))){
      psdf <- psdf %>%
        mutate(Sample = factor(Sample, levels = sample_order))
    } else {
      stop("Error: not all(sample_order %in% sample_names(ps_obj)).")
    }

  }

  # Generate a bar plot
  p <- ggnested(psdf,
                aes_string(main_group = top_level,
                           sub_group = nested_level,
                           x = "Sample",
                           y = "Abundance"),
                ...,
                main_palette = pal) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_nested(theme_light) +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
  if (relative_abundances){
    p <- p + geom_col(position = position_fill())
  } else {
    p <- p + geom_col()
  }
  return(p)

}
