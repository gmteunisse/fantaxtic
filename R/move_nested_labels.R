#' Relevel a nested factor in a psmelt dataframe
#'
#' Moves a nested factor level in a fantaxtic psmelt dataframe to the desired position.
#'
#' @param psdf A psmelt dataframe
#' @param top_level The name of the top level factor column
#' @param nested_level The name of the nested level factor column
#' @param top_merged_label The label for the top level merged taxon
#' @param nested_label The nested label to be moved
#' @param pos The position to move the level to.
#' @return A data.frame
#' @import dplyr forcats
#' @importFrom magrittr %>%
#' @examples
#' data(GlobalPatterns)
#'
#' # Move "Other" to position Inf
#' top <- nested_top_taxa(GlobalPatterns, top_tax_level = "Phylum", nested_tax_level = "Species", n_top_taxa = 3, n_nested_taxa = 3)
#' psdf <- psmelt(top$ps_obj)
#' psdf <- move_label(psdf, col_name = "Phylum", label =  "Other", pos = 0)
#' psdf <- move_nested_labels(psdf, top_level = "Phylum", nested_level = "Species", pos = Inf)
#' levels(psdf$Species)
#' @export
move_nested_labels <- function(psdf,
                               top_level,
                               nested_level,
                               top_merged_label = "Other",
                               nested_label = "Other",
                               pos = Inf){

  # Split into dataframe per top level
  df_list <- psdf %>%
    mutate(!!sym(top_level) := as.character(!!sym(top_level)),
           !!sym(nested_level) := as.character(!!sym(nested_level))) %>%
    filter(!!sym(top_level) != top_merged_label) %>%
    group_by(!!sym(top_level)) %>%
    group_split()

  # Reorder the labels so that the label comes at the end for each group
  nested_levels <- lapply(df_list, function(df){

    # Get the entire label
    lab <- df %>%
      filter(grepl(nested_label, !!sym(nested_level))) %>%
      pull(nested_level)

    # If the merged group was found, move the label to the end
    if(length(lab) > 0){
      df %>%
        pull(nested_level) %>%
        fct_relevel(lab[1], after = pos) %>%
        levels()

      # If not, return the other levels in alphabetical order
    } else {
      df %>%
        pull(nested_level) %>%
        unique() %>%
        sort()
    }
  }) %>% unlist()

  # Add the label of the merged top level
  nested_levels <- c(nested_levels, top_merged_label)

  # Refactor and sort
  psdf %>%
    mutate(!!sym(nested_level) := factor(!!sym(nested_level), levels = nested_levels, ordered = T)) %>%
    arrange(top_level, nested_level)
}
