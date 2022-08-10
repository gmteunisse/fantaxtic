#' Relevel a factor in a psmelt dataframe
#'
#' Moves a factor level in a psmelt dataframe to the desired position.
#'
#' @param psdf A psmelt dataframe
#' @param col_name The name of the factor column
#' @param label The factor level to be moved
#' @param pos The position to move the level to.
#' @return A data.frame
#' @import dplyr forcats
#' @importFrom magrittr %>%
#' @examples
#' data(GlobalPatterns)
#'
#' # Move "Other" to position 0
#' top <- top_taxa(GlobalPatterns)
#' psdf <- psmelt(top$ps_obj)
#' psdf <- move_label(psdf, col_name = "Phylum", label =  "Other", pos = 0)
#' levels(psdf$Phylum)
#' @export
move_label <- function(psdf, col_name, label = "Other", pos = 0){
  psdf %>%
    mutate(!!sym(col_name) := factor(!!sym(col_name), ordered = T)) %>%
    mutate(!!sym(col_name) := fct_relevel(.f = !!sym(col_name), label, after = pos)) %>%
    arrange(col_name)
}
