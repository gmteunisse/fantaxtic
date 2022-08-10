#' Add unique labels to identical taxa
#'
#' This function takes a phyloseq object and adds a label (\code{<id>}) to ASVs with the same
#' taxonomic annotation (\code{<tax>}). The label is either the ASV name (taken from row.names)
#' or the count of each duplicated taxon.
#'
#' @details
#' By default, this function finds the most specific non-NA annotation for each
#' ASV and adds labels to that level. However, it is also possible to add the
#' labels at a specific taxonomic level by specifying \code{tax_level = <rank>}.
#'
#' @param ps_obj Phyloseq object
#' @param tax_level The taxonomic level at which to identify duplicates. Can be
#' set to \code{NULL}. See details.
#' @param asv_as_id use ASVs as ids, rather than counts?
#' @param duplicate_label String that specifies how to construct the label. Must
#' contain at least \code{<tax>} and \code{<id>}
#' @return A phyloseq object
#' @import phyloseq dplyr tidyr stringr
#' @importFrom magrittr %>%
#' @examples
#' data(GlobalPatterns)
#'
#' # Label the lowest non-NA level
#' ps_tmp <- label_duplicate_taxa(GlobalPatterns)
#' View(tax_table(ps_tmp))
#'
#' # Use ASVs as ids rather than counts
#' ps_tmp <- label_duplicate_taxa(GlobalPatterns, asv_as_id = T, duplicate_label = "<tax> ASV <id>")
#' View(tax_table(ps_tmp))
#'
#' # Label at species leve lonly
#' ps_tmp <- label_duplicate_taxa(GlobalPatterns, tax_level = "Species")
#' View(tax_table(ps_tmp))
#' @export
label_duplicate_taxa <- function(ps_obj, tax_level = NULL, asv_as_id = F,
                                 duplicate_label = "<tax> <id>"){

  # Check arguments
  if(!grepl("<tax>", duplicate_label)){
    stop("Error: include '<tax>' in the duplicate_label")
  }
  # Check arguments
  if(!grepl("<id>", duplicate_label)){
    stop("Error: include '<id>' in the duplicate_label")
  }

  # Find duplicate taxa and flag for labeling
  tax_tab <- tax_table(ps_obj) %>%
    data.frame(row_name = row.names(.)) %>%
    group_by(across(-row_name)) %>%
    mutate(n = row_number(),
           dupl = n > 1) %>%
    mutate(add_label = sum(dupl) > 0)

  # Add label to the lowest non-NA taxonomic level
  if(is.null(tax_level)){

    # Convert to long
    tax_tab <- tax_tab %>%
      pivot_longer(!c(row_name, n, dupl, add_label),
                   names_to = "rank",
                   values_to = "tax") %>%
      filter(!is.na(tax)) %>%
      group_by(row_name)

    # Add label and fill in tax field
    tax_tab <- tax_tab %>%
      mutate(label = ifelse(add_label & row_number() == n(),
                            duplicate_label,
                            tax),
             tax = str_replace(label, "<tax>", tax))

    # Fill in id field
    if(asv_as_id){
      tax_tab <- tax_tab %>%
        mutate(tax = str_replace(tax, "<id>", as.character(row_name)))

    } else {
      tax_tab <- tax_tab %>%
        mutate(tax = str_replace(tax, "<id>", as.character(n)))
    }

    # convert to wide
    tax_tab <- tax_tab %>%
      select(!c(n, dupl, add_label, label)) %>%
      pivot_wider(names_from = rank, values_from = tax)

  # Add a label to a specific taxonomic rank
  } else {

    # Add label and fill in tax field
    tax_tab <- tax_tab %>%
      mutate(label = ifelse(add_label,
                            duplicate_label,
                            !!sym(tax_level)),
             !!sym(tax_level) := str_replace(label, "<tax>", !!sym(tax_level)))

    # Fill in id field
    if(asv_as_id){
      tax_tab <- tax_tab %>%
        mutate(!!sym(tax_level) := str_replace(!!sym(tax_level), "<id>", row_name))

    } else {
      tax_tab <- tax_tab %>%
        mutate(!!sym(tax_level) := str_replace(!!sym(tax_level), "<id>", as.character(n)))
    }
    tax_tab <- tax_tab %>%
      select(!c(n, dupl, add_label, label))
  }

  # Convert back to matrix and insert in ps_obj
  tax_mat <- tax_tab %>%
    as.matrix()
  row.names(tax_mat) <- tax_mat[,"row_name"]
  tax_mat <- tax_mat[,colnames(tax_mat) != "row_name"]
  tax_table(ps_obj) <- tax_mat

  return(ps_obj)




}
