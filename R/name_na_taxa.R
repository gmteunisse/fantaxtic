#' Name NA taxa
#'
#' This function fills in any NA annotation with the lowest known taxonomic annotation (\code{<tax>})
#' from a higher rank (\code{<rank>}).
#'
#' @param ps_obj A phyloseq object
#' @param include_rank Whether to include the rank of the lowest known annotation
#' @param na_label The label to assign to each filled in NA annotation.
#' @return A phyloseq object
#' @import phyloseq dplyr tidyr
#' @importFrom magrittr %>%
#' @examples
#' data(GlobalPatterns)
#'
#' # Fill in names for NA taxa, including their rank
#' ps_tmp <- name_na_taxa(GlobalPatterns)
#' View(tax_table(ps_tmp))
#'
#' # Leave the rank out and alter the label
#' ps_tmp <- name_na_taxa(GlobalPatterns, include_rank = F, na_label = "Unknown <tax>")
#' View(tax_table(ps_tmp))
#' @export
name_na_taxa <- function(ps_obj, include_rank = T, na_label = "Unknown <tax> (<rank>)"){

  # Check arguments
  if(!grepl("<tax>", na_label)){
    stop("Error: include '<tax>' in the na_label")
  }
  if (include_rank){
    if(!grepl("<rank>", na_label)){
      stop("Error: include_rank = TRUE; include '<rank>' in the na_label")
    }
  } else {
    if(grepl("<rank>", na_label)){
      stop("Error: include_rank = FALSE; remove '<rank>' from the na_label")
    }
  }

  # Convert to long data
  taxa_long <- tax_table(ps_obj) %>%
    data.frame(row_name = row.names(.)) %>%
    pivot_longer(!row_name,
                 names_to = "rank",
                 values_to = "tax")

  # Fill in NAs using the value above
  taxa_long <- taxa_long %>%
    mutate(na = is.na(tax)) %>%
    group_by(row_name) %>%
    fill(tax)

  # Create na_labels
  taxa_long <- taxa_long %>%
    mutate(expr = ifelse(na,
                         na_label,
                         tax),
           na_label = str_replace(expr, "<tax>", tax))

  # Add the last annotated rank
  if (include_rank){
    taxa_long <- taxa_long %>%
      mutate(last_rank = ifelse(na,
                                NA,
                                rank)) %>%
      fill(last_rank) %>%
      mutate(na_label = str_replace(na_label, "<rank>", last_rank))
  }

  # Convert back to tax_table
  taxa_mat <- taxa_long %>%
    select(row_name, rank, na_label) %>%
    pivot_wider(names_from = rank, values_from = na_label) %>%
    as.matrix()
  row.names(taxa_mat) <- taxa_mat[,"row_name"]
  taxa_mat <- taxa_mat[,colnames(taxa_mat) != "row_name"]
  tax_table(ps_obj) <- taxa_mat

  return(ps_obj)
}
