#' Subset a phyloseq object to a set of taxa.
#'
#' This function takes a phyloseq object and a list of taxon ids to be kept,
#' and discards or merges all other taxa.
#'
#' @details
#' This function is essentially a wrapper around \link[phyloseq]{prune_taxa} and
#' \link[phyloseq]{merge_taxa}.
#'
#' This function, together with \link[fantaxtic]{top_taxa} replaces
#' \link[fantaxtic]{get_top_taxa}. Identical output can be obtained by setting
#' \code{FUN = sum} in \link[fantaxtic]{top_taxa}.
#'
#' @param ps_obj Phyloseq object
#' @param taxa_to_keep taxon ids (taxids) to be kept. These taxon ids need to
#' be part of \code{taxa_names(ps_obj)}
#' @param discard_other if \code{TRUE}, any taxon not in \code{taxa_to_keep} is
#' discard. If \code{FALSE}, they are collapsed into a single taxon labeled
#' \code{other_label}
#' @param merged_label Label for the new merged taxon
#' @return A phyloseq object
#' @import phyloseq
#' @examples
#' data(GlobalPatterns)
#'
#' # Top 10 most abundant ASVs over all samples, collapse other ASVs into 'Low abundance'
#' top <- top_taxa(GlobalPatterns, 10)
#' ps_collapsed <- collapse_taxa(GlobalPatterns,
#'                               taxa_to_keep = top$top_taxa$taxid,
#'                               merged_label = "Low abundance")
#'
#' # Top 10 most abundant ASVs over all samples, discard other taxa
#' top <- top_taxa(GlobalPatterns, 10)
#' ps_collapsed <- collapse_taxa(GlobalPatterns, taxa_to_keep = top$taxid,
#'                               discard_other = TRUE)
#'
#' # Keep genus Clostridium, collapse all others into Other genera
#' ps_tmp <- subset_taxa(GlobalPatterns, Genus == "Clostridium")
#' taxids <- taxa_names(ps_tmp)
#' ps_collapsed <- collapse_taxa(GlobalPatterns, taxa_to_keep = taxids,
#'                               merged_label = "Other genera")
#'
#' @export
collapse_taxa <- function(ps_obj, taxa_to_keep, discard_other = FALSE, merged_label = "Other"){

  # Make sure taxa are rows
  if (!phyloseq::taxa_are_rows(ps_obj)) {
    phyloseq::otu_table(ps_obj) <- phyloseq::otu_table(t(phyloseq::otu_table), taxa_are_rows = T)
  }

  # Merge or discard all other taxa
  if (discard_other) {
    ps_obj <- phyloseq::prune_taxa(taxa_to_keep, ps_obj)
  } else {

    # Merge taxa
    to_merge <- phyloseq::taxa_names(ps_obj)
    to_merge <- to_merge[!(to_merge %in% taxa_to_keep)]
    ps_obj <- merge_taxa(ps_obj, to_merge, archetype = 1)

    # Update the taxon name to merged_label
    tax_tbl <- phyloseq::tax_table(ps_obj)
    indx <- which(row.names(tax_tbl) %in% to_merge)
    for (i in indx){
      tax_tbl[i, is.na(tax_tbl[i,])] <- merged_label
    }
    phyloseq::tax_table(ps_obj) <- tax_tbl
  }
  return(ps_obj)
}
