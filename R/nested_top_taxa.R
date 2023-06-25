
#' Get the most abundant taxa over two taxonomic levels
#'
#' This function identifies the top \eqn{n} named taxa and the top \eqn{m} named taxa at a
#' nested level in a phyloseq object. Users specify the
#' summary statistic that is used to rank the taxa, e.g. \code{sum}, \code{mean} or
#' \code{median}. Furthermore, it is possible to add one or more grouping
#' factors from the \code{tax_table} to get group-specific top \eqn{n,m} taxa.
#'
#' @details
#' This function first finds the top \eqn{n} named taxa at the top level, after
#' which it merges all other top_taxa level into a single taxon with the
#' \code{merged_label} annotation. Next, it loops through each remaining top level taxon,
#' and identifies the top \eqn{m} named taxa at the nested level. If
#' \eqn{\le m} taxa are available, it will only return those taxa. If more are
#' available, it will merge all non-top-taxa into a single taxon with the
#' \code{merged_label} annotation, together with its top level annotation.
#' If no named taxa are available at the nested_level, all taxa will be merged
#' into a single taxon with \code{merged_label} annotation, together with its
#' \code{top_level} annotation. Thus, the \code{merged_label} taxon overall and
#' in each group represents the combination of taxa without an annotation and
#' taxa with an annotation that were not in the top \eqn{(n,m)} abundant taxa.
#'
#' If \code{nested_tax_level = "ASV"}, \code{row.names(tax_table(ps_obj))} will
#' be added as an ASV column to the tax_table, unless this column already exists.
#'
#' The top taxa can be identified based on the absolute abundances or proportions.
#' When using absolute abundances, please make sure to normalize or rarefy the data
#' before using this function. If \code{by_proportion = TRUE}, abundances will
#' be converted to relative abundance before applying \code{FUN}.
#'
#' @param ps_obj A phyloseq object with an \code{otu_table} and a
#' \code{tax_table}.
#' @param top_tax_level The name of the top taxonomic rank in the phyloseq object
#' @param nested_tax_level The name of the nested taxonomic rank in the phyloseq object
#' @param n_top_taxa The number of top taxa to identify at the top level.
#' @param n_nested_taxa The number of top taxa to identify at the nested level.
#' For ASVs, specify "ASV"
#' @param top_merged_label Label to assign to the merged top_tax_level taxa
#' @param nested_merged_label Label to assign to the merged nested_tax_level taxa
#' @param by_proportion Converts absolute abundances to proportions before
#' calculating the summary statistic (default = \code{TRUE}).
#' @param ... Additional arguments to be passed \code{\link{top_taxa}}
#' (e.g. \code{grouping = <string>, FUN = mean, na.rm = TRUE},).
#' @return A list in which \code{top_taxa} is a tibble with the rank, taxon id, grouping
#' factors, abundance summary statistic and taxonomy of the top taxa and \code{ps_obj}
#' is the phyloseq object after collapsing all non-top taxa.
#' @import phyloseq dplyr tidyr
#' @importFrom magrittr %>%
#' @examples
#' data(GlobalPatterns)
#'
#' # Top 3 most abundant orders, top 3 most abundant families over all samples,
#' # using the mean as the aggregation function
#' nested_top_taxa(GlobalPatterns, top_tax_level = "Order", nested_tax_level
#' = "Family", n_top_taxa = 3, n_nested_taxa = 3,
#' FUN = mean, na.rm = T)
#'
#' #' # Top 1 most abundant genera, top 2 most abundant species per SampleType,
#' # using the median as the aggregation function
#' nested_top_taxa(GlobalPatterns, top_tax_level = "Genus", nested_tax_level
#' = "Species", n_top_taxa = 1, n_nested_taxa = 2, grouping = "SampleType",
#' FUN = median)
#' @export
nested_top_taxa <- function(ps_obj,
                            top_tax_level,
                            nested_tax_level,
                            n_top_taxa = 1,
                            n_nested_taxa = 1,
                            top_merged_label = "Other",
                            nested_merged_label = "Other <tax>",
                            by_proportion = T,
                            ...){

  # Check arguments
  tax_ranks <- c(rank_names(ps_obj))
  if (nested_tax_level == "ASV" & !"ASV" %in% tax_ranks){
    tax_ranks <- c(tax_ranks, "ASV")
  }
  if (is.null(top_tax_level)){
    stop("Error: no top_tax_level provided")
  } else if (!top_tax_level %in% tax_ranks){
    msg <- sprintf("Error: top_tax_level %s not in rank_names(ps_obj)", top_tax_level)
    stop(msg)
  }
  if (is.null(nested_tax_level)){
    stop("Error: no nested_tax_level provided")
  } else if (!nested_tax_level %in% tax_ranks){
    msg <- sprintf("Error: nested_tax_level %s not in rank_names(ps_obj)", nested_tax_level)
    stop(msg)
  }
  if(top_tax_level == nested_tax_level){
    msg <- sprintf("Error: top_tax_level and nested_tax_level cannot both be %s", top_tax_level)
    stop(msg)
  }
  if(which(tax_ranks == top_tax_level) > which(tax_ranks == nested_tax_level)){
    msg <- sprintf("Error: top_tax_level needs to be a higher taxonomic rank than nested_tax_level.
                    Tax ranks in this object are %s", paste(tax_ranks, collapse = ", "))
    stop(msg)
  }
  if(!grepl("<tax>", nested_merged_label)){
    stop("Error: nested merged label must include <tax>")
  }

  # Make sure taxa are rows
  if (!phyloseq::taxa_are_rows(ps_obj)) {
    otu_table(ps_obj) <- phyloseq::otu_table(t(otu_table(ps_obj)), taxa_are_rows = T)
  }

  # Get the top n top_tax_level taxa
  top_top <- top_taxa(ps_obj, tax_level = top_tax_level, n_taxa = n_top_taxa,
                      by_proportion = by_proportion, ...)$top_taxa
  ranks <- tax_ranks[1:which(tax_ranks == top_tax_level)]
  top_top_taxa <- apply(top_top[,ranks,drop=FALSE], 1, paste, collapse = '_')

  # Extract the taxonomy of the top_tax_level and collapse all others
  taxa <- tax_table(ps_obj) %>%
    data.frame() %>%
    select(ranks) %>%
    unite(col = 'full_taxonomy') %>%
    filter(full_taxonomy %in% top_top_taxa)
  taxids <- row.names(taxa)
  ps_obj_nest <- collapse_taxa(ps_obj, taxids, merged_label = top_merged_label)

  # Add an ASV column if working at ASV level
  if (nested_tax_level == "ASV"){
    ranks <- rank_names(ps_obj_nest)
    if(!"ASV" %in% ranks){
      tax_table(ps_obj_nest) <- cbind(tax_table(ps_obj_nest),
                                      ASV = row.names(tax_table(ps_obj_nest)))
    }
  } else {

    # If not ASVs, glom at the nested level.
    ps_obj_nest <- tax_glom(ps_obj_nest, taxrank = nested_tax_level, NArm = F)
  }

  # Loop through each top_tax_level and get the top n nested_tax_level taxa
  lvls <- unique(taxa$full_taxonomy)
  top_nest <- list()
  merged_ids <- c()
  for (lvl in lvls){

    # Get the taxids of all nested_taxa within the current top_level taxon.
    all_taxids <- taxa %>%
      filter(full_taxonomy == lvl) %>%
      row.names(.)
    all_taxids <- all_taxids[all_taxids %in% taxa_names(ps_obj_nest)]

    # Select named taxa only
    named_taxids <- tax_table(ps_obj_nest) %>%
     data.frame(taxid = row.names(.)) %>%
      filter(taxid %in% all_taxids,
             !is.na(!!as.symbol(nested_tax_level))
             ) %>%
      pull(taxid)

    # Get the top taxa
    if (length(named_taxids) != 0){

      # Remove all other taxa
      ps_obj_tmp <- collapse_taxa(ps_obj_nest, named_taxids, discard_other = T) %>%
        suppressWarnings()

      # Get the top n nested_tax_level taxa
      top_nest[[lvl]] <- top_taxa(ps_obj_tmp, n_taxa = n_nested_taxa, by_proportion = by_proportion, ...)$top_taxa %>%
        suppressWarnings()

      # Find all taxa to merge
      to_merge <- all_taxids[!all_taxids %in% top_nest[[lvl]]$taxid]
    } else {
      to_merge <- all_taxids
    }


    # Merge the non-top taxa within this top-level taxon
    ps_obj_nest <- merge_taxa(ps_obj_nest, to_merge, 1) %>%
      suppressWarnings()

    # Store the merged ids
    if(length(to_merge) > 0 ){
      merged_ids <- c(merged_ids, to_merge[1])
    }
  }

  # Combine the results
  top_nest <- do.call("rbind", top_nest)

  # Add top abundances
  top <- top_top %>%
    select(where(~!all(is.na(.x)))) %>%
    select(!c(taxid)) %>%
    rename(top_abundance = abundance,
           top_tax_rank = tax_rank) %>%
    left_join(top_nest, .) %>%
    rename(nested_abundance = abundance,
           nested_tax_rank = tax_rank) %>%
    relocate(taxid, top_abundance, nested_abundance, top_tax_rank, nested_tax_rank) %>%
    suppressMessages()

  # Update the taxon name to nested_merged_label
  tax_tbl <- phyloseq::tax_table(ps_obj_nest) %>%
    as.data.frame() %>%
    as.matrix()
  for (i in merged_ids){
    lab <- gsub("<tax>", tax_tbl[i,top_tax_level], nested_merged_label)
    tax_tbl[i, is.na(tax_tbl[i,])] <- lab
  }
  phyloseq::tax_table(ps_obj_nest) <- tax_table(tax_tbl)

  # Return a list of top and merged values
  return(list(ps_obj = ps_obj_nest,
              top_taxa = top))

}




