
#' Get the most abundant taxa over two nested levels
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
#' \code{merged_label} annotation. Next, loops through each remaining top level taxon,
#' and identifies the top  \eqn{m} named taxa at the nested level. If
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
#' @param by_proportion Converts absolute abundances to proportions before
#' calculating the summary statistic (default = \code{TRUE}).
#' @param merged_label Prefix that will be added after merging non-top taxa.
#' @param ... Additional arguments to be passed \code{top_taxa}
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
#' = "Family", n_top_taxa = 3, n_nested_taxa = 3, merged_label = "Other",
#' FUN = mean, na.rm = T)
#'
#' #' # Top 3 most abundant genera, top 3 most abundant species over all samples,
#' # using the mean as the aggregation function
#' nested_top_taxa(GlobalPatterns, top_tax_level = "Genus", nested_tax_level
#' = "Species", n_top_taxa = 3, n_nested_taxa = 3, merged_label = "Other",
#' FUN = mean, na.rm = T)
#' @export
nested_top_taxa <- function(ps_obj, top_tax_level, nested_tax_level,
                            n_top_taxa = 1, n_nested_taxa = 1,
                            top_merged_label = "Other",
                            nested_merged_label = "Other",
                            include_nested_tax = TRUE,
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

  # Make sure taxa are rows
  if (!phyloseq::taxa_are_rows(ps_obj)) {
    otu_table(ps_obj) <- phyloseq::otu_table(t(otu_table(ps_obj)), taxa_are_rows = T)
  }

  # Get the top n top_tax_level taxa
  top_top <- top_taxa(ps_obj, tax_level = top_tax_level, n_taxa = n_top_taxa,
                      by_proportion = by_proportion, ...)$top_taxa

  # Extract the taxonomy of to the top_tax_level and collapse all others
  ranks <- tax_ranks[1:which(tax_ranks == top_tax_level)]
  if (length(ranks) == 1){

    # Check only the top_level to
    # get the taxids for the nested levels and collapse all other taxa
    taxa <- tax_table(ps_obj) %>%
      data.frame(taxid = row.names(.)) %>%
      filter(!!as.symbol(top_tax_level) %in% top_top[[top_tax_level]])
    ps_obj_nest <- collapse_taxa(ps_obj, taxa$taxid, merged_label = top_merged_label)

  } else {

    # Also check a level above the top level to
    # get the taxids for the nested levels and collapse all other taxa
    higher_tax_level <- ranks[length(ranks) - 1]
    taxa <- apply(top_top[,ranks], 1, function(tax){
      tax_table(ps_obj) %>%
        data.frame(taxid = row.names(.)) %>%
        filter(!!as.symbol(top_tax_level) %in% tax[length(tax)],
               !!as.symbol(higher_tax_level) %in% tax[length(tax) - 1])
    })
    taxa <- do.call("rbind", taxa)
    ps_obj_nest <- collapse_taxa(ps_obj, taxa$taxid, merged_label = top_merged_label)
  }

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
  lvls <- unique(taxa[,top_tax_level])
  top_nest <- list()
  for (lvl in lvls){

    # Get the taxids of all nested_taxa within the current top_level taxon.
    # Do not include taxa where nested_level is NA
    taxids <- tax_table(ps_obj_nest) %>%
      data.frame(taxid = row.names(.)) %>%
      filter(!!as.symbol(top_tax_level) == lvl,
             !is.na(!!as.symbol(nested_tax_level))) %>%
      pull(taxid)

    # If no nested_level annotations are available, do include taxa
    # where nested_level is NA
    if (length(taxids) == 0){
      taxids <- tax_table(ps_obj_nest) %>%
        data.frame(taxid = row.names(.)) %>%
        filter(!!as.symbol(top_tax_level) == lvl) %>%
        pull(taxid)
    }

    # Remove all other taxa
    ps_obj_tmp <- collapse_taxa(ps_obj_nest, taxids, discard_other = T) %>%
      suppressWarnings()

    # Get the top n nested_tax_level taxa
    top_nest[[lvl]] <- top_taxa(ps_obj_tmp, n_taxa = n_nested_taxa, by_proportion = by_proportion, ...)$top_taxa %>%
      suppressWarnings()

    # Merge all other taxa, including the NA taxa
    to_merge <- taxa_names(ps_obj_tmp) %>%
      .[!. %in% top_nest[[lvl]]$taxid]
    na_taxids <- tax_table(ps_obj_nest) %>%
      data.frame(taxid = row.names(.)) %>%
      filter(!!as.symbol(top_tax_level) == lvl,
             is.na(!!as.symbol(nested_tax_level))) %>%
      pull(taxid)
    to_merge <- c(to_merge, na_taxids)
    ps_obj_nest <- merge_taxa(ps_obj_nest, to_merge, 1) %>%
      suppressWarnings()
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
  if (include_nested_tax){
    tax_tbl <- phyloseq::tax_table(ps_obj_nest) %>%
      data.frame() %>%
      mutate(!!nested_tax_level := ifelse(is.na(!!sym(nested_tax_level)),
                                          paste(nested_merged_label, !!sym(top_tax_level)),
                                          !!sym(nested_tax_level))
      ) %>%
      as.matrix()
  } else {
    tax_tbl <- phyloseq::tax_table(ps_obj_nest) %>%
      data.frame() %>%
      mutate(!!nested_tax_level := ifelse(is.na(!!sym(nested_tax_level)),
                                          nested_merged_label,
                                          !!sym(nested_tax_level))
      ) %>%
      as.matrix()
  }
  phyloseq::tax_table(ps_obj_nest) <- tax_tbl

  # Return a list of top and merged values
  return(list(ps_obj = ps_obj_nest,
              top_taxa = top))

}
