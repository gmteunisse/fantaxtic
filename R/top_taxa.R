#' Get the most abundant taxa from a phyloseq object
#'
#' This function identifies the top \eqn{n} taxa in a phyloseq object. Users specify the
#' summary statistic that is used to rank the taxa, e.g. \code{sum}, \code{mean} or
#' \code{median}. Furthermore, it is possible to add one or more grouping
#' factors from the \code{tax_table} to get group-specific top \eqn{n} taxa.
#'
#' @details
#' When \code{tax_level = NULL}, the analysis will be done at the ASV level. If a
#' \code{tax_level} is specified, the object will first be glommed using
#' \code{tax_glom(ps_obj, tax_rank = tax_level, NArm = F)} at the specified level.
#' This can lead to taxa with NA annotations at the specified \code{tax_level}. By default,
#' these taxa will not be considered for the analysis, but they can be included
#' by setting \code{include_na_taxa = T}.
#'
#' This function, together with \code{\link[fantaxtic]{collapse_taxa}}, replaces
#' \code{\link[fantaxtic]{get_top_taxa}}. Identical output can be obtained  by setting
#' \code{FUN = sum}.
#'
#' The top taxa can be identified based on the absolute abundances or proportions.
#' When using absolute abundances, please make sure to normalize or rarefy the data
#' before using this function. If \code{by_proportion = TRUE}, abundances will
#' be converted to relative abundance before applying \code{FUN}.
#'
#' @param ps_obj A phyloseq object with an \code{otu_table} and a
#' \code{tax_table}.
#' @param tax_level Optional taxonomic level at which to get the top taxa.
#' @param n_taxa The number of top taxa to identify.
#' @param grouping A character vector with the names of one or more grouping
#' factors found in the \code{sample_data}. To group by sample, specify \code{sample_id}.
#' @param by_proportion Converts absolute abundances to proportions before
#' calculating the summary statistic (default = \code{TRUE}).
#' @param include_na_taxa When \code{tax_level} is specified, include NA taxa? See details.
#' @param merged_label The label to assign to merged taxa
#' @param FUN Function that returns a single summary statistic from an input vector,
#' e.g. \code{sum}, \code{mean} (default) or \code{median}
#' @param ... Additional arguments to be passed to \code{FUN}.
#' @return A tibble with the rank, taxon id, grouping factors, abundance summary
#' statistic and taxonomy.
#' @import phyloseq dplyr tidyr
#' @importFrom magrittr %>%
#' @examples
#' data(GlobalPatterns)
#'
#' # Top 10 most abundant ASVs over all samples
#' top_taxa(GlobalPatterns, n_taxa = 10)
#'
#' # Top 10 most abundant ASVs over all samples by median abundance
#' top_taxa(GlobalPatterns, n_taxa = 10, FUN = median, na.rm = T)
#'
#' # Top 10 most abundant ASVs over all samples using absolute abundances
#' top_taxa(GlobalPatterns, n_taxa = 10, by_proportion = FALSE)
#'
#' # Top 2 most abundant ASVs per sample
#' top_taxa(GlobalPatterns, n_taxa = 2, grouping = "sample_id")
#'
#' # Top 2 most abundant ASVs per sample type
#' top_taxa(GlobalPatterns, n_taxa = 2, grouping = "SampleType")
#'
#' # Top 2 most abundant ASVs per sample type and group
#' set.seed(1)
#' sample_data(GlobalPatterns)$group <- as.factor(rbinom(nsamples(GlobalPatterns), 1, .5))
#' top_taxa(GlobalPatterns, n_taxa = 2, grouping = c("SampleType", "group"))
#'
#' # Top 20 most abundant genera
#' top_taxa(GlobalPatterns, n_taxa = 20, tax_level = "Genus")
#'
#' #' # Top 20 most abundant genera including NAs
#' top_taxa(GlobalPatterns, n_taxa = 20, tax_level = "Genus", include_na_taxa = T)
#'
#' @export
top_taxa <- function(ps_obj, tax_level = NULL, n_taxa = 1, grouping = NULL,
                     by_proportion = TRUE, include_na_taxa = F, merged_label = "Other",
                     FUN = mean, ...){

  # Make sure taxa are rows
  if (!phyloseq::taxa_are_rows(ps_obj)) {
    otu_table(ps_obj) <- phyloseq::otu_table(t(otu_table(ps_obj)), taxa_are_rows = T)
  }

  # Check arguments
  tax_ranks <- rank_names(ps_obj)
  if(!is.null(tax_level)){
    if (tax_level == "ASV" & !"ASV" %in% tax_ranks){
      tax_ranks <- c(tax_ranks, "ASV")
    }
    if (!tax_level %in% tax_ranks){
      msg <- sprintf("Error: tax_level %s not in rank_names(ps_obj)", tax_rank)
      stop(msg)
    }
  }

  # Check the grouping
  if(!is.null(grouping)){

    # Check whether groups exist
    group_in_vars <- grouping %in% c(sample_variables(ps_obj), "sample_id")
    if(sum(group_in_vars) != length(grouping)){
      msg <- sprintf("Error: group(s) %s not found in sample_data of phyloseq object. Check your grouping argument.",
                     paste(grouping[!group_in_vars], collapse = ", "))
      stop(msg)
    }
  }

  # Glom and merge taxa when tax_level is specified
  merged_taxid <- NULL
  if(!is.null(tax_level)){
    if(tax_level != "ASV"){

      # Glom taxa
      ps_tmp <- tax_glom(ps_obj, taxrank = tax_level, NArm = FALSE)

      # Merge NA taxa
      if (!include_na_taxa){
        taxids <- tax_table(ps_tmp) %>%
          data.frame(taxid = row.names(.)) %>%
          filter(!is.na(!!sym(tax_level))) %>%
          pull(taxid)
        ps_tmp <- collapse_taxa(ps_tmp, taxa_to_keep = taxids, discard_other = F, merged_label = merged_label)

        #Get the taxid
        merged_taxid <- tax_table(ps_tmp) %>%
          data.frame(taxid = row.names(.)) %>%
          filter(!!sym(tax_level) == merged_label) %>%
          pull(taxid)
      }
    } else {
      ps_tmp <- ps_obj
    }
  } else {
    ps_tmp <- ps_obj
  }

  # Convert to long
  ps_long <- otu_table(ps_tmp) %>%
    data.frame(taxid = row.names(.)) %>%
    rename_with(function(x){gsub("X", "", x)}) %>%
    pivot_longer(!taxid, names_to = "sample_id", values_to = "abundance") %>%
    left_join(sample_data(ps_tmp) %>%
                data.frame(sample_id = row.names(.)),
              by = "sample_id")

  # Calculate relative abundances
  if (by_proportion){
    ps_long <- ps_long %>%
      group_by(sample_id) %>%
      mutate(abundance = abundance / sum(abundance))

    # Set NaNs to 0
    ps_long[is.nan(ps_long$abundance),'abundance'] <- 0
  }

  # Get the top taxa
  top <- ps_long %>%
    group_by(across(all_of(!!grouping))) %>%
    group_by(taxid, .add = T) %>%
    summarise(abundance = FUN(abundance, ...), .groups = "keep") %>%
    group_by(across(all_of(!!grouping))) %>%
    filter(!taxid %in% merged_taxid) %>%
    slice_max(abundance, n = n_taxa) %>%
    arrange(desc(abundance), .by_group = TRUE) %>%
    mutate(tax_rank = rank(x = -abundance, ties.method = "first", na.last = T))

  # Add taxonomic annotations
  top <- ps_tmp %>%
    tax_table() %>%
    data.frame(taxid = row.names(.)) %>%
    right_join(top, by = "taxid")

  # Reorder
  top <- top %>%
    relocate(all_of(!!grouping), tax_rank, taxid, abundance)

  # Update the phyloseq object
  ps_tmp <- collapse_taxa(ps_tmp, taxa_to_keep = top$taxid, discard_other = F,
                          merged_label = merged_label)

  return(list(ps_obj = ps_tmp,
              top_taxa = top))

}


