#' Get the most abundant taxa from a phyloseq object
#'
#' This function identifies the top \eqn{n} taxa in a phyloseq object. Users specify the
#' summary statistic that is used to rank the taxa, e.g. \code{sum}, \code{mean} or
#' \code{median}. Furthermore, it is possible to add one or more grouping
#' factors from the \code{tax_table} to get group-specific top \eqn{n} taxa.
#'
#' @details
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
#' @param n_taxa The number of top taxa to identify.
#' @param grouping A character vector with the names of one or more grouping
#' factors found in the \code{sample_data}. To group by sample, specify \code{sample_id}.
#' @param by_proportion Converts absolute abundances to proportions before
#' calculating the summary statistic (default = \code{TRUE}).
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
#' top_taxa(GlobalPatterns, 10)
#'
#' # Top 10 most abundant ASVs over all samples by median abundance
#' top_taxa(GlobalPatterns, 10, FUN = median, na.rm = T)
#'
#' # Top 10 most abundant ASVs over all samples using absolute abundances
#' top_taxa(GlobalPatterns, 10, by_proportion = FALSE)
#'
#' # Top 2 most abundant ASVs per sample
#' top_taxa(GlobalPatterns, 2, grouping = "sample_id")
#'
#' # Top 2 most abundant ASVs per sample type
#' top_taxa(GlobalPatterns, 2, grouping = "SampleType")
#'
#' # Top 2 most abundant ASVs per sample type and group
#' set.seed(1)
#' sample_data(GlobalPatterns)$group <- as.factor(rbinom(nsamples(GlobalPatterns), 1, .5))
#' top_taxa(GlobalPatterns, 2, grouping = c("SampleType", "group"))
#'
#' # Top 2 most abundant genera per sample type
#' ps_genus <- tax_glom(GlobalPatterns, taxrank = "Genus")
#' top_taxa(ps_genus, 2)
#' @export
top_taxa <- function(ps_obj, n_taxa = 1, grouping = NULL, by_proportion = TRUE, FUN = mean, ...){

  # Make sure taxa are rows
  if (!phyloseq::taxa_are_rows(ps_obj)) {
    otu_table(ps_obj) <- phyloseq::otu_table(t(otu_table(ps_obj)), taxa_are_rows = T)
  }

  # Objects with a single taxon cause errors
  if (ntaxa(ps_obj) > 1){
    #Check for empty samples
    smpl_sms <- phyloseq::sample_sums(ps_obj)
    if (0 %in% smpl_sms){
      msg <- sprintf("Warning: some samples contain 0 reads. The following samples have been removed from the analysis:\n%s",
                     paste(names(smpl_sms)[smpl_sms == 0], collapse = "\n"))
      warning(msg)
      ps_obj <- phyloseq::subset_samples(ps_obj, smpl_sms > 0)
    }

    #Use relative abundances if requested
    if (by_proportion){
      otu_tab <- apply(otu_table(ps_obj), 2, function(x){
        x / sum (x)
      })
      phyloseq::otu_table(ps_obj) <- phyloseq::otu_table(otu_tab, taxa_are_rows = T)
    }
  } else {
    if (by_proportion){
      phyloseq::otu_table(ps_obj)[1,] <- 1
    }
  }


  # Get the top taxa
  if(!is.null(grouping)){

    # Check whether groups exist
    group_in_vars <- grouping %in% c(sample_variables(ps_obj), "sample_id")
    if(sum(group_in_vars) != length(grouping)){
      msg <- sprintf("Error: group(s) %s not found in sample_data of phyloseq object. Check your grouping argument.",
                     paste(grouping[!group_in_vars], collapse = ", "))
      stop(msg)
    }
  }

  # Get the top taxa per group
  top <- otu_table(ps_obj) %>%
    data.frame(taxid = row.names(.)) %>%
    pivot_longer(!taxid, names_to = "sample_id", values_to = "abundance") %>%
    left_join(sample_data(ps_obj) %>%
                data.frame(sample_id = row.names(.)),
              by = "sample_id") %>%
    group_by(across(all_of(!!grouping))) %>%
    group_by(taxid, .add = T) %>%
    summarise(abundance = FUN(abundance, ...), .groups = "keep") %>%
    group_by(across(all_of(!!grouping))) %>%
    slice_max(abundance, n = n_taxa) %>%
    arrange(desc(abundance), .by_group = TRUE) %>%
    mutate(tax_rank = rank(x = -abundance, ties.method = "first", na.last = T))

  # Add taxonomic annotations
  top <- ps_obj %>%
    tax_table() %>%
    data.frame(taxid = row.names(.)) %>%
    right_join(top, by = "taxid")

  # Reorder
  top <- top %>%
    relocate(all_of(!!grouping), tax_rank, taxid, abundance)

  return(top)
}
