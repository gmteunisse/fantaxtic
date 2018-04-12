#' Get the most abundant taxa from a phyloseq object
#'
#' This function subsets a phyloseq object to the top \eqn{n} most abundant
#' taxa. Users can choose whether to select by absolute counts or relative
#' counts, and whether to discard or collapse remaining taxa.
#'
#' @param physeq_obj A phyloseq object with an \code{otu_table} and a
#' \code{tax_table}.
#' @param n The number of top taxa to subset to.
#' @param relative Select taxa based on relative abundance?
#' (default = \code{TRUE})
#' @param discard_other Discard remainig taxa? (default = \code{FALSE})
#' @param other_label Label to give to collapsed taxa.
#' @return A phyloseq object
#' @examples
#' data(GlobalPatterns)
#' get_top_taxa(GlobalPatterns, 10)
#' get_top_taxa(GlobalPatterns, 10, relative = FALSE)
#' get_top_taxa(GlobalPatterns, 10, relative = TRUE, discard_other = TRUE)
#' get_top_taxa(GlobalPatterns, 10, relative = TRUE, other_label = "Non-abundant taxa")
#' @export
get_top_taxa <- function(physeq_obj, n, relative = TRUE, discard_other = FALSE, other_label = "Other"){

  #Define a temporary physeq object
  ps_tmp <- physeq_obj

  #Check for 0 entries
  smpl_sms <- phyloseq::sample_sums(ps_tmp)
  if (0 %in% smpl_sms){
    stop("Error: some samples contain 0 reads. These have to be removed to avoid
            downstream problems.")
  }

  #Extract the otu_table as a data.frame
  otu_tbl <- phyloseq::otu_table(ps_tmp)
  if (!phyloseq::taxa_are_rows(ps_tmp)){
    otu_tbl <- t(otu_tbl)
  }

  #Use relative abundances if requested
  if (relative){
    otu_tbl <- apply(otu_tbl, 2, function(x){
      x / sum (x)
    })
  }

  #Update the phyloseq object
  phyloseq::otu_table(ps_tmp) <- phyloseq::otu_table(otu_tbl, taxa_are_rows = T)

  #Get the top taxa names and discard or merge other taxa
  abun_taxa <- names(sort(phyloseq::taxa_sums(ps_tmp), decreasing = TRUE)[1:n])
  if (discard_other){
    physeq_obj <- phyloseq::prune_taxa(abun_taxa, physeq_obj)
  } else {
    to_merge <- phyloseq::taxa_names(physeq_obj)
    to_merge <- to_merge[!(to_merge %in% abun_taxa)]
    physeq_obj <- merge_taxa(physeq_obj, to_merge)
    tax_tbl <- phyloseq::tax_table(physeq_obj)
    indx <- which(row.names(tax_tbl) %in% to_merge)
    tax_tbl[indx,] <- other_label
    phyloseq::tax_table(physeq_obj) <- tax_tbl
  }
  return(physeq_obj)
}


