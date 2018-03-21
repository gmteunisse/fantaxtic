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

  #Rotate table if necessary
  if (!taxa_are_rows(ps_tmp)){
    tbl <- t(as.matrix(otu_table(ps_tmp)))
    taxa_are_rows(tbl) <- TRUE
    otu_table(ps_tmp) <- otu_table(tbl, taxa_are_rows = TRUE)
  }

  #Check for 0 entries
  smpl_sms <- sample_sums(ps_tmp)
  if (0 %in% smpl_sms){
    warning("Warning: some samples contain 0 reads. These have been removed to avoid
            downstream problems.")
    ps_tmp <- prune_samples(smpl_sms != 0, ps_tmp)
  }

  #Use relative abundances if requested
  if (relative){
    otu_table(ps_tmp) <- otu_table(ps_tmp) / sample_sums(ps_tmp)
  }

  #Get the top taxa names and discard or merge other taxa
  abun_taxa <- names(sort(taxa_sums(ps_tmp), decreasing = TRUE)[1:n])
  if (discard_other){
    ps_tmp <- prune_taxa(abun_taxa, ps_tmp)
  } else {
    to_merge <- rownames(otu_table(ps_tmp))
    to_merge <- to_merge[!(to_merge %in% abun_taxa)]
    ps_tmp <- merge_taxa(ps_tmp, to_merge)
    tax_tbl <- tax_table(ps_tmp)
    indx <- which(row.names(tax_tbl) %in% to_merge)
    tax_tbl[indx,] <- other_label
    tax_table(ps_tmp) <- tax_tbl
  }
  return(ps_tmp)
}


