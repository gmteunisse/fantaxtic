#' Add names for missing taxon names in phyloseq objects
#'
#' This function adds names for OTUs/ASVs with incomplete taxonomic
#' annotations, i.e. annotation is only available up to a certain
#' taxonomic rank. It replaces \code{NA} values with the lowest
#' available taxonomic annotation for an OTU/ASV and a label indicating
#' that the annotation is unknown. Species can be renamed to genus +
#' species if desired. To avoid downstream problem, it numbers
#' taxa with the same taxonomic annotation but different sequences.
#'
#' @param physeq_obj A phyloseq object with a \code{tax_table}.
#' @param label Label to prepend the taxon name with (default =
#' \code{"Unknown"}).
#' @param other_label The label(s) of samples whose names should not be altered.
#' @param species Generate a 'Genus species' (i.e. 'Escherichia Coli') label
#' for the species level?
#' @param unique_rank The taxonomic rank by which to generate unique labels
#' (added number) if desired (default = \code{"Unannotated"}).
#' @param unique_sep The text character(s) by which to separate the annotation
#' and the unique number (only when \code{!is.null(unique_rank)}).
#' @examples
#' data(GlobalPatterns)
#' name_na_taxa(GlobalPatterns)
#' name_na_taxa(GlobalPatterns, unknown_label = "Unannotated")
#' @export
name_taxa <- function(physeq_obj, label = "Unannotated", other_label = NULL, species = FALSE, unique_rank = NULL, unique_sep = " "){

  #Get the tax table
  tax_tbl <- phyloseq::tax_table(physeq_obj)

  #Store the names
  tax_names <- colnames(tax_tbl)

  #Change any NA value to the lowest available taxonomic
  #annotation of that OTU/ASV
  tax_tbl <- t(apply(tax_tbl, 1, function(x){
    n <- length(x)
    if (sum(is.na(x)) == n){
      tax_ranks <- rep("Unknown", n)
    } else {
      if (sum(is.na(x)) != 0){
        i <- max(which(!is.na(x)))
        rank <- x[i] #The last known rank
        x[which(is.na(x))] <- sprintf("%s %s (%s)", label, rank, names(x)[i])
      } else {
        if (!is.null(other_label)){
          if (sum(other_label %in%  x) > 0){
            tax_ranks <- x
          } else {
            if (species){
              x[n] <- sprintf("%s %s", x[n-1], x[n])
            }
          }
        } else {
          if (species){
            x[n] <- sprintf("%s %s", x[n-1], x[n])
          }
        }
      }
      tax_ranks <- x
      return(tax_ranks)
    }
  }))

  #Generate unique labels, i.e. add a number on the desired taxonomic level
  if (!is.null(unique_rank)){
    ind <- which(colnames(tax_tbl) == unique_rank)
    tax_tbl[,ind] <- as.character(tax_tbl[,ind])
    tax_tbl[,ind] <- as.character(gen_uniq_lbls(tax_tbl[,ind], sep_char = unique_sep))
  }

  #Update the phyloseq object
  phyloseq::tax_table(physeq_obj) <- tax_tbl
  colnames(phyloseq::tax_table(physeq_obj)) <- tax_names
  return(physeq_obj)
}
