#' Add names for missing taxon names in phyloseq objects
#'
#' This function adds names for OTUs/ASVs with incomplete taxonomic
#' annotations, i.e. annotation is only available up to a certain
#' taxonomic rank. It replaces \code{NA} values with the lowest
#' available taxonomic annotation for an OTU/ASV and a label indicating
#' that the annotation is unknown. Species can be renamed to genus +
#' species if desired.
#'
#' @param physeq_obj A phyloseq object with a \code{tax_table}.
#' @param unknown_label Label to prepend the taxon name with (default =
#' \code{"Unknown"}).
#' @return A phyloseq object
#' @examples
#' data(GlobalPatterns)
#' name_na_taxa(GlobalPatterns)
#' name_na_taxa(GlobalPatterns, unknown_label = "Unannotated")
#' @export
name_taxa <- function(physeq_obj, label = "Unknown", other_label = NULL, species = FALSE){

  #Get the tax table
  tax_tbl <- tax_table(physeq_obj)

  #Change any NA value to the lowest available taxonomic
  #annotation of that OTU/ASV
  tax_tbl <- t(apply(tax_tbl, 1, function(x){
    n <- length(x)
    if (sum(is.na(x)) == n){
      tax_ranks <- rep(label, n)
    } else {
      if (sum(is.na(x)) != 0){
        i <- min(which(is.na(x)))
        rank <- x[i-1] #The last known rank
        x[which(is.na(x))] <- sprintf("%s %s", label, rank)
      } else {
        if (!is.null(other_label)){
          if (other_label %in%  x){
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

  #Update the phyloseq object
  tax_table(physeq_obj) <- tax_tbl
  return(physeq_obj)
}
