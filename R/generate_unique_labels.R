#' Generate unique labels.
#'
#' This function generates unique labels for an input vector by adding a
#' count to labels of which multiple occurrences are found.
#'
#' @param lbls A vector of labels for which to generate unique names.
#' @param sep_char A (combination of) characters with which to separate labels
#' and their unique counts.
#' @return A factor of unique labels.
#' @examples
#' lbls <- as.character(c("A", "B", "B", "C", "C", "C", "D"))
#' gen_uniq_lbls(lbls)
#' gen_uniq_lbls(lbls, sep_char = "_")
#'
#' #Also works with factors or numerics
#' gen_uniq_lvls(as.factor(lbls))
#' lbls <- as.numeric(c(1, 2, 2, 3, 3, 3, 4))
#' gen_uniq_lbls(lbls)
#' @export
gen_uniq_lbls <- function(lbls, sep_char = " "){
  lbls <- as.character(lbls)
  lbl_tbl <- as.data.frame(table(lbls))
  dupl <- which(lbl_tbl$Freq > 1)
  names(dupl) <- lbl_tbl$lbls[dupl]
  if(length(dupl) > 0){
    ind <- unlist(sapply(dupl, function(x){
      which(lbls %in% lbl_tbl$lbls[x])
    }))
    lbls_new <- unlist(sapply(dupl, function(x){
      paste(lbl_tbl$lbls[x], 1:lbl_tbl$Freq[x], sep = sep_char)
    }))
    lbls[ind] <- lbls_new
  }
  lbls <- factor(lbls, unique(lbls), ordered = T)
  return(lbls)
}
