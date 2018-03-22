## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = T, fig.align = "center", fig.width = 8-------------------
library(fantaxtic)
data(GlobalPatterns)

## ---- include = T, fig.align = "center", fig.width = 8-------------------

ps_tmp <- get_top_taxa(physeq_obj = GlobalPatterns, n = 10, relative = TRUE,
                       discard_other = FALSE, other_label = "Other")

## ---- include = T, fig.align = "center", fig.width = 8-------------------
ps_tmp <- name_taxa(ps_tmp, label = "Unkown", species = T, other_label = "Other")

## ---- include = T, fig.align = "center", fig.width = 8-------------------
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Species", other_label = "Other")

## ---- include = T, fig.align = "center", fig.width = 8-------------------
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Family", other_label = "Other",
               gen_uniq_lbls = TRUE)

## ---- include = T, fig.align = "center", fig.width = 8-------------------
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Family", other_label = "Other",
               order_alg = "other.abnd")

## ---- include = T, fig.align = "center", fig.width = 8-------------------
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Family", facet_by = "SampleType", other_label = "Other")

