# Fantaxtic
This is Fantaxtic, a package that contains a set of wrapper functions for the phyloseq and ggplot2 packages in the R software that turns ordinary taxonomic count data into fantaxtic, publication ready plots. Fantaxtic figures come as completely stylized figures with pleasing colors and maximal plotting control. As the output figures are simply ggplot2 objects, expert users can also manipulate the final output to customize figures up to their desire.

For now the package contains a sole plotting function, `fantaxtic_bar`, and some helper functions, but Fantaxtic will be extended with more functionality over time. The function `fantaxtic_bar` creates a bar plot of the relative abundances of selected taxa. Main colors are selected based on a user specified taxonomic rank (for example Phylum), and subcoloring (shades and tints) are generated automatically based on a user specified lower taxonomic rank (for example Genus). See the example below for details.

Keywords: _taxonomy, phyloseq, relative, abundance, bar, ggplot2, R_

# Installation
Make sure that `phyloseq` is installed:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
```
Now install `Fantaxtic` with `devtools`:
```
if(!"devtools" %in% installed.packages()){
  install.packages("devtools", dependencies = T, )
}

devtools::install_github("gmteunisse/Fantaxtic")
```

# Loading Fantaxtic
In the R console, type:
```
library(fantaxtic)
```

# Example: Barplot of relative OTU/ASV abundances
First load the data.
```{r, include = T, fig.align = "center", fig.width = 8}
library(fantaxtic)
data(GlobalPatterns)
```

Now subset the data to the 10 most abundants OTUs or ASVs using the `get_top_taxa` function. All other OTUs/ASVs will be collapsed into `"Other"` 
```{r, include = T, fig.align = "center", fig.width = 8}
ps_tmp <- get_top_taxa(physeq_obj = GlobalPatterns, n = 10, relative = TRUE,
                       discard_other = FALSE, other_label = "Other")
```

Next, create labels for missing taxonomic ranks using the `name_taxa` function.
```{r, include = T, fig.align = "center", fig.width = 8}
ps_tmp <- name_taxa(ps_tmp, label = "Unkown", species = T, other_label = "Other")
```

You are now ready to generate a `fantaxtix_bar`! Generate a plot that is colored by Phylum and labeled by Species, coloring the collapsed taxa grey.
The output is a translation table of phylum to color, and a `ggplot2` plot.
```{r, include = T, fig.align = "center", fig.width = 8}
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Species", other_label = "Other")
```

Perhaps a different sample ordering is more meaningful. Let's order the samples by the abundance of the collapsed category, labeled as "Other".
```{r, include = T, fig.align = "center", fig.width = 8}
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Family", other_label = "Other",
               order_alg = "other.abnd")
```

Does the story change when we facet the plot by sample type? Add a faceting factor.
```{r, include = T, fig.align = "center", fig.width = 8, fig.height = 10}
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Family", facet_by = "SampleType", other_label = "Other")
