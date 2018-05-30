# Fantaxtic
This is Fantaxtic, a package that contains a set of wrapper functions for the phyloseq and ggplot2 packages in the R software that turns ordinary taxonomic count data into fantaxtic, publication ready plots. Fantaxtic figures come as completely stylized figures with pleasing colors and maximal plotting control. As the output figures are simply ggplot2 objects, expert users can also manipulate the final output to customize figures up to their desire.

For now the package contains a sole plotting function, `fantaxtic_bar`, and some helper functions, but Fantaxtic will be extended with more functionality over time. The function `fantaxtic_bar` creates a bar plot of the relative abundances of selected taxa. Main colors are selected based on a user specified taxonomic rank (for example Phylum), and subcoloring (shades and tints) are generated automatically based on a user specified lower taxonomic rank (for example Genus). See the example below for details.

Keywords: _taxonomy, phyloseq, relative, abundance, bar, ggplot2, R_

# Installation
In the R console, type:
```
#install.packages(devtools)
devtools::install_github("gmteunisse/Fantaxtic")
```

# Loading Fantaxtic
In the R console, type:
```
library(fantaxtic)
```

# An example
As an example, you can create a barplot using the GlobalPatterns data from the Phyloseq package.
```
#Load data
data(GlobalPatterns)

#Get the 10 most abundant OTUs / ASVs
ps_tmp <- get_top_taxa(physeq_obj = GlobalPatterns, n = 10, relative = TRUE,
                       discard_other = FALSE, other_label = "Other")

#Create labels for missing taxonomic ranks
ps_tmp <- name_taxa(ps_tmp, label = "Unkown", species = T, other_label = "Other")

#Generate a barplot that is colored by Phylum and labeled by Species, coloring
#collapsed taxa grey.
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Species", other_label = "Other")
