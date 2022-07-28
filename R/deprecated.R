#' @title Deprecated functions in package \pkg{fantaxtic}.
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality are also mentioned. Help pages for deprecated functions are
#'   available at \code{help("<function>-deprecated")}.
#' @name fantaxtic-deprecated
#' @keywords internal
NULL


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
#'
#' @name get_top_taxa-deprecated
#' @seealso \code{\link{fantaxtic-deprecated}}
#' @usage get_top_taxa(physeq_obj, n, relative = T, discard_other, other_label)
#' @examples
#' data(GlobalPatterns)
#' get_top_taxa(GlobalPatterns, 10)
#' get_top_taxa(GlobalPatterns, 10, relative = FALSE)
#' get_top_taxa(GlobalPatterns, 10, relative = TRUE, discard_other = TRUE)
#' get_top_taxa(GlobalPatterns, 10, relative = TRUE, other_label = "Non-abundant taxa")
NULL

#' @rdname fantaxtic-deprecated
#' @section get_top_taxa:
#' For \code{get_top_taxa}, use \code{\link{top_taxa}} or \code{\link{nested_top_taxa}}.
#'
#' @export
get_top_taxa <- function(physeq_obj, n, relative = TRUE, discard_other = FALSE, other_label = "Other"){

  .Deprecated(c("top_taxa", "nested_top_taxa"))

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

#' Generate a color palette.
#'
#' This function generates \eqn{i_1, i_2, ... i_n} shades and tints specified
#' in \code{clr.tbl} for \eqn{n} colors \code{clr.pal}. If no \code{clr.tbl}
#' is specified, this function generates \eqn{n} maximally separated colors
#' with equal saturation and brightness from a base color.
#'
#' @param clr_tbl A \code{data.frame} or \code{matrix} with two columns: the
#' first column gives the labels to which colors need to be assigned; the
#' second column gives the number of tints and shades that need to be generated
#' for that color.
#' @param clr_pal A color palette with with \eqn{n} colors, where
#' \eqn{n =} \code{nrow(clr_tbl)}.
#' @param base_clr The base color from which to generate a color palette.
#' @return A list of \eqn{n} items, each containing \eqn{i_n} shades.
#'
#' @name gen_palette-deprecated
#' @seealso \code{\link{fantaxtic-deprecated}}
#' @usage gen_palette(clr_tbl, clr_pal = NULL, base_clr = "#6495ed")
#'
#' @examples
#' clr_tbl <- data.frame(Var = c("Phylum A", "Phylum B"),
#'                       Freq = c(5, 8))
#' gen_palette(clr_tbl)
#' gen_palette(clr_tbl, clr_pal = c("red", "blue"))
#' gen_palette(clr_tbl, base_clr = "blue")
NULL

#' @rdname fantaxtic-deprecated
#' @section gen_palette:
#' For \code{gen_palette}, use \code{\link[ggnested]{nested_palette}} in the ggnested package.
#'
#' @export
gen_palette <- function(clr_tbl, clr_pal = NULL, base_clr = "#6495ed"){

  .Deprecated(new = "nested_palette", package = "ggnested")

  #Define a palette
  base_pal <- c("#6495ed", "#ff7256", "#edbc64", "#8470ff", "#8ee5ee", "#EE8DD6")

  #Define the number of variants required per color
  clr_vars <- clr_tbl[,2]
  n_clrs <- length(clr_vars)

  #Check arguments and generate colors if required
  if (is.null(clr_pal)){
    if (base_clr != "#6495ed"){
      clr_pal <- gen_colors(n_clrs, base_clr)
    } else {
      if (n_clrs <= length(base_pal)){
        clr_pal <- base_pal[1:n_clrs]
      } else {
        clr_pal <- gen_colors(n_clrs, base_clr)
      }
    }
  } else {
    if (length(clr_pal) < n_clrs){
      stop(sprintf("Error: %d values required in clr.pal, %d provided.", n_clrs, length(clr_pal)))
    }
    clr_pal <- clr_pal[1:n_clrs]
  }

  #Create the translation table from color.by to their central colors
  clr_mtrx <- cbind(clr_tbl, clr_pal)
  names(clr_mtrx) <- c("Level", "N.color.shades", "Central.color")
  clr_mat <- clr_mtrx

  #Generate i vars for each clr in  clr.pal
  clr_mtrx <- cbind(clr_vars, clr_pal)
  clr_list <- split(clr_mtrx, seq(nrow(clr_mtrx)))
  palette <- lapply(clr_list, function(v){
    i <- as.numeric(v[1])
    clr <- v[2]
    gen_shades_tints(i, clr)
  })
  return(list(palette = palette,
              colours = clr_mat))
}

#' Generate a barplot of relative taxon abundances
#'
#' This function generates a \code{ggplot2} barplot for the relative abundances
#' of taxa in a phyloseq object.
#'
#' Coloring occurs by a user specified taxonomic
#' level, and subcoloring according to another level can be added if desired (e.g.
#' color by phylum, subcolor by genus). In addition, one or more taxon names can
#' be specified as "other" that receive a specific color (e.g. outliers,
#' collapsed taxa).
#'
#' By default, unique labels per taxon will be generated in the case that
#' multiple taxa with identical labels exist, unless the user chooses to suppress
#' this. Moreover, \code{NA} values in the \code{tax_table} of the phyloseq
#' object will be renamed to \code{"Unknown"} to avoid confusion. WARNING:
#' duplicate labels in the data lead to incorrect displaying of data and
#' labels.
#'
#' To facilitate visualisation and/or interpretation, samples can be reordered
#' according alphabetically, by the abundance of a certain taxon, by
#' hierarhical clustering or by the abundance of "other" taxa.
#'
#' @param physeq_obj A phyloseq object with an \code{otu_table}, a \code{tax_table}
#' and, in case of facetting, \code{sample_data}.
#' @param color_by The name of the taxonomic level by which to color the bars.
#' @param label_by The name of the taxonomic level by which to label the bars and
#' generate subcolors.
#' @param facet_by The name of the factor in the \code{sample_data} by which to
#' facet the plots.
#' @param grid_by The name of a second factor in the \code{sample_data} by which to
#' facet to plots, resulting in a grid.
#' @param facet_type The type of faceting from ggplot2 to use, either \code{grid}
#' or \code{wrap} (default).
#' @param facet_cols The number of columns to use for faceting.
#' @param gen_uniq_lbls Generate unique labels (default = \code{TRUE})?
#' @param other_label A character vector specifying the names of taxa in
#' \code{label_by} to use a specific color for.
#' @param order_alg The algorithm by which to order samples, or one or more taxa
#' found in \code{label_by}. Algorithms can be one of \code{hclust}
#' (hierarhical clustering; default), \code{as.is} (current order) or \code{alph}
#' (alphabetical).
#' @param color_levels Character vector containing names of levels. Useful to
#' enforce identical colors for levels across different plots or to pair levels
#' with colors.
#' @param base_color The base color from which to generate colors.
#' @param other_color The base color from which to generate shades for "other"
#' taxa.
#' @param palette A user specified palette to color the bars with. Replaces
#' \code{base_color}.
#' @param bar_width The width of the bars as a fraction of the binwidth
#' (default = 0.9).
#' @return A \code{ggplot2} object.
#' @import ggplot2 phyloseq reshape2
#'
#' @name fantaxtic_bar-deprecated
#' @seealso \code{\link{fantaxtic-deprecated}}
#' @usage fantaxtic_bar(physeq_obj, color_by, label_by = NULL, facet_by = NULL, grid_by = NULL, facet_type = "wrap", bar_width = 0.9, facet_cols =  1, gen_uniq_lbls = TRUE, other_label= NULL, order_alg = "hclust", color_levels = NULL, base_color = "#6495ed", other_color = "#f3f3f3", palette = NULL)
#'
#' @examples
#' #Load data
#' data(GlobalPatterns)
#'
#' #Get the 10 most abundant OTUs / ASVs
#' ps_tmp <- get_top_taxa(physeq_obj = GlobalPatterns, n = 10, relative = TRUE,
#'                        discard_other = FALSE, other_label = "Other")
#'
#' #Create labels for missing taxonomic ranks
#' ps_tmp <- name_taxa(ps_tmp, label = "Unkown", species = T, other_label = "Other")
#'
#' #Generate a barplot that is colored by Phylum and labeled by Species, coloring
#' #collapsed taxa grey.
#' fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Species", other_label = "Other")
#'
#' #Generate a barplot that is colored by Phylum and lebeled by Species. As multiple
#' ASVs have the same family annotation, generate unique labels.
#' fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Family", other_label = "Other",
#'               gen_uniq_lbls = TRUE)
#'
#' #Change the sample ordering
#' fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Family", other_label = "Other",
#'               order_alg = "other.abnd")
#'
#' #Add faceting by sample type
#' fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Family",
#'               facet_by = "SampleType", other_label = "Other")
NULL

#' @rdname fantaxtic-deprecated
#' @section fantaxtic_bar:
#' For \code{fantaxtic_bar}, use \code{\link{fantaxtic}}.
#'
#' @export
fantaxtic_bar <- function(physeq_obj, color_by, label_by = NULL, facet_by = NULL,
                          grid_by = NULL, facet_type = "wrap", bar_width = 0.9,
                          facet_cols =  1, gen_uniq_lbls = TRUE, other_label= NULL,
                          order_alg = "hclust", color_levels = NULL,
                          base_color = "#6495ed",
                          other_color = "#f3f3f3", palette = NULL){

  .Deprecated(new = "fantaxtic")

  #Check for subcoloring
  if (is.null(label_by)){
    label_by <- color_by
  }

  #Extract tax_tbl and add OTU names
  tax_tbl <- as.data.frame(phyloseq::tax_table(physeq_obj))
  tax_tbl$otu_name <- row.names(tax_tbl)

  #Replace NAs with Unknown
  tax_tbl <- as.data.frame(apply(tax_tbl, 2, function(x){
    x[is.na(x)] <- "Unknown"
    return(x)
  }))

  #Move Other taxa to the beginning and alter taxonomic annotations
  #of Other taxa
  if(!is.null(other_label)){
    main_ind <- which(!tax_tbl[[label_by]] %in% other_label)
    other_ind <- which(tax_tbl[[label_by]] %in% other_label)
    new_color_by <- as.character(tax_tbl[[color_by]])
    new_color_by[other_ind] <- as.character(tax_tbl[[label_by]][other_ind])
    tax_tbl[[color_by]] <- as.factor(new_color_by)
    ordr <- c(other_ind, main_ind)
    tax_tbl <- tax_tbl[ordr,]
  }

  #Refactor for legend ordering and order
  if (is.null(color_levels)){
    tax_levels <- unique(tax_tbl[[color_by]])
  } else {
    if (is.null(other_label)){
      tax_levels <- color_levels
    } else {
      tax_levels <- c(other_label, color_levels)
    }
  }
  tax_tbl[[label_by]] <- factor(tax_tbl[[label_by]], unique(tax_tbl[[label_by]]), ordered = T)
  tax_tbl[[color_by]] <- factor(tax_tbl[[color_by]], tax_levels, ordered = T)
  tax_tbl <- tax_tbl[order(tax_tbl[[color_by]]),]

  #Get the tax and OTU tables
  otu_tbl <- as.data.frame(phyloseq::otu_table(physeq_obj))

  #Check the orientation of the otu_tbl and change if required
  if (!taxa_are_rows(phyloseq::otu_table(physeq_obj))){
    otu_tbl <- as.data.frame(t(otu_tbl))
  }

  #Calculate the number of colors and color variations required
  clr_tbl <- as.data.frame(table(tax_tbl[[color_by]], useNA = "ifany"), stringsAsFactors = F)
  if(!is.null(other_label)){
    ind <- which(!clr_tbl$Var1 %in% other_label)
    clr_tbl <- clr_tbl[ind,]
  }

  #Generate the required color palette
  clr_pal <- gen_palette(clr_tbl = clr_tbl, clr_pal = palette, base_clr = base_color)$palette
  names(clr_pal) <- clr_tbl$Var1
  clr_pal <- as.vector(unlist(clr_pal))
  if(!is.null(other_label)){
    n_other <- length(other_label)
    other_pal <- gen_shades_tints(n_other, clr = other_color)
    clr_pal <- c(other_pal, clr_pal)
  }

  #Generate unique label names if required
  if(gen_uniq_lbls){
    tax_tbl[[label_by]] <- gen_uniq_lbls(tax_tbl[[label_by]])
  }

  #Transform absolute taxon counts to relative values
  otu_tbl <- as.data.frame(apply(otu_tbl, 2, function(x){
    if (sum(x) > 0){x/sum(x)}
    else(x)
  }))

  #Match order of tax.tbl and otu.tbl
  ord <- match(tax_tbl$otu_name, row.names(otu_tbl))
  otu_tbl <- otu_tbl[ord,]

  #Order the samples according to the specified algorithm
  #Order according to selected taxonomies
  if (sum(order_alg %in% c("alph", "hclust", "as.is")) == 0){

    #Get the summed abundances
    sums <- list()
    i <- 0
    for (lvl in order_alg){
      i <- i + 1
      sums[[i]] <- round(colSums(otu_tbl[which(tax_tbl[[label_by]] == lvl),]), digits = 3)
    }

    #Sort
    cmd <- paste(sprintf("sums[[%d]]", 1:i), collapse = ", ")
    smpl_ord <- eval(parse(text = sprintf("order(%s)", cmd)))
    otu_tbl <- otu_tbl[,smpl_ord]

    #Order according to selected algorithm
  }else{
    if (order_alg == "alph"){
      otu_tbl <- otu_tbl[,order(names(otu_tbl))]
    } else {
      if(order_alg == "hclust"){
        hc <- hclust(dist(x = t(otu_tbl), method = "euclidian", upper = F))
        smpl_ord <- hc$order
        otu_tbl <- otu_tbl[,smpl_ord]
      } else {
        if (order_alg == "as.is"){
          #do nothing
        }
      }
    }
  }


  #Join labels and counts and transform to a long data format
  counts <- cbind(tax_tbl[[color_by]], tax_tbl[[label_by]], otu_tbl)
  names(counts) <- c("color_by", "label_by", colnames(otu_tbl))
  counts_long <- reshape2::melt(counts,
                                id.vars = c("color_by", "label_by"),
                                variable.name = "Sample",
                                value.name = "Abundance")

  #Add facet levels if needed and transform to a long data format
  if(is.null(facet_by) & !is.null(grid_by)){
    facet_by <- grid_by
    grid_by <- NULL
  }
  if (!is.null(facet_by)){
    facet <- as.data.frame(phyloseq::sample_data(physeq_obj))[[facet_by]]
    names(facet) <- row.names(phyloseq::sample_data(physeq_obj))
    ord <- match(counts_long$Sample, names(facet))
    facet <- facet[ord]
    counts_long$facet <- facet
  }
  if (!is.null(grid_by)){
    grid <- as.data.frame(phyloseq::sample_data(physeq_obj))[[grid_by]]
    names(grid) <- row.names(phyloseq::sample_data(physeq_obj))
    ord <- match(counts_long$Sample, names(grid))
    grid <- grid[ord]
    counts_long$grid <- grid
  }

  #Generate a plot
  p <- ggplot2::ggplot(counts_long, aes(x = Sample, y = Abundance, fill = label_by)) +
    ggplot2::geom_bar(position = "stack", stat = "identity", width = bar_width) +
    ggplot2::guides(fill=guide_legend(title = label_by, ncol = 1)) +
    ggplot2::scale_fill_manual(values = clr_pal) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::theme(axis.line.x = element_line(colour = 'grey'),
                   axis.line.y = element_line(colour = 'grey'),
                   axis.ticks = element_line(colour = 'grey'),
                   axis.text.x = element_text(angle = 90, family = "Helvetica",
                                              size = 6, hjust = 1, vjust = 0.5),
                   legend.background = element_rect(fill = 'transparent', colour = NA),
                   legend.key = element_rect(fill = "transparent"),
                   legend.key.size = unit(0.4, "cm"),
                   panel.background = element_rect(fill = 'transparent', colour = NA),
                   panel.grid.major.x = element_blank(),
                   panel.grid.major.y = element_line(colour = adjustcolor('grey', 0.2)),
                   panel.grid.minor = element_line(colour = NA),
                   plot.background = element_rect(fill = 'transparent', colour = NA),
                   plot.title = element_text(hjust = 0.5),
                   strip.background = element_blank(),
                   strip.text = element_text(family = "Helvetica", size = 8, face = "bold"),
                   text = element_text(family = "Helvetica", size = 8))

  if (!is.null(facet_by)) {
    if (facet_type == "wrap"){
      if (is.null(grid_by)){
        p <- p + ggplot2::facet_wrap(~facet, scales = "free", ncol = facet_cols)
      }else{
        p <- p + ggplot2::facet_wrap(~grid + facet, scales = "free", ncol = facet_cols)
      }
    }else{
      if (facet_type == "grid"){
        if (is.null(grid_by)){
          p <- p + ggplot2::facet_grid(~facet, scales = "free", space = "free")
        }else{
          p <- p + ggplot2::facet_grid(facet ~ grid, scales = "free", space = "free")
        }
      }
    }
  }

  return(p)
}

#' Generate a color palette.
#'
#' This function generates \code{n} maximally separated colors with equal
#' saturation and brightness from a base color.
#'
#' @param n The number of colors to generate.
#' @param clr The base color from which to generate other colors.
#' @return A vector of \code{n} colors.
#'
#' @name gen_colors-deprecated
#' @seealso \code{\link{fantaxtic-deprecated}}
#' @usage gen_colors(n, clr = "#6495ed")
#'
#' @examples
#' gen_colors(5)
#' gen_colors(5, "blue")
NULL

#' @rdname fantaxtic-deprecated
#' @section gen_colors:
#' For \code{gen_colors}, use \code{\link[ggnested]{nested_palette}} in the ggnested package.
#'
#' @export
gen_colors <- function(n, clr = "#6495ed"){

  .Deprecated(new = "nested_palette", package = "ggnested")
  if (n == 1){
    return(clr)
  }
  clr.rgb <- t(col2rgb(clr))
  clr.hsv <- t(rgb2hsv(clr.rgb[1], clr.rgb[2], clr.rgb[3], 255))
  clrs <- 0:(n-1)
  clrs <- sapply(clrs, function(x){
    offset <- x/(n)
    h <- (clr.hsv[1] + offset) %% 1
    s <- clr.hsv[2]
    v <- clr.hsv[3]
    clr <- hsv(h, s, v)
    return(clr)
  })
  clrs <- shuffle_colors(clrs)
  return(clrs)
}

#' Reshuffle a list of colors.
#'
#' This function reshuffles a vector of colours so that in a list of \code{n}
#' colors, every color \code{i} in \code{\{1, n / 2\}} will be placed next
#' to color \code{i + n / 2}. This is useful for visualization of
#' automatically generated color palettes based upon for example gradients
#' or HSL color circles.
#'
#' @param clrs A vector of colors to reshuffle.
#' @return A reshuffled vector of \code{clrs}.
#'
#' @name shuffle_colors-deprecated
#' @seealso \code{\link{fantaxtic-deprecated}}
#' @usage shuffle_colors(clrs)
#'
#' @examples
#' shuffle_colors(gen_colors(6))
NULL

#' @rdname fantaxtic-deprecated
#' @section shuffle_colors:
#' For \code{shuffle_colors}, use \code{\link[ggnested]{nested_palette}} in the ggnested package.
#'
#' @export
shuffle_colors <- function(clrs){

  .Deprecated(new = "nested_palette", package = "ggnested")

  n.col <- length(clrs)
  ordr <- rep(1:ceiling(n.col/2), each = 2, length.out = n.col)
  indx <- seq(2, length(ordr), 2)
  ordr[indx] <- ordr[indx] + ceiling(n.col/2)
  clrs <- clrs[ordr]
  return(clrs)
}

#' Generate shade and tint variants for a base color.
#'
#' \code{gen_shades_tints} returns \code{i/2} shades and \code{i/2} tints
#' of a base color.
#'
#' @param i The number of tints and/or shades to be generated.
#' @param clr The base color to use (hex or named color).
#' @param incl.base Should the base color be part of the generated palette?
#' @return \code{i} variants of a base color.
#'
#' @name gen_shades_tints-deprecated
#' @seealso \code{\link{fantaxtic-deprecated}}
#' @usage gen_shades_tints(i, clr = "#6495ed")
#'
#' @examples
#' gen_shades(3)
#' gen_tints(3)
#' gen_shades_tints(3)
#' gen_shades(3, "blue")
#' gen_tints(3, "blue")
#' gen_shades_tints(3, "blue")
NULL

#' @rdname fantaxtic-deprecated
#' @section gen_shades_tints:
#' For \code{gen_shades_tints}, use \code{\link[ggnested]{nested_palette}} in the ggnested package.
#'
#' @export
gen_shades_tints <- function(i, clr = "#6495ed"){

  .Deprecated(new = "nested_palette", package = "ggnested")
  if (i == 0){
    return(c())
  }
  if (i == 1){
    return(clr)
  }
  if (i == 2){
    n_tnts <- 2
    n_shds <- 1
    tnts <- gen_tints(n_tnts, clr, incl.base = F)[1]
  } else {
    n_tnts <- ceiling(i/2)
    n_shds <- i - n_tnts
    tnts <- gen_tints(n_tnts, clr, incl.base = F)
  }
  shds <- gen_shades(n_shds, clr, incl.base = T)
  clrs <- c(shds, tnts)
  return(clrs)
}

gen_shades <- function(i, clr = "#6495ed", incl.base = FALSE){
  if (i == 0){
    return(c())
  }
  if (i == 1){
    return(clr)
  }
  clr_rgb <- t(col2rgb(clr))
  shds <- (1 + round(i)):(2 * round(i) + 1)
  shds <- sapply(shds, function(x){
    shd_rgb <- clr_rgb * x / (2 * round(i) + 1)
    shd <- rgb(shd_rgb[1], shd_rgb[2], shd_rgb[3], maxColorValue = 255)
    return(shd)
  })
  if (incl.base){
    return(shds[1:i+1])
  }
  return(shds[1:i])
}

gen_tints <- function(i, clr = "#6495ed", incl.base = FALSE){
  if (i == 0){
    return(c())
  }
  if (i == 1){
    return(clr)
  }
  clr_rgb <- t(col2rgb(clr))
  tnts <- 0:i
  tnts <- sapply(tnts, function(x){
    tnt_rgb <- clr_rgb + (255 - clr_rgb) * x / (round(i) * 1.5)
    tnt <- rgb(tnt_rgb[1], tnt_rgb[2], tnt_rgb[3], maxColorValue = 255)
    return(tnt)
  })
  if (incl.base){
    return(tnts[1:i])
  }
  return(tnts[1:i+1])
}


