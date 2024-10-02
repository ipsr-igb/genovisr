#' Plot Graphical Genotypes
#'
#' This function plots graphical genotypes for either haplotype or dosage data from a `genovis` object.
#'
#' @param object A `genovis` class object.
#' @param data A character string specifying which data to use, either "haplotype" or "dosage". Default is "haplotype".
#' @param sample A character or numeric vector specifying samples to include in the plot. Default is `NULL`.
#' @param direction A character string specifying the plot direction, either "h" (horizontal) or "v" (vertical). Default is "h".
#' @param width A numeric value specifying the width of the plot elements. Default is 0.1.
#' @return A ggplot object representing the graphical genotypes.
#' @export
#' @import ggplot2
plotGraphGeno <- function(object,
                          data = "haplotype",
                          sample = NULL,
                          direction = "h",
                          margin = 0.1) {
  # Check if the input object is of class 'genovis'
  if (!inherits(x = object, what = "genovis")) {
    stop("The input object should be a genovis class object.")
  }

  # Check if segments data is available
  if (is.null(object$segments)) {
    stop('Run evalSegments() to prepare data to plot graphical genotypes.')
  }

  # If no samples are specified, select all samples
  if (is.null(sample)) {
    sample <- rep(TRUE, nrow(object$sample_info))
  } else {
    # Determine samples based on the type of input
    if (is.character(sample)) {
      sample <- object$sample_info$id %in% sample
    } else if (is.numeric(sample)) {
      sample <- seq_along(object$sample_info$id) %in% sample
    }
  }

  # Limit to one sample if plotting vertically
  if (direction == "v" & sum(sample) > 1) {
    message("When 'v' was specified to 'direction', ",
            "a graphical genotype plot for only one sample can be drawn.")
    message("The first sample was chosen for plotting.")
    sample[-which(sample)[1]] <- FALSE
  }

  # Match the data argument with allowed choices
  data <- match.arg(arg = data, choices = c("haplotype", "dosage"))

  # Process haplotype data if specified
  if (data == "haplotype") {
    df <- object$segments$haplotype
    n_hap <- dim(object$haplotype)[1]
    sample_labels <- paste(rep(object$sample_info$id[sample], each = n_hap),
                           paste0("hap", seq_len(n_hap)), sep = "_")
    legend <- "Haplotype"
    scale_breaks <- attributes(object$haplotype)$scale_breaks
    scale_labels <- attributes(object$haplotype)$scale_labels

  } else if (data == "dosage") {   # Process dosage data if specified
    df <- object$segments$dosage
    sample_labels <- object$sample_info$id[sample]
    legend <- "Dosage"
    scale_breaks <- attributes(object$dosage)$scale_breaks
    scale_labels <- attributes(object$dosage)$scale_labels
  }

  # Filter segments data for the specified samples
  df <- subset(df, subset = name %in% sample_labels)
  df$class <- factor(df$class, levels = scale_breaks)
  df$name <- factor(df$name, levels = sample_labels)

  # Plot based on the specified direction
  if (direction == "v") {
    p <- .plotVertical(df = df,
                       legend = legend,
                       scale_breaks = scale_breaks,
                       scale_labels = scale_labels,
                       sample_name = object$sample_info$id[sample],
                       margin = margin)

  } else {
    p <- .plotHorizontal(df = df,
                         legend = legend,
                         scale_breaks = scale_breaks,
                         scale_labels = scale_labels,
                         margin = margin)
  }

  return(p)
}

#' Plot Vertical Graphical Genotypes
#'
#' @param df A dataframe containing the segments data.
#' @param legend A character string for the legend title.
#' @param scale_breaks A numeric vector for scale breaks.
#' @param scale_labels A character vector for scale labels.
#' @param sample_name A character string for the sample name.
#' @param margin A numeric value specifying the margin of the plot elements.
#' @return A ggplot object representing the vertical graphical genotypes.
.plotVertical <- function(df, legend, scale_breaks, scale_labels,
                          sample_name, margin) {
  df$ymin <- df$start_pos * 1e-6
  df$ymax <- df$end_pos * 1e-6
  name_index <- as.numeric(df$name)
  col_index <- max(name_index) * max(df$chr)
  col_val <- matrix(data = seq_len(col_index), nrow = max(name_index))
  col_val <- col_val - 1
  col_val <- col_val + margin * col_val
  col_val <- col_val + rep(c(0, margin * seq_len(max(df$chr) - 1)), each = max(name_index))
  df$xmin <- sapply(seq_len(nrow(df)), function(i){
    return(col_val[name_index[i], df$chr[i]])
  })
  df$xmax <- df$xmin + 1
  df$xmid <- (df$xmax + df$xmin) / 2
  x_breaks <- sort(unique(df$xmid))
  x_breaks <- matrix(data = x_breaks, nrow = max(name_index))
  x_breaks <- apply(x_breaks, 2, median)
  x_labels <- unique(df$chr)

  max_pos <- max(df$ymax)
  pow <- as.integer(max_pos / 50)
  by <- 5 * 2^pow
  major_breaks <- seq(0, max_pos, by)
  major_breaks <- c(major_breaks, max(major_breaks) + by)
  minor_breaks <- seq(0, max(major_breaks), pow + 1)
  minor_breaks <- minor_breaks[!minor_breaks %in% major_breaks]

  hit <- match(df$class, scale_breaks)
  df$Class <- scale_labels[hit]

  p <- ggplot(df, aes(text = paste('</br>ID: ', name,
                                   '</br>Chromosome: ', chr,
                                   '</br>Start pos (bp): ', start_pos,
                                   '</br>End pos (bp): ', end_pos))) +
    geom_rect(aes(ymin = ymin, ymax = ymax,
                  xmin = xmin, xmax = xmax, fill = Class),
              color = "black", size = 0.2) +
    scale_y_reverse(breaks = major_breaks, minor_breaks = minor_breaks) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, position = "top") +
    scale_fill_viridis_d(name = legend) +
    labs(title = sample_name) +
    ylab("Physical position (Mb)") +
    xlab("Chromosome") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x.top = element_text(margin = margin(b = -rel(10), unit = "pt"), size = rel(1.2)),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour = "gray20", size = 0.2, linetype = 1),
          panel.grid.minor.y = element_line(colour = "gray50", size = 0.2, linetype = 1),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x.top = element_text(vjust = -2),
          legend.position = "top",
          panel.spacing = unit(0, "lines"))
  return(p)
}

#' Plot Horizontal Graphical Genotypes
#'
#' @param df A dataframe containing the segments data.
#' @param legend A character string for the legend title.
#' @param scale_breaks A numeric vector for scale breaks.
#' @param scale_labels A character vector for scale labels.
#' @return A ggplot object representing the horizontal graphical genotypes.
.plotHorizontal <- function(df, legend, scale_breaks, scale_labels, margin) {
  hit <- match(df$class, scale_breaks)
  df$Class <- scale_labels[hit]
  row_val <- rev(as.numeric(df$name)) - 1
  col_val <- as.numeric(df$chr)
  df$ymin <- row_val + margin * row_val
  df$ymax <- df$ymin + 1
  df$ymid <- (df$ymax + df$ymin) / 2
  y_breaks <- unique(df$ymid)
  y_labels <- unique(df$name)

  max_x <- tapply(df$end_pos, df$chr, max)
  x_offset <- c(0, head(cumsum(max_x), -1)) + (seq_along(max_x) - 1) * margin * 3e7
  df$xmin <- (df$start_pos + x_offset[col_val]) * 1e-6
  df$xmax <- (df$end_pos + x_offset[col_val]) * 1e-6
  x_breaks <- (x_offset + x_offset + max_x) / 2 * 1e-6
  x_labels <- names(max_x)

  p <- ggplot(df, aes(text = paste('</br>ID: ', name,
                                   '</br>Chromosome: ', chr,
                                   '</br>Start pos (bp): ', start_pos,
                                   '</br>End pos (bp): ', end_pos))) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax, fill = Class), color = "gray", linewidth = 0.3) +
    scale_fill_viridis_d(name = legend) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, position = "top") +
    scale_y_continuous(breaks = y_breaks, labels = y_labels) +
    ylab("") +
    xlab("Chromosome") +
    theme(axis.ticks = element_blank(),
          axis.text.y.left = element_text(margin = margin(r = -rel(15), unit = "pt"), size = rel(1.2)),
          axis.text.x.top = element_text(margin = margin(b = -rel(10), unit = "pt"), size = rel(1.2)),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "top",
          panel.spacing = unit(0, "lines"))

  return(p)
}
