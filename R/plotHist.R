#' Plot Histogram of Segment Lengths
#'
#' This function plots a histogram of segment lengths for either haplotype or dosage data from a `genovis` object.
#'
#' @param object A `genovis` class object.
#' @param data A character string specifying which data to use, either "haplotype" or "dosage". Default is "haplotype".
#' @param sample A character or numeric vector specifying samples to include in the plot. Default is `NULL`.
#' @param binwidth A numeric value specifying the bin width for the histogram. Default is 1e6 (1 Mb).
#' @param fill A character string specifying the fill color for the histogram bars. Default is "darkgreen".
#' @param chrwise A logical value indicating whether to plot histograms chromosome-wise. Default is `FALSE`.
#' @param samplewise A logical value indicating whether to plot histograms sample-wise. Default is `FALSE`.
#' @param ncol A numeric value specifying the number of columns for facet wrapping. Default is `NULL`.
#' @return A ggplot object representing the histogram of segment lengths.
#' @export
#' @import ggplot2
plotHist <- function(object,
                     data = "haplotype",
                     sample = NULL,
                     binwidth = 1e6,
                     fill = "darkgreen",
                     chrwise = FALSE,
                     samplewise = FALSE,
                     ncol = NULL) {
  # Check if the input object is of class 'genovis'
  if (!inherits(x = object, what = "genovis")) {
    stop("The input object should be a genovis class object.")
  }

  # Check if segments data is available
  if (is.null(object$segments)) {
    stop('Run evalSegments() to prepare data to plot graphical genotypes.')
  }

  # Check if stats data is available
  if (is.null(object$stats)) {
    stop('Run statsGeno() to prepare data to plot graphical genotypes.')
  }

  # Get segments data
  df <- getSegments(object = object, data = data, sample = sample)
  df$segment_len <- df$segment_len * 1e-6  # Convert segment length to Mb

  # Create histogram plot
  p <- ggplot(data = df) +
    geom_histogram(mapping = aes(x = segment_len), fill = fill,
                   binwidth = binwidth * 1e-6) +
    xlab("Segment length (Mb)") +
    ylab("Count")

  # Add facets if specified
  if (chrwise) {
    if (is.null(ncol)) {
      p <- p + facet_wrap(facets = ~ chr)
    } else {
      p <- p + facet_wrap(facets = ~ chr, ncol = ncol)
    }
  } else if (samplewise) {
    p <- p + facet_wrap(facets = ~ name)
  }

  return(p)
}
