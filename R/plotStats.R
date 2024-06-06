#' Plot Statistics for Genovis Data
#'
#' This function plots statistics for haplotype, dosage, or recombination data from a `genovis` object.
#'
#' @param object A `genovis` class object.
#' @param data A character string specifying which data to use, either "haplotype", "dosage", or "recomb". Default is "haplotype".
#' @param group A character string specifying the grouping for the statistics, either "hap", "sample", "chr", or "genome". Default is "hap".
#' @param value A character string specifying the value to be plotted, either "class" or "segment". Default is "class".
#' @return A ggplot object representing the statistics plot.
#' @export
plotStats <- function(object,
                      data = c("haplotype", "dosage", "recomb"),
                      group = c("hap", "sample", "chr", "genome"),
                      value = c("class", "segment")) {
  # Check if the input object is of class 'genovis'
  if (!inherits(x = object, what = "genovis")) {
    stop("The input object should be a genovis class object.")
  }

  # Check if stats data is available
  if (is.null(object$stats)) {
    stop('Run statsGeno() to prepare data to plot graphical genotypes.')
  }

  # Match the data, group, and value arguments with allowed choices
  data <- match.arg(arg = data, choices = c("haplotype", "dosage", "recomb"), several.ok = FALSE)
  group <- match.arg(arg = group, choices = c("hap", "sample", "chr", "genome"), several.ok = FALSE)
  value <- match.arg(arg = value, choices = c("class", "segment"), several.ok = FALSE)

  # Plot recombination data if specified
  if (data == "recomb") {
    if (group == "genome") {
      stop("Whole-genome summary for recombination frequency is not available.")
    }
    df <- object$stats$haplotype[[paste(group, "recomb", sep = "_")]]
    mean_recomb <- mean(df$value, na.rm = TRUE)
    median_recomb <- median(df$value, na.rm = TRUE)
    p <- ggplot(data = df) +
      geom_bar(mapping = aes(x = name, y = value), fill = "darkgreen", stat = "identity") +
      geom_hline(yintercept = c(mean_recomb, median_recomb), color = c("magenta", "blue"), linetype = 1:2) +
      xlab("") +
      ylab("Number of recombinations")

    # Plot haplotype or dosage data if specified
  } else {
    df <- object$stats[[data]][[paste(group, value, sep = "_")]]

    # Handle class value
    if (value == "class") {
      if (data == "haplotype") {
        legend <- "Haplotype"
        scale_breaks <- attributes(object$haplotype)$scale_breaks
        scale_labels <- attributes(object$haplotype)$scale_labels
      } else if (data == "dosage") {
        legend <- "Dosage"
        scale_breaks <- attributes(object$dosage)$scale_breaks
        scale_labels <- attributes(object$dosage)$scale_labels
      }

      hit <- match(df$class, scale_breaks)
      df$class <- scale_labels[hit]

      p <- ggplot(data = df) +
        geom_bar(mapping = aes(x = name, y = value, fill = class), stat = "identity") +
        scale_fill_viridis_d(name = legend) +
        xlab("") +
        ylab("Proportion")

      # Handle segment value
    } else if (value == "segment") {
      p <- ggplot(data = df) +
        geom_bar(mapping = aes(x = name, y = value, fill = stats), stat = "identity") +
        scale_fill_viridis_d() +
        facet_grid(rows = vars(stats), scales = "free_y") +
        xlab("") +
        ylab("Stats of segment lengths") +
        theme(legend.position = "none")
    }
  }

  return(p)
}
