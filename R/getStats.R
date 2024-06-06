#' Get Statistics Data
#'
#' This function retrieves statistics data for haplotype, dosage, or recombination
#' from a `genovis` object. The data can then be used for plotting graphical genotypes.
#'
#' @param object A `genovis` class object.
#' @param data A character string specifying which data to use, either "haplotype", "dosage", or "recomb". Default is "haplotype".
#' @param group A character string specifying the grouping for the statistics, either "hap", "sample", "chr", or "genome". Default is "hap".
#' @param value A character string specifying the value to be retrieved, either "class" or "segment". Default is "class".
#' @return A dataframe containing the statistics data.
#' @export
getStats <- function(object,
                     data = c("haplotype", "dosage", "recomb"),
                     group = c("hap", "sample", "chr", "genome"),
                     value = c("class", "segment")) {
  # Check if the input object is of class 'genovis'
  if (!inherits(x = object, what = "genovis")) {
    stop("The input object should be a genovis class object.")
  }

  # Check if stats data is available
  if (is.null(object$stats)) {
    stop('Run statsGeno().')
  }

  # Match the data argument with allowed choices
  data <- match.arg(arg = data, choices = c("haplotype", "dosage", "recomb"), several.ok = FALSE)

  # Match the group argument with allowed choices
  group <- match.arg(arg = group, choices = c("hap", "sample", "chr", "genome"), several.ok = FALSE)

  # Match the value argument with allowed choices
  value <- match.arg(arg = value, choices = c("class", "segment"), several.ok = FALSE)

  # Retrieve statistics data based on the specified arguments
  if (data == "recomb") {
    if (group == "genome") {
      stop("Whole-genome summary for recombination frequency is not available.")
    }
    out <- object$stats$haplotype[[paste(group, "recomb", sep = "_")]]
  } else {
    out <- object$stats[[data]][[paste(group, value, sep = "_")]]

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

      hit <- match(out$class, scale_breaks)
      out$class <- scale_labels[hit]
    }
  }

  return(out)
}
