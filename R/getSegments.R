#' Get Segments Data
#'
#' This function retrieves segments data for either haplotype or dosage from a `genovis` object.
#' The data can then be used for plotting graphical genotypes.
#'
#' @param object A `genovis` class object.
#' @param data A character string specifying which data to use, either "haplotype" or "dosage". Default is "haplotype".
#' @param sample A character or numeric vector specifying samples to include in the retrieval. Default is `NULL`.
#' @return A dataframe containing the segments data.
#' @export
getSegments <- function(object, data = "haplotype", sample = NULL) {
  # Check if the input object is of class 'genovis'
  if (!inherits(x = object, what = "genovis")) {
    stop("The input object should be a genovis class object.")
  }

  # Check if segments data is available
  if (is.null(object$segments)) {
    stop('Run evalSegments().')
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

  # Match the data argument with allowed choices
  data <- match.arg(arg = data, choices = c("haplotype", "dosage"))

  # Process haplotype data if specified
  if (data == "haplotype") {
    out <- object$segments$haplotype
    n_hap <- dim(object$haplotype)[1]
    sample_labels <- paste(rep(object$sample_info$id[sample], each = n_hap),
                           paste0("hap", seq_len(n_hap)), sep = "_")
    legend <- "Haplotype"
    scale_breaks <- attributes(object$haplotype)$scale_breaks
    scale_labels <- attributes(object$haplotype)$scale_labels
    
  } else if (data == "dosage") {
    out <- object$segments$dosage
    sample_labels <- object$sample_info$id[sample]
    legend <- "Dosage"
    scale_breaks <- attributes(object$dosage)$scale_breaks
    scale_labels <- attributes(object$dosage)$scale_labels
  }

  # Filter segments data for the specified samples
  out <- subset(out, subset = name %in% sample_labels)
  out$class <- factor(out$class, levels = scale_breaks)
  out$name <- factor(out$name, levels = sample_labels)

  return(out)
}
