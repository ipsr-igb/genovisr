#'
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#'

getSegments <- function(object,
                          data = "haplotype",
                          sample = NULL){
  if(!inherits(x = object, what = "genovis")){
    stop("The input object should be a genovis class object.")
  }

  if(is.null(object$segments)){
    stop('Run evalSegments() to prepare data to plot graphical genotypes.')
  }
  if(is.null(sample)){
    sample <- rep(TRUE, nrow(object$sample_info))

  } else {
    if(is.character(sample)){
      sample <- object$sample_info$id %in% sample

    } else if(is.character(sample)){
      sample <- object$sample_info$id %in% sample

    } else if(is.numeric(sample)){
      sample <- seq_along(object$sample_info$id) %in% sample
    }
  }

  data <- match.arg(arg = data, choices = c("haplotype", "dosage"))

  if(data == "haplotype"){
    out <- object$segments$haplotype
    n_hap <- dim(object$haplotype)[1]
    sample_lables <- paste(rep(object$sample_info$id[sample], each = n_hap),
                           paste0("hap", seq_len(n_hap)), sep = "_")
    legend <- "Haplotype"
    scale_breaks <- attributes(object$haplotype)$scale_breaks
    scale_labels <- attributes(object$haplotype)$scale_labels

  } else if(data == "dosage"){
    out <- object$segments$dosage
    sample_lables <- object$sample_info$id[sample]
    legend <- "Dosage"
    scale_breaks <- attributes(object$dosage)$scale_breaks
    scale_labels <- attributes(object$dosage)$scale_labels
  }
  out <- subset(out, subset = name %in% sample_lables)
  out$class <- factor(out$class, levels = scale_breaks)
  out$name <- factor(out$name, levels = sample_lables)

  return(out)
}
