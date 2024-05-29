#'
#'
#'
#'
#' @export
#'
getStats <- function(object = object,
                      data = c("haplotype", "dosage", "recomb"),
                      group = c("hap", "sample", "chr", "genome"),
                      value = c("class", "segment")){
  if(!inherits(x = object, what = "genovis")){
    stop("The input object should be a genovis class object.")
  }

  if(is.null(object$stats)){
    stop('Run statsGeno() to prepare data to plot graphical genotypes.')
  }

  data <- match.arg(arg = data,
                    choices = c("haplotype", "dosage", "recomb"),
                    several.ok = FALSE)

  group <- match.arg(arg = group,
                     choices = c("hap", "sample", "chr", "genome"),
                     several.ok = FALSE)

  value <- match.arg(arg = value,
                     choices = c("class", "segment"),
                     several.ok = FALSE)

  if(data == "recomb"){
    if(group == "genome"){
      stop("Whole-genome summary for recombination frequency is not available.")
    }
    out <- object$stats$haplotype[[paste(group, "recomb", sep = "_")]]

  } else {
    out <- object$stats[[data]][[paste(group, value, sep = "_")]]

    if(value == "class"){
      if(data == "haplotype"){
        legend <- "Haplotype"
        scale_breaks <- attributes(object$haplotype)$scale_breaks
        scale_labels <- attributes(object$haplotype)$scale_labels

      } else if(data == "dosage"){
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
