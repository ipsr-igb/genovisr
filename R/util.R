#' Print a summary of a genovis object
#'
#' This function print a summary of a `genovis` object.
#'
#' @param object A `genovis` class object.
#' @export
print.genovis <- function(object){
  message("Markers: ", nrow(object$marker_info))
  message("Samples: ", nrow(object$sample_info))
  if(is.null(object$haplotype)){
    message("Haplotype classes: no haplotype information")
  } else {
    message("Haplotype classes: ", paste(attributes(object$haplotype)$scale_labels, collapse = " "))
  }
  if(is.null(object$dosage)){
    message("Dosage classes: no dosage information")
  } else {
    message("Dosage classes: ", paste(attributes(object$dosage)$scale_labels, collapse = " "))
  }
}
