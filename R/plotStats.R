#'
#'
#'
#'
#' @export
#'
plotStats <- function(object = object,
                      data = c("haplotype", "dosage"),
                      group = c("hap", "sample", "chr", "genome"),
                      value = c("class", "block")){
  if(!inherits(x = object, what = "genovis")){
    stop("The input object should be a genovis class object.")
  }

  if(is.null(object$stats)){
    stop('Run statsGeno() to prepare data to plot graphical genotypes.')
  }

  data <- match.arg(arg = data,
                    choices = c("haplotype", "dosage"),
                    several.ok = FALSE)

  group <- match.arg(arg = group,
                     choices = c("hap", "sample", "chr", "genome"),
                     several.ok = FALSE)
  if(data == "haplotype"){
  }
}
