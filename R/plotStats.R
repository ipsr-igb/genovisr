#'
#'
#'
#'
#' @export
#'
plotStats(object = object,
          sample = 1,
          data = c("haplotype", "dosage"),
          group = c("hapwise", "samplewise", "chrwise", "genome"),
          stats = c("prop", "block" ),
          type = c("histgram", "barplot", "boxplot")){
  if(!inherits(x = object, what = "genovis")){
    stop("The input object should be a genovis class object.")
  }

  if(is.null(object$stats)){
    stop('Run statsGeno() to prepare data to plot graphical genotypes.')
  }

  data <- match.arg(arg = data,
                    choices = c("haplotype", "dosage"),
                    several.ok = FALSE)

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

  stats <- match.arg(arg = stats,
                     choices = c("haplotype", "dosage"),
                     several.ok = FALSE)
  if(data == "haplotype"){
    if()
  }
}
