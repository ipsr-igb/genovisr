#'
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#'

plotHist <- function(object,
                     data = "haplotype",
                     sample = NULL,
                     binwidth = 1e6,
                     fill = "darkgreen",
                     chrwise = FALSE,
                     samplewise = FALSE,
                     ncol = NULL){
  if(!inherits(x = object, what = "genovis")){
    stop("The input object should be a genovis class object.")
  }

  if(is.null(object$segments)){
    stop('Run evalSegments() to prepare data to plot graphical genotypes.')
  }

  if(is.null(object$stats)){
    stop('Run statsGeno() to prepare data to plot graphical genotypes.')
  }

  df <- getSegments(object = object, data = data, sample = sample)

  p <- ggplot(data = df) +
    geom_histogram(mapping = aes(x = segment_len * 1e-6), fill = fill,
                   binwidth = binwidth * 1e-6) +
    xlab("Segment length (Mb)") +
    ylab("Count")

  if(chrwise){
    if(is.null(ncol)){
      p <- p + facet_wrap(facets = ~ chr)

    } else {
      p <- p + facet_wrap(facets = ~ chr, ncol = ncol)
    }
  } else if(samplewise){
    p <- p + facet_wrap(facets = ~ name)
  }

  return(p)
}
