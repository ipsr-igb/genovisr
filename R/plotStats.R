#'
#'
#'
#'
#' @export
#'
plotStats <- function(object = object,
                      data = c("haplotype", "dosage", "recomb"),
                      group = c("hap", "sample", "chr", "genome"),
                      value = c("class", "block")){
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
                     choices = c("class", "block"),
                     several.ok = FALSE)

  if(data == "recomb"){
    if(group == "genome"){
      stop("Whole-genome summary for recombination frequency is not available.")
    }
    df <- object$stats$haplotype[[paste(group, "recomb", sep = "_")]]
    mean_recomb <- mean(df$value, na.rm = TRUE)
    median_recomb <- median(df$value, na.rm = TRUE)
    p <- ggplot(data = df) +
      geom_bar(mapping = aes(x = name, y = value), fill = "darkgreen",
               stat = "identity") +
      geom_hline(yintercept = mean_recomb, color = "magenta") +
      geom_hline(yintercept = median_recomb, color = "blue") +
      xlab("") +
      ylab("Number of recombinations")

  } else {
    df <- object$stats[[data]][[paste(group, value, sep = "_")]]

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

      hit <- match(df$class, scale_breaks)
      df$class <- scale_labels[hit]

      p <- ggplot(data = df) +
        geom_bar(mapping = aes(x = name, y = value, fill = factor(class)),
                 stat = "identity") +
        scale_fill_viridis_d(name = legend) +
        xlab("") +
        ylab("Proportion")

    } else if(value == "block"){
      p <- ggplot(data = df)  +
        geom_bar(mapping = aes(x = name, y = value, fill = factor(stats)),
                 stat = "identity") +
        scale_fill_viridis_d() +
        facet_grid(rows = vars(stats), scales = "free_y") +
        xlab("") +
        ylab("Stats of block lengths") +
        theme(legend.position = "none")
    }
  }

  return(p)
}
