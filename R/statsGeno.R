#' Calculate Statistics for Genovis Data
#'
#' This function calculates statistics for haplotype and dosage data from a `genovis` object.
#'
#' @param object A `genovis` class object.
#' @param chr_len A named numeric vector specifying the lengths of each chromosome. Default is `NULL`, which will calculate lengths based on marker positions.
#' @return A `genovis` class object with calculated statistics.
#' @export
statsGeno <- function(object, chr_len = NULL) {
  # Check if the input object is of class 'genovis'
  if (!inherits(x = object, what = "genovis")) {
    stop("The input object should be a genovis class object.")
  }

  # Check if segments data is available
  if (is.null(object$segments)) {
    stop('Run evalSegments() to prepare data to plot graphical genotypes.')
  }

  # Calculate chromosome lengths if not provided
  if (is.null(chr_len)) {
    chr_len <- tapply(X = object$marker_info$pos,
                      INDEX = object$marker_info$chr,
                      FUN = max)
  }
  if (is.null(names(chr_len))) {
    names(chr_len) <- seq_along(chr_len)
  }

  object <- .extendBorders(object = object, chr_len = chr_len)

  if (is.null(object$haplotype)) {
    message("Haplotype information is not available for the given dataset.")
    hap_stats <- NULL
  } else {
    hap_stats <- .evalStats(object = object, data = "haplotype", chr_len = chr_len)
  }

  if (is.null(object$dosage)) {
    message("Dosage information is not available for the given dataset.")
    dos_stats <- NULL
  } else {
    dos_stats <- .evalStats(object = object, data = "dosage", chr_len = chr_len)
  }

  object$stats <- list(haplotype = hap_stats, dosage = dos_stats)
  return(object)
}

#' Extend Segment Borders
#'
#' @param object A `genovis` class object.
#' @param chr_len A named numeric vector specifying the lengths of each chromosome.
#' @return A `genovis` class object with extended segment borders.
.extendBorders <- function(object, chr_len) {
  for (data in c("haplotype", "dosage")) {
    if (is.null(object$segments[[data]])) {
      next
    }

    chr_first <- which(!duplicated(object$marker_info$chr))
    to_1 <- object$segments[[data]]$start_index %in% chr_first
    object$segments[[data]]$start_border <- NA
    object$segments[[data]]$start_border[to_1] <- 1

    chr_end <- c(chr_first[-1] - 1, nrow(object$marker_info))
    to_end <- object$segments[[data]]$end_index %in% chr_end
    hit <- match(object$segments[[data]]$chr[to_end], names(chr_len))
    object$segments[[data]]$end_border <- NA
    object$segments[[data]]$end_border[to_end] <- chr_len[hit]

    na_start_border <- is.na(object$segments[[data]]$start_border)
    pos1_index <- object$segments[[data]]$start_index[na_start_border] - 1
    pos1 <- object$marker_info$pos[pos1_index]
    pos2 <- object$segments[[data]]$start_pos[na_start_border]
    mid_pos <- round((pos2 - pos1) / 2 + pos1)
    object$segments[[data]]$start_border[na_start_border] <- mid_pos

    na_end_border <- is.na(object$segments[[data]]$end_border)
    pos2_index <- object$segments[[data]]$end_index[na_end_border] + 1
    pos1 <- object$segments[[data]]$end_pos[na_end_border]
    pos2 <- object$marker_info$pos[pos2_index]
    mid_pos <- round((pos2 - pos1) / 2 + pos1)
    object$segments[[data]]$end_border[na_end_border] <- mid_pos
    object$segments[[data]]$segment_len <- object$segments[[data]]$end_border - object$segments[[data]]$start_border + 1

    check <- any(object$segments[[data]]$segment_len < 0)
    if (check) {
      stop('Some segment(s) have lengths of less than zero.',
           "\nPlease confirm the values specified to 'chr_len'.")
    }
  }
  return(object)
}

#' Evaluate Statistics
#'
#' @param object A `genovis` class object.
#' @param data A character string specifying which data to use, either "haplotype" or "dosage".
#' @param chr_len A named numeric vector specifying the lengths of each chromosome.
#' @return A list of evaluated statistics.
.evalStats <- function(object, data, chr_len) {
  out1 <- .evalClassProp(object = object, data = data, chr_len = chr_len)
  out2 <- .evalBlockLen(object = object, data = data, chr_len = chr_len)

  out <- c(out1, out2)
  return(out)
}

#' Evaluate Class Proportions
#'
#' @param object A `genovis` class object.
#' @param data A character string specifying which data to use, either "haplotype" or "dosage".
#' @param chr_len A named numeric vector specifying the lengths of each chromosome.
#' @return A list of class proportions.
.evalClassProp <- function(object, data, chr_len){
  classes <- attributes(object[[data]])$scale_breaks
  target_data <- object$segments[[data]]

  # Haploid-wise class ratio summary
  hap_class_prop <- aggregate(x = segment_len ~ name + class, FUN = sum,
                              drop = FALSE, data = target_data)
  hap_class_prop$segment_len[is.na(hap_class_prop$segment_len)] <- 0
  hap_total_len <- aggregate(x = segment_len ~ name, FUN = sum,
                             drop = FALSE, data = hap_class_prop)
  hit <- match(hap_class_prop$name, hap_total_len$name)
  hap_class_prop$value <- hap_class_prop$segment_len / hap_total_len$segment_len[hit]

  # Sample-wise class ratio summary
  sample_class_prop <- aggregate(x = segment_len ~ sample_id + class, FUN = sum,
                                 drop = FALSE, data = target_data)
  sample_class_prop$segment_len[is.na(sample_class_prop$segment_len)] <- 0
  sample_total_len <- aggregate(x = segment_len ~ sample_id, FUN = sum,
                                drop = FALSE, data = sample_class_prop)
  hit <- match(sample_class_prop$sample_id, sample_total_len$sample_id)
  sample_class_prop$value <- sample_class_prop$segment_len / sample_total_len$segment_len[hit]
  names(sample_class_prop)[1] <- "name"

  # Chromosome-wise class ratio summary
  chr_class_prop <- aggregate(x = segment_len ~ chr + class, FUN = sum,
                              drop = FALSE, data = target_data)
  chr_class_prop$segment_len[is.na(chr_class_prop$segment_len)] <- 0
  chr_total_len <- aggregate(x = segment_len ~ chr, FUN = sum,
                             drop = FALSE, data = chr_class_prop)
  hit <- match(chr_class_prop$chr, chr_total_len$chr)
  chr_class_prop$value <- chr_class_prop$segment_len / chr_total_len$segment_len[hit]
  names(chr_class_prop)[1] <- "name"

  # Whole-genome class ratio summary
  genome_class_prop <- aggregate(x = segment_len ~ class, FUN = sum,
                                 drop = FALSE, data = target_data)
  genome_class_prop$segment_len[is.na(genome_class_prop$segment_len)] <- 0
  genome_class_prop$value <- genome_class_prop$segment_len / sum(genome_class_prop$segment_len)
  genome_class_prop$name <- "genome"

  if(data == "haplotype"){
    # Haploid-wise recombination frequency summary
    hap_recomb <- aggregate(x = class ~ name + chr, FUN = .getRecombNumber,
                            drop = FALSE, data = target_data)
    hap_recomb <- aggregate(x = class ~ name, FUN = sum,
                            drop = FALSE, data = hap_recomb)
    names(hap_recomb) <- c("name", "value")

    # Sample-wise recombination frequency summary
    sample_recomb <- aggregate(x = class ~ name + chr + sample_id, FUN = .getRecombNumber,
                               drop = TRUE, data = target_data)
    sample_recomb <- aggregate(x = class ~ sample_id, FUN = sum,
                               drop = FALSE, data = sample_recomb)
    names(sample_recomb) <- c("name", "value")

    # Chromosome-wise recombination frequency summary
    chr_recomb <- aggregate(x = class ~ name + chr + sample_id, FUN = .getRecombNumber,
                               drop = TRUE, data = target_data)
    chr_recomb <- aggregate(x = class ~ sample_id + chr, FUN = sum,
                            drop = FALSE, data = target_data)
    chr_recomb <- aggregate(x = class ~ chr, FUN = sum,
                            drop = FALSE, data = chr_recomb)
    names(chr_recomb) <- c("name", "value")
  } else {
    hap_recomb <- sample_recomb <- chr_recomb <- NULL
  }

  out <- list(hap_class = hap_class_prop, sample_class = sample_class_prop,
              chr_class = chr_class_prop, genome_class = genome_class_prop,
              hap_recomb = hap_recomb, sample_recomb = sample_recomb,
              chr_recomb = chr_recomb)
  return(out)
}

.getRecombNumber <- function(x){
  return(sum(diff(na.omit(x)) != 0))
}

.evalBlockLen <- function(object, data, chr_len){
  target_data <- object$segments[[data]]
  target_data <- subset(target_data, subset = !is.na(class))

  hap_segment <- rbind(data.frame(stats = "mean",
                                  aggregate(x = segment_len ~ name,
                                            FUN = mean,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "sd",
                                  aggregate(x = segment_len ~ name,
                                            FUN = sd,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "var",
                                  aggregate(x = segment_len ~ name,
                                            FUN = var,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "min",
                                  aggregate(x = segment_len ~ name,
                                            FUN = min,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "Q1",
                                  aggregate(x = segment_len ~ name,
                                            FUN = quantile,
                                            prob = 0.25,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "Q2",
                                  aggregate(x = segment_len ~ name,
                                            FUN = quantile,
                                            prob = 0.5,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "Q3",
                                  aggregate(x = segment_len ~ name,
                                            FUN = quantile,
                                            prob = 0.75,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "max",
                                  aggregate(x = segment_len ~ name,
                                            FUN = max,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)))
  names(hap_segment)[3] <- "value"

  sample_segment <- rbind(data.frame(stats = "mean",
                                     aggregate(x = segment_len ~ sample_id,
                                               FUN = mean,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "sd",
                                     aggregate(x = segment_len ~ sample_id,
                                               FUN = sd,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "var",
                                     aggregate(x = segment_len ~ sample_id,
                                               FUN = var,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "min",
                                     aggregate(x = segment_len ~ sample_id,
                                               FUN = min,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "Q1",
                                     aggregate(x = segment_len ~ sample_id,
                                               FUN = quantile,
                                               prob = 0.25,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "Q2",
                                     aggregate(x = segment_len ~ sample_id,
                                               FUN = quantile,
                                               prob = 0.5,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "Q3",
                                     aggregate(x = segment_len ~ sample_id,
                                               FUN = quantile,
                                               prob = 0.75,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "max",
                                     aggregate(x = segment_len ~ sample_id,
                                               FUN = max,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)))
  names(sample_segment)[2:3] <- c("name", "value")

  chr_segment <- rbind(data.frame(stats = "mean",
                                  aggregate(x = segment_len ~ chr,
                                            FUN = mean,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "sd",
                                  aggregate(x = segment_len ~ chr,
                                            FUN = sd,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "var",
                                  aggregate(x = segment_len ~ chr,
                                            FUN = var,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "min",
                                  aggregate(x = segment_len ~ chr,
                                            FUN = min,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "Q1",
                                  aggregate(x = segment_len ~ chr,
                                            FUN = quantile,
                                            prob = 0.25,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "Q2",
                                  aggregate(x = segment_len ~ chr,
                                            FUN = quantile,
                                            prob = 0.5,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "Q3",
                                  aggregate(x = segment_len ~ chr,
                                            FUN = quantile,
                                            prob = 0.75,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)),
                       data.frame(stats = "max",
                                  aggregate(x = segment_len ~ chr,
                                            FUN = max,
                                            na.rm = TRUE,
                                            drop = FALSE,
                                            data = target_data)))
  names(chr_segment)[2:3] <- c("name", "value")

  target_data$all <- "genome"
  genome_segment <- rbind(data.frame(stats = "mean",
                                     aggregate(x = segment_len ~ all,
                                               FUN = mean,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "sd",
                                     aggregate(x = segment_len ~ all,
                                               FUN = sd,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "var",
                                     aggregate(x = segment_len ~ all,
                                               FUN = var,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "min",
                                     aggregate(x = segment_len ~ all,
                                               FUN = min,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "Q1",
                                     aggregate(x = segment_len ~ all,
                                               FUN = quantile,
                                               prob = 0.25,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "Q2",
                                     aggregate(x = segment_len ~ all,
                                               FUN = quantile,
                                               prob = 0.5,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "Q3",
                                     aggregate(x = segment_len ~ all,
                                               FUN = quantile,
                                               prob = 0.75,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)),
                          data.frame(stats = "max",
                                     aggregate(x = segment_len ~ all,
                                               FUN = max,
                                               na.rm = TRUE,
                                               drop = FALSE,
                                               data = target_data)))
  names(genome_segment)[2:3] <- c("name", "value")
  out <- list(hap_segment = hap_segment, sample_segment = sample_segment,
              chr_segment = chr_segment, genome_segment = genome_segment)
}
