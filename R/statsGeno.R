#'
#'
#'
#'
#' @export
#'
statsGeno <- function(object, chr_len = NULL){
  if(!inherits(x = object, what = "genovis")){
    stop("The input object should be a genovis class object.")
  }

  if(is.null(object$segments)){
    stop('Run evalSegments() to prepare data to plot graphical genotypes.')
  }

  if(is.null(chr_len)){
    chr_len <- tapply(X = object$marker_info$pos,
                      INDEX = object$marker_info$chr,
                      FUN = max)
  }
  if(is.null(names(chr_len))){
    names(chr_len) <- seq_along(chr_len)
  }

  object <- .extendBorders(object = object, chr_len = chr_len)

  if(is.null(object$haplotype)){
    message("Haplotype information is not available for the given dataset.")
    hap_stats <- NULL

  } else {
    hap_stats <- .evalStats(object = object, data = "haplotype", chr_len = chr_len)
  }

  if(is.null(object$dosage)){
    message("Dosage information is not available for the given dataset.")
    dos_stats <- NULL

  } else {
    dos_stats <- .evalStats(object = object, data = "dosage", chr_len = chr_len)
  }

  object$stats <- list(haplotype = hap_stats, dosage = dos_stats)
  return(object)
}

.extendBorders <- function(object, chr_len){
  for(data in c("haplotype", "dosage")){
    if(is.null(object$segments[[data]])){
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
    if(check){
      stop('Some segment(s) have lengths of less than zero.',
           "\nPlease confirm the values specified to 'chr_len'.")
    }
  }
  return(object)
}

.evalStats <- function(object, data, chr_len){
  object$segments[[data]]$sample <- sub("_hap.*", "", object$segments[[data]]$name)
  out1 <- .evalClassProp(object = object, data = data, chr_len = chr_len)
  out2 <- .evalBlockLen(object = object, data = data, chr_len = chr_len)

  out <- c(out1, out2)
  return(out)
}

.evalClassProp <- function(object, data, chr_len){
  classes <- attributes(object[[data]])$scale_breaks

  # Haploid-wise class ratio summary
  hap_class_prop <- aggregate(x = segment_len ~ name + class, FUN = sum,
                              drop = FALSE, data = object$segments[[data]])
  hap_class_prop$segment_len[is.na(hap_class_prop$segment_len)] <- 0
  hap_total_len <- aggregate(x = segment_len ~ name, FUN = sum,
                             drop = FALSE, data = hap_class_prop)
  hit <- match(hap_class_prop$name, hap_total_len$name)
  hap_class_prop$value <- hap_class_prop$segment_len / hap_total_len$segment_len[hit]

  # Sample-wise class ratio summary
  sample_class_prop <- aggregate(x = segment_len ~ sample + class, FUN = sum,
                                 drop = FALSE, data = object$segments[[data]])
  sample_class_prop$segment_len[is.na(sample_class_prop$segment_len)] <- 0
  sample_total_len <- aggregate(x = segment_len ~ sample, FUN = sum,
                                drop = FALSE, data = sample_class_prop)
  hit <- match(sample_class_prop$sample, sample_total_len$sample)
  sample_class_prop$value <- sample_class_prop$segment_len / sample_total_len$segment_len[hit]
  names(sample_class_prop)[1] <- "name"

  # Chromosome-wise class ratio summary
  chr_class_prop <- aggregate(x = segment_len ~ chr + class, FUN = sum,
                              drop = FALSE, data = object$segments[[data]])
  chr_class_prop$segment_len[is.na(chr_class_prop$segment_len)] <- 0
  chr_total_len <- aggregate(x = segment_len ~ chr, FUN = sum,
                             drop = FALSE, data = chr_class_prop)
  hit <- match(chr_class_prop$chr, chr_total_len$chr)
  chr_class_prop$value <- chr_class_prop$segment_len / chr_total_len$segment_len[hit]
  names(chr_class_prop)[1] <- "name"

  # Whole-genome class ratio summary
  genome_class_prop <- aggregate(x = segment_len ~ class, FUN = sum,
                                 drop = FALSE, data = object$segments[[data]])
  genome_class_prop$segment_len[is.na(genome_class_prop$segment_len)] <- 0
  genome_class_prop$value <- genome_class_prop$segment_len / sum(genome_class_prop$segment_len)
  genome_class_prop$name <- "genome"

  if(data == "haplotype"){
    # Haploid-wise recombination frequency summary
    hap_recomb <- aggregate(x = segment_len ~ name + chr, FUN = length,
                            drop = FALSE, data = object$segments[[data]])
    hap_recomb$segment_len <- hap_recomb$segment_len - 1
    hap_recomb <- aggregate(x = segment_len ~ name, FUN = sum,
                            drop = FALSE, data = hap_recomb)
    names(hap_recomb) <- c("name", "value")

    # Sample-wise recombination frequency summary
    sample_recomb <- aggregate(x = segment_len ~ name + chr + sample, FUN = length,
                               drop = FALSE, data = object$segments[[data]])
    sample_recomb$segment_len <- sample_recomb$segment_len - 1
    sample_recomb <- aggregate(x = segment_len ~ sample, FUN = sum,
                               drop = FALSE, data = sample_recomb)
    names(sample_recomb) <- c("name", "value")

    # Chromosome-wise recombination frequency summary
    chr_recomb <- aggregate(x = segment_len ~ sample + chr, FUN = length,
                            drop = FALSE, data = object$segments[[data]])
    chr_recomb$segment_len <- chr_recomb$segment_len - 1
    chr_recomb <- aggregate(x = segment_len ~ chr, FUN = sum,
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

.evalBlockLen <- function(object, data, chr_len){
  hap_segment <- rbind(data.frame(stats = "mean",
                                aggregate(x = segment_len ~ name,
                                          FUN = mean,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "sd",
                                aggregate(x = segment_len ~ name,
                                          FUN = sd,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "var",
                                aggregate(x = segment_len ~ name,
                                          FUN = var,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "min",
                                aggregate(x = segment_len ~ name,
                                          FUN = min,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q1",
                                aggregate(x = segment_len ~ name,
                                          FUN = quantile,
                                          prob = 0.25,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q2",
                                aggregate(x = segment_len ~ name,
                                          FUN = quantile,
                                          prob = 0.5,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q3",
                                aggregate(x = segment_len ~ name,
                                          FUN = quantile,
                                          prob = 0.75,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "max",
                                aggregate(x = segment_len ~ name,
                                          FUN = max,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])))
  names(hap_segment)[3] <- "value"

  sample_segment <- rbind(data.frame(stats = "mean",
                                   aggregate(x = segment_len ~ sample,
                                             FUN = mean,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "sd",
                                   aggregate(x = segment_len ~ sample,
                                             FUN = sd,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "var",
                                   aggregate(x = segment_len ~ sample,
                                             FUN = var,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "min",
                                   aggregate(x = segment_len ~ sample,
                                             FUN = min,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q1",
                                   aggregate(x = segment_len ~ sample,
                                             FUN = quantile,
                                             prob = 0.25,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q2",
                                   aggregate(x = segment_len ~ sample,
                                             FUN = quantile,
                                             prob = 0.5,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q3",
                                   aggregate(x = segment_len ~ sample,
                                             FUN = quantile,
                                             prob = 0.75,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "max",
                                   aggregate(x = segment_len ~ sample,
                                             FUN = max,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])))
  names(sample_segment)[2:3] <- c("name", "value")

  chr_segment <- rbind(data.frame(stats = "mean",
                                aggregate(x = segment_len ~ chr,
                                          FUN = mean,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "sd",
                                aggregate(x = segment_len ~ chr,
                                          FUN = sd,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "var",
                                aggregate(x = segment_len ~ chr,
                                          FUN = var,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "min",
                                aggregate(x = segment_len ~ chr,
                                          FUN = min,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q1",
                                aggregate(x = segment_len ~ chr,
                                          FUN = quantile,
                                          prob = 0.25,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q2",
                                aggregate(x = segment_len ~ chr,
                                          FUN = quantile,
                                          prob = 0.5,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q3",
                                aggregate(x = segment_len ~ chr,
                                          FUN = quantile,
                                          prob = 0.75,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "max",
                                aggregate(x = segment_len ~ chr,
                                          FUN = max,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])))
  names(chr_segment)[2:3] <- c("name", "value")

  object$segments[[data]]$all <- "genome"
  genome_segment <- rbind(data.frame(stats = "mean",
                                   aggregate(x = segment_len ~ all,
                                             FUN = mean,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "sd",
                                   aggregate(x = segment_len ~ all,
                                             FUN = sd,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "var",
                                   aggregate(x = segment_len ~ all,
                                             FUN = var,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "min",
                                   aggregate(x = segment_len ~ all,
                                             FUN = min,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q1",
                                   aggregate(x = segment_len ~ all,
                                             FUN = quantile,
                                             prob = 0.25,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q2",
                                   aggregate(x = segment_len ~ all,
                                             FUN = quantile,
                                             prob = 0.5,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q3",
                                   aggregate(x = segment_len ~ all,
                                             FUN = quantile,
                                             prob = 0.75,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "max",
                                   aggregate(x = segment_len ~ all,
                                             FUN = max,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])))
  names(genome_segment)[2:3] <- c("name", "value")
  out <- list(hap_segment = hap_segment, sample_segment = sample_segment,
              chr_segment = chr_segment, genome_segment = genome_segment)
}
