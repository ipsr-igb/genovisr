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
    hap_stats <- .statsHaplotype(object = object)
  }

  if(is.null(object$dosage)){
    message("Dosage information is not available for the given dataset.")
    dos_stats <- NULL

  } else {
    dos_stats <- .statsDosage(object = object)
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
    object$segments[[data]]$block_len <- object$segments[[data]]$end_border - object$segments[[data]]$start_border + 1
  }
  return(object)
}

.statsHaplotype <- function(object){
  classes <- attributes(object$haplotype)$scale_breaks
  sample_id <- unique(object$segments$haplotype$name)

  hapwise_class_prop <- vapply(X = sample_id,
                               FUN = .propClasses,
                               df = object$segments$haplotype,
                               classes = classes,
                               chr_len = chr_len,
                               FUN.VALUE = list(1))
  hapwise_class_prop <- do.call(what = "rbind", args = hapwise_class_prop)
  rownames(hapwise_class_prop) <- NULL
  out <- list(hapwise_class = hapwise_class_prop)

  object$segments$haplotype$name <- sub(pattern = "_hap.+",
                                      replacement = "",
                                      x = object$segments$haplotype$name)
  sample_id <- sub(pattern = "_hap.+", replacement = "", x = sample_id)
  sample_id <- unique(sample_id)
  samplewise_class_prop <- vapply(X = sample_id,
                                  FUN = .propClasses,
                                  df = object$segments$haplotype,
                                  classes = classes,
                                  chr_len = chr_len * 2,
                                  FUN.VALUE = list(1))
  samplewise_class_prop <- do.call(what = "rbind", args = samplewise_class_prop)
  rownames(samplewise_class_prop) <- NULL
  out <- c(out, list(samplewise_class = samplewise_class_prop))

  chr <- sort(unique(object$segments$haplotype$chr))
  chrwise_block_len <- vapply(X = chr,
                              FUN = .blockLen,
                              df = object$segments$haplotype,
                              FUN.VALUE = list(1))
  chrwise_block_len <- do.call(what = "rbind", args = chrwise_block_len)
  rownames(chrwise_block_len) <- NULL
  out <- c(out, list(chrwise_block = chrwise_block_len))

  total_block_len <- .blockLen(chr = chr, df = object$segments$haplotype)
  out <- c(out, list(all_block = total_block_len[[1]]))
  return(out)
}

.statsDosage <- function(object){
  classes <- attributes(object$dosage)$scale_breaks
  sample_id <- unique(object$segments$dosage$name)
  samplewise_class_prop <- vapply(X = sample_id,
                                  FUN = .propClasses,
                                  df = object$segments$dosage,
                                  classes = classes,
                                  chr_len = chr_len,
                                  FUN.VALUE = list(1))
  samplewise_class_prop <- do.call(what = "rbind", args = samplewise_class_prop)
  rownames(samplewise_class_prop) <- NULL
  out <- list(samplewise_class = samplewise_class_prop)

  chr <- sort(unique(object$segments$dosage$chr))
  chrwise_block_len <- vapply(X = chr,
                              FUN = .blockLen,
                              df = object$segments$dosage,
                              FUN.VALUE = list(1))
  chrwise_block_len <- do.call(what = "rbind", args = chrwise_block_len)
  rownames(chrwise_block_len) <- NULL
  out <- c(out, list(chrwise_block = chrwise_block_len))

  total_block_len <- .blockLen(chr = chr, df = object$segments$dosage)
  out <- c(out, list(all_block = total_block_len[[1]]))
  return(out)
}

.propClasses <- function(sample_id, df, classes, chr_len){
  out1 <- aggregate(x = block_len ~ chr + value, FUN = sum, drop = FALSE,
                    data = df, subset = name == sample_id)
  hit <- match(out1$chr, seq_along(chr_len))
  out1$prop_classes <- out1$block_len / chr_len[hit]
  out1$prop_classes[is.na(out1$prop_classes)] <- 0
  out1 <- data.frame(seq_along(chr_len),
                     matrix(out1$prop_classes, nrow = length(chr_len)))
  names(out1) <- c("chr", seq(1, ncol(out1) - 1))
  out2 <- aggregate(x = block_len ~ value, FUN = sum, drop = FALSE,
                    data = df, subset = name == sample_id)
  out2 <- data.frame("total", t(out2$block_len / sum(chr_len)))
  names(out2) <- c("chr", seq(1, ncol(out2) - 1))
  out <- rbind(out1, out2)
  out <- cbind(name = sample_id, out)
  return(list(out))
}

.blockLen <- function(chr, df){
  chr_block_len <- df$block_len[df$chr %in% chr]
  value <- c(mean = mean(chr_block_len),
             sd = sd(chr_block_len),
             var = var(chr_block_len),
             min = min(chr_block_len),
             q1 = quantile(x = chr_block_len, probs = 0.25, names = FALSE),
             q2 = quantile(x = chr_block_len, probs = 0.5, names = FALSE),
             q3 = quantile(x = chr_block_len, probs = 0.75, names = FALSE),
             max = max(chr_block_len))
  if(length(chr) != 1){
    out <- data.frame(chr = "all", stats = names(value), value = value)
  } else {
    out <- data.frame(chr = chr, stats = names(value), value = value)
  }
  return(list(out))
}
