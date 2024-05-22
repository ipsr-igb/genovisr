#' 
#' 
#' 
#' 
#' @export
#' 
statsGeno <- function(object, chr_len){
  if(!inherits(x = object, what = "genovis")){
    stop("The input object should be a genovis class object.")
  }
  
  if(is.null(object$segments)){
    stop('Run evalSegments() to prepare data to plot graphical genotypes.')
  }
  
  object <- .extendBorders(object = object, chr_len = chr_len)
  hap_stats <- .evalStats(object = object, data = "haplotype", chr_len = chr_len)
  dos_stats <- .evalStats(object = object, data = "dosage", chr_len = chr_len)
  object$stats <- list(haplotype = hap_stats, dosage = dos_stats)
  return(object)
}

.extendBorders <- function(object, chr_len){
  for(data in c("haplotype", "dosage")){
    chr_first <- which(!duplicated(object$marker_info$chr))
    to_1 <- object$segments[[data]]$start_index %in% chr_first
    object$segments[[data]]$start_border <- NA
    object$segments[[data]]$start_border[to_1] <- 1
    
    chr_end <- c(chr_first[-1] - 1, nrow(object$marker_info))
    to_end <- object$segments[[data]]$end_index %in% chr_end
    hit <- match(object$segments[[data]]$chr[to_end], seq_along(chr_len))
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

.evalStats <- function(object, data, chr_len){
  classes <- attributes(object[[data]])$scale_breaks
  
  hap_class_prop <- aggregate(x = block_len ~ name + value, FUN = sum, 
                              drop = FALSE, data = object$segments[[data]])
  hap_class_prop$prop_classes <- hap_class_prop$block_len / sum(chr_len)
  hap_class_prop$prop_classes[is.na(hap_class_prop$prop_classes)] <- 0
  
  object$segments[[data]]$sample <- sub("_hap.*", "", object$segments[[data]]$name)
  sample_class_prop <- aggregate(x = block_len ~ sample + value, FUN = sum, 
                                 drop = FALSE, data = object$segments[[data]])
  sample_class_prop$prop_classes <- sample_class_prop$block_len / sum(chr_len * 2)
  sample_class_prop$prop_classes[is.na(sample_class_prop$prop_classes)] <- 0
  
  n_hap <- nrow(hap_class_prop) / 2
  chr_class_prop <- aggregate(x = block_len ~ chr + value, FUN = sum, 
                              drop = FALSE, data = object$segments[[data]])
  chr_class_prop$prop_classes <- chr_class_prop$block_len / (chr_len * n_hap)
  chr_class_prop$prop_classes[is.na(chr_class_prop$prop_classes)] <- 0
  
  genome_class_prop <- aggregate(x = block_len ~ value, FUN = sum, 
                                 drop = FALSE, data = object$segments[[data]])
  genome_class_prop$prop_classes <- genome_class_prop$block_len / sum(chr_len * n_hap)
  genome_class_prop$prop_classes[is.na(genome_class_prop$prop_classes)] <- 0
  
  hap_block <- rbind(data.frame(stats = "mean", 
                                aggregate(x = block_len ~ name, 
                                          FUN = mean, 
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "sd", 
                                aggregate(x = block_len ~ name, 
                                          FUN = sd,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "var", 
                                aggregate(x = block_len ~ name, 
                                          FUN = var, 
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "min", 
                                aggregate(x = block_len ~ name, 
                                          FUN = min,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q1", 
                                aggregate(x = block_len ~ name, 
                                          FUN = quantile,
                                          prob = 0.25,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q2", 
                                aggregate(x = block_len ~ name, 
                                          FUN = quantile,
                                          prob = 0.5,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q3", 
                                aggregate(x = block_len ~ name, 
                                          FUN = quantile,
                                          prob = 0.75,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "max", 
                                aggregate(x = block_len ~ name, 
                                          FUN = max,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])))
  
  sample_block <- rbind(data.frame(stats = "mean", 
                                   aggregate(x = block_len ~ sample, 
                                             FUN = mean, 
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "sd", 
                                   aggregate(x = block_len ~ sample, 
                                             FUN = sd,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "var", 
                                   aggregate(x = block_len ~ sample, 
                                             FUN = var, 
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "min", 
                                   aggregate(x = block_len ~ sample, 
                                             FUN = min,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q1", 
                                   aggregate(x = block_len ~ sample, 
                                             FUN = quantile,
                                             prob = 0.25,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q2", 
                                   aggregate(x = block_len ~ sample, 
                                             FUN = quantile,
                                             prob = 0.5,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q3", 
                                   aggregate(x = block_len ~ sample, 
                                             FUN = quantile,
                                             prob = 0.75,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "max", 
                                   aggregate(x = block_len ~ sample, 
                                             FUN = max,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])))
  names(sample_block)[3] <- "value"
  
  chr_block <- rbind(data.frame(stats = "mean", 
                                aggregate(x = block_len ~ chr, 
                                          FUN = mean, 
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "sd", 
                                aggregate(x = block_len ~ chr, 
                                          FUN = sd,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "var", 
                                aggregate(x = block_len ~ chr, 
                                          FUN = var, 
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "min", 
                                aggregate(x = block_len ~ chr, 
                                          FUN = min,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q1", 
                                aggregate(x = block_len ~ chr, 
                                          FUN = quantile,
                                          prob = 0.25,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q2", 
                                aggregate(x = block_len ~ chr, 
                                          FUN = quantile,
                                          prob = 0.5,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "Q3", 
                                aggregate(x = block_len ~ chr, 
                                          FUN = quantile,
                                          prob = 0.75,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])),
                     data.frame(stats = "max", 
                                aggregate(x = block_len ~ chr, 
                                          FUN = max,
                                          na.rm = TRUE,
                                          drop = FALSE,
                                          data = object$segments[[data]])))
  names(chr_block)[3] <- "value"
  
  object$segments[[data]]$all <- 1
  genome_block <- rbind(data.frame(stats = "mean", 
                                   aggregate(x = block_len ~ all, 
                                             FUN = mean, 
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "sd", 
                                   aggregate(x = block_len ~ all, 
                                             FUN = sd,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "var", 
                                   aggregate(x = block_len ~ all, 
                                             FUN = var, 
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "min", 
                                   aggregate(x = block_len ~ all, 
                                             FUN = min,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q1", 
                                   aggregate(x = block_len ~ all, 
                                             FUN = quantile,
                                             prob = 0.25,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q2", 
                                   aggregate(x = block_len ~ all, 
                                             FUN = quantile,
                                             prob = 0.5,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "Q3", 
                                   aggregate(x = block_len ~ all, 
                                             FUN = quantile,
                                             prob = 0.75,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])),
                        data.frame(stats = "max", 
                                   aggregate(x = block_len ~ all, 
                                             FUN = max,
                                             na.rm = TRUE,
                                             drop = FALSE,
                                             data = object$segments[[data]])))
  
  out <- list(hap_class = hap_class_prop, sample_class = sample_class_prop, 
              chr_class = chr_class_prop, genome_class = genome_class_prop,
              hap_block = hap_block, sample_block = sample_block,
              chr_block = chr_block, genome_block = genome_block)
  return(out)
}
