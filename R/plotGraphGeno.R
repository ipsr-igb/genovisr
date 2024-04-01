#'
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#'

plotGraphGeno <- function(object,
                          data = "genotype",
                          sample = NULL,
                          marker = NULL,
                          direction = "h"){
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
  
  if(direction == "v" & sum(sample) > 1){
    message("When 'v' was specified to 'direction', ",
            "a graphical genotype plot for only one sample can be drawn.")
    message("The first sample was chosen for plotting.")
    sample[-which(sample)[1]] <- FALSE
  }
  
  if(is.null(marker)){
    marker <- rep(TRUE, nrow(object$marker_info))
    
  } else {
    if(is.character(marker)){
      marker <- object$marker_info$id %in% marker
      
    } else if(is.character(marker)){
      marker <- object$marker_info$id %in% marker
      
    } else if(is.numeric(marker)){
      marker <- seq_along(object$marker_info$id) %in% marker
    }
  }
  
  data <- match.arg(arg = data, choices = c("genotype", "haplotype", "dosage"))
  if(data == "genotype"){
    dat <- object$genotype[sample, marker]
    rownames(dat) <- object$sample_info$id[sample]
    legend <- "Genotype"
    sample_lables <- rownames(dat)
    scale_breaks <- 0:2
    scale_labels <- c("Ref", "Het", "Alt")
    
  } else if(data == "haplotype"){
    dat <- object$haplotype[, sample, marker]
    n_hap <- dim(dat)[1]
    n_sample <- dim(dat)[2]
    dat <- apply(X = dat, MARGIN = 3, c)
    rownames(dat) <- paste(rep(object$sample_info$id[sample], each = n_hap),
                           paste0("hap", seq_len(n_hap)), sep = "_")
    legend <- "Haplotype"
    sample_lables <- rownames(dat)
    scale_breaks <- attributes(object$haplotype)$scale_breaks
    scale_labels <- attributes(object$haplotype)$scale_labels
    
  } else if(data == "dosage"){
    dat <- object$dosage[sample, marker]
    dosage_levels <- sort(unique(as.vector(dat)))
    rownames(dat) <- object$sample_info$id[sample]
    legend <- "Dosage"
    sample_lables <- rownames(dat)
    scale_breaks <- dosage_levels
    scale_labels <- paste0("Plex", dosage_levels)
  }
  
  df <- cbind(subset(object$marker_info[marker, ], select = chr:pos), t(dat))
  df <- .getRanges(df = df)
  df <- pivot_longer(data = df, cols = -(chr:pos))
  df$value <- factor(df$value, levels = scale_breaks)
  df$pos <- factor(df$pos)
  
  if(direction == "v"){
    df$name <- factor(df$name, levels = sample_lables)
    p <- ggplot(df, aes(x = name, y = pos, fill = value)) +
      facet_wrap(facets = ~ chr, nrow = 1, scales = "free_y")
    
  } else {
    df$name <- factor(df$name, levels = rev(sample_lables))
    p <- ggplot(df, aes(x = pos, y = name, fill = value)) +
      facet_wrap(facets = ~ chr, nrow = 1, scales = "free_x")
  }
  
  p <- p + geom_tile() +
    scale_fill_viridis_d(name = legend,
                         breaks = scale_breaks,
                         labels = scale_labels) +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          legend.position = "top",
          panel.spacing = unit(0, "lines"))
}

.getRanges <- function(df){
  xo <- apply(X = subset(df, select = -(chr:pos)), MARGIN = 2, diff)
  xo_i <- apply(X = xo, MARGIN = 2, function(x){x != 0})
  chr_boundary <- diff(as.numeric(factor(df$chr))) != 0
  blocks <- xo_i | chr_boundary
  out <- lapply(X = seq_len(ncol(blocks)), function(i){
    blocks_i <- which(blocks[, i])
    blocks_df <- data.frame(start = c(0, ))
  })
}
