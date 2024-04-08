#'
#'
#'
#'
#'
#'

getRanges <- function(object, marker = NULL, 
                      data = c("haplotype", "dosage")){
  
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
  
  data <- match.arg(arg = data,
                    choices = c("haplotype", "dosage"),
                    several.ok = TRUE)
  
  object$ranges <- NULL
  if("haplotype" %in% data){
    out <- .getRangesHaplotype(object = object, marker = marker)
    object$ranges <- c(object$ranges, list(haplotype = out))
  }
  if("dosage" %in% data){
    out <- .getRangesDosage(object = object, marker = marker)
    object$ranges <- c(object$ranges, list(dosage = out))
  }
  return(object)
}

.getRangesHaplotype <- function(object, marker){
  dat <- object$haplotype[,, marker]
  n_hap <- dim(dat)[1]
  dat <- apply(X = dat, MARGIN = 3, c)
  rownames(dat) <- paste(rep(object$sample_info$id, each = n_hap),
                         paste0("hap", seq_len(n_hap)), sep = "_")
  df <- cbind(subset(object$marker_info[marker, ], select = chr:pos), t(dat))
  out <- .getRangesEngine(df = df)
  return(out)
}

.getRangesDosage <- function(object, marker){
  dat <- object$dosage[, marker]
  rownames(dat) <- object$sample_info$id
  df <- cbind(subset(object$marker_info[marker, ], select = chr:pos), t(dat))
  out <- .getRangesEngine(df = df)
  return(out)
}

.getRangesEngine <- function(df){
  xo <- apply(X = subset(df, select = -(chr:pos)), MARGIN = 2, diff)
  xo_i <- apply(X = xo, MARGIN = 2, function(x){x != 0})
  chr_boundary <- diff(as.numeric(factor(df$chr))) != 0
  blocks <- xo_i | chr_boundary
  out <- lapply(X = seq_len(ncol(blocks)), function(i){
    blocks_i <- which(blocks[, i])
    blocks_s <- c(0, blocks_i) + 1
    blocks_e <- c(blocks_i, nrow(df))
    blocks_df <- data.frame(name = colnames(df)[i + 2],
                            start_index = blocks_s, 
                            end_index = blocks_e,
                            chr = df$chr[blocks_e],
                            start_pos = df$pos[blocks_s],
                            end_pos = df$pos[blocks_e],
                            value = df[c(blocks_i, nrow(df)), i + 2])
  })
  out <- do.call("rbind", out)
  return(out)
}
