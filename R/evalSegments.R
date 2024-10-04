#' Evaluate Segments
#'
#' This function evaluates the segments for haplotype and dosage data in a given `genovis` object.
#'
#' @param object A `genovis` class object.
#' @param marker A character or numeric vector specifying markers to include in the evaluation. Default is `NULL`.
#' @param data A character vector specifying which data to use, either "haplotype" or "dosage". Default is c("haplotype", "dosage").
#' @return The input `genovis` object with evaluated segments.
#' @export
evalSegments <- function(object, marker = NULL, data = c("haplotype", "dosage")) {

  # Check if the input object is of class 'genovis'
  if (!inherits(x = object, what = "genovis")) {
    stop("The input object should be a genovis class object.")
  }

  # If no markers are specified, select all markers
  if (is.null(marker)) {
    marker <- rep(TRUE, nrow(object$marker_info))
  } else {
    # Determine markers based on the type of input
    if (is.character(marker)) {
      marker <- object$marker_info$id %in% marker
    } else if (is.numeric(marker)) {
      marker <- seq_along(object$marker_info$id) %in% marker
    }
  }

  # Match the data argument with allowed choices
  data <- match.arg(arg = data, choices = c("haplotype", "dosage"), several.ok = TRUE)

  options(scipen = 10)  # Set scientific notation penalty

  object$segments <- NULL

  # Process haplotype data if specified
  if ("haplotype" %in% data) {
    if (is.null(object$haplotype)) {
      message("Haplotype information is not available for the given dataset.")
    } else {
      out <- .getSegmetsHaplotype(object = object, marker = marker)
      object$segments <- c(object$segments, list(haplotype = out))
    }
  }

  # Process dosage data if specified
  if ("dosage" %in% data) {
    if (is.null(object$dosage)) {
      message("Dosage information is not available for the given dataset.")
    } else {
      out <- .getSegmetsDosage(object = object, marker = marker)
      object$segments <- c(object$segments, list(dosage = out))
    }
  }

  return(object)
}

#' Get Segments for Haplotype Data
#'
#' @param object A `genovis` class object.
#' @param marker A logical vector indicating selected markers.
#' @return A dataframe containing the segments for haplotype data.
.getSegmetsHaplotype <- function(object, marker) {
  dat <- object$haplotype[,, marker]
  n_hap <- dim(dat)[1]
  dat <- apply(X = dat, MARGIN = 3, c)
  rownames(dat) <- paste(rep(object$sample_info$id, each = n_hap),
                         paste0("hap", seq_len(n_hap)), sep = "_")
  df <- cbind(subset(object$marker_info[marker, ], select = chr:pos), t(dat))
  out <- .getSegmetsEngine(df = df)
  return(out)
}

#' Get Segments for Dosage Data
#'
#' @param object A `genovis` class object.
#' @param marker A logical vector indicating selected markers.
#' @return A dataframe containing the segments for dosage data.
.getSegmetsDosage <- function(object, marker) {
  dat <- object$dosage[, marker]
  rownames(dat) <- object$sample_info$id
  df <- cbind(subset(object$marker_info[marker, ], select = chr:pos), t(dat))
  out <- .getSegmetsEngine(df = df)
  return(out)
}

#' Get Segments Engine
#'
#' @param df A dataframe containing marker and data information.
#' @return A dataframe with evaluated segments.
.getSegmetsEngine <- function(df) {
  df[is.na(df)] <- 100
  xo <- apply(X = subset(df, select = -(chr:pos)), MARGIN = 2, diff)
  xo_i <- apply(X = xo, MARGIN = 2, function(x) {x != 0})
  chr_boundary <- diff(as.numeric(factor(df$chr))) != 0
  blocks <- xo_i | chr_boundary
  out <- lapply(X = seq_len(ncol(blocks)), function(i) {
    blocks_i <- which(blocks[, i])
    blocks_s <- c(0, blocks_i) + 1
    blocks_e <- c(blocks_i, nrow(df))
    name <- colnames(df)[i + 2]
    sample_id <- sub("_hap[0-9]", "", name)
    hap_id <- sub(paste0(sample_id, "_"), "", name)
    blocks_df <- data.frame(name = colnames(df)[i + 2],
                            sample_id = sample_id,
                            hap_id = hap_id,
                            start_index = blocks_s,
                            end_index = blocks_e,
                            chr = df$chr[blocks_e],
                            start_pos = df$pos[blocks_s],
                            end_pos = df$pos[blocks_e],
                            class = df[blocks_s, i + 2])
  })
  out <- do.call("rbind", out)
  out$class[out$class == 100] <- NA
  return(out)
}
