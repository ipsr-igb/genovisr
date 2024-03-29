#' Create a genovis class object storing dummy data
#'
#' @param n_sample An integer to specify the number of samples.
#' @param n_marker An integer to specify the number of markers.
#' @param n_ploidy An integer to specify the ploidy.
#' @param n_haplotype An integer to specify the number of haplotypes.
#'
#' @examples
#' # Create simulation data with the default settings.
#' genovis <- simGenovis()
#'
#' # Create simulation data with tweaking settings.
#' genovis <- simGenovis(n_sample = 100, n_amrker = 2500, n_chr = 5)
#'
#' @export
#'

simGenovis <- function(n_sample = 10,
                       n_marker = 1000,
                       n_chr = 12,
                       chr_len = sample(20:40, n_chr, replace = TRUE),
                       n_ploidy = 2,
                       n_haplotype = 2){
  # Simulate genotypes of each haplotype
  geno4hap <- .simGeno4Hap(n_haplotype = n_haplotype,
                           n_marker = n_marker,
                           n_chr = n_chr)

  # Simulate halpotypes of samples
  sample_hap <- .simHap(n_haplotype = n_haplotype,
                        n_marker = n_marker,
                        n_chr = n_chr,
                        n_sample = n_sample,
                        n_ploidy = n_ploidy,
                        chr_len = chr_len)

  attributes(sample_hap)$scale_breaks <- seq_len(n_haplotype)
  attributes(sample_hap)$scale_labels <- paste0("Hap", seq_len(n_haplotype))

  # Get genotypes of samples based on the simulated haplotypes
  sample_geno <- apply(X = sample_hap,
                       MARGIN = 2,
                       FUN = .mapGeno,
                       geno4hap = geno4hap,
                       n_ploidy = n_ploidy)

  # Calculate dosages at markers in samples if the number of haplotypes is 2.
  if(n_haplotype == 2){
    sample_ds <- apply(X = sample_hap - 1,
                       MARGIN = c(2, 3),
                       FUN = sum)

  } else {
    sample_ds <- NULL
  }

  # Generate marker information.
  chr <- rep(x = seq_len(n_chr), each = n_marker)
  pos <- sapply(seq_len(n_chr), function(i){
    return(sort(sample(x = seq_len(chr_len[i] * 1e6), size = n_marker)))
  })
  pos <- as.vector(pos)

  # Organize an output as a genovis class object.
  out <- list(marker_info = data.frame(id = paste(chr, pos, sep = "_"),
                                       chr = chr,
                                       pos= pos,
                                       ref_allele = "G",
                                       alt_allele = "A"),
              sample_info = data.frame(id = paste("sample",
                                                  seq_len(n_sample),
                                                  sep = "_")),
              genotype = t(sample_geno),
              haplotype = sample_hap,
              dosage = sample_ds)

  # Set the genovis class to the output.
  class(out) <- c(class(out), "genovis")

  # return the genovis class object.
  return(out)
}

.simGeno4Hap <- function(n_haplotype, n_marker, n_chr){
  out <- matrix(data = NA, nrow = n_haplotype, ncol = n_marker * n_chr)
  for(i in seq_len(n_haplotype)){
    out[i, ] <- sample(x = c(0, 1), size = n_marker * n_chr, replace = TRUE)
  }
  return(out)
}

.simHap <- function(n_haplotype, n_marker, n_chr, n_sample, n_ploidy, chr_len){
  n_xo <- rpois(n = n_sample * n_ploidy, lambda = sum(chr_len) * 0.04)
  xo_pos <- lapply(X = seq_along(n_xo),
                   FUN = .getXoPos,
                   n_xo = n_xo,
                   n_marker = n_marker,
                   n_chr = n_chr)
  xo_pos <- do.call("c", xo_pos)
  hap_block <- cbind(c(0, xo_pos),
                     c(xo_pos, n_marker * n_chr * n_sample * n_ploidy),
                     sample(x = seq_len(n_haplotype),
                            size = length(xo_pos) + 1,
                            replace = TRUE))
  out <- apply(X = hap_block,
               MARGIN = 1,
               FUN = function(x)rep(x[3], x[2] - x[1]))
  out <- unlist(out)
  dim(out) <- c(n_sample, n_marker * n_chr, n_ploidy)
  out <- apply(X = out, MARGIN = 2, FUN = t)
  dim(out) <- c(n_sample, n_ploidy, n_marker * n_chr)
  out <- apply(X = out, MARGIN = 3, FUN = t)
  dim(out) <- c(n_ploidy, n_sample, n_marker * n_chr)
  return(out)
}

.getXoPos <- function(i, n_xo, n_marker, n_chr, n_ploidy){
  out <- sort(sample(x = seq_len(n_marker * n_chr), size = n_xo[i]))
  out <- out + n_marker * n_chr * (i - 1)
  return(out)
}

.mapGeno <- function(sample_hap, geno4hap, n_ploidy){
  sample_hap <- t(sample_hap)
  sample_hap  <- sample_hap + 0:(nrow(sample_hap) - 1) * n_ploidy
  out <- geno4hap[as.vector(sample_hap)]
  dim(out) <- c(ncol(geno4hap), n_ploidy)
  out <- rowSums(out)
  return(out)
}

