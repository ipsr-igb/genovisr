#' Import markers, genotypes, and haplotypes information from a GDS file
#'
#' Chromosome IDs, physical positions, and alleles of markers stored in an input
#' GDS file are retrieved as well as genotypes and haplotypes information at the
#' markers. The obtained information is then organized as a genovis class
#' object that can be handled by the functions in the genovisr package.
#'
#' @param gbsr A GbsrGenotypeData object. See [GBScleanR::gbsrVCF2GDS()] and
#' [GBScleanR::loadGDS()].
#' @param verbose A logical value to indicate whether to show messages from
#' the function.
#'
#' @importClassesFrom GBScleanR GbsrGenotypeData
#' @importMethodsFrom GBScleanR getGenotype getHaplotype getRead getAllele
#' @importFrom gdsfmt exist.gdsn
#' @importFrom GBScleanR loadGDS
#'
#' @examples
#' # Not run
#' # genovis <- gdscleanr2genovis(gds_fn = "path/to/gds/file.gds",
#' #                              gbscleanr_filter = TRUE)
#'
#' @export
gbscleanr2genovis <- function(gbsr, verbose = TRUE) {
  if (verbose) {
    message('Loading GDS file.')
  }

  # Check if the input is a GbsrGenotypeData object
  if (!inherits(x = gbsr, what = "GbsrGenotypeData")) {
    stop("A path to a GDS file should be specified to 'gds_fn'.")
  }

  parents <- getParents(object = gbsr, verbose = FALSE)

  if(is.null(parents)){
    gbsr <- .setParentsFromGDS(gbsr = gbsr)
    parents <- getParents(object = gbsr, verbose = FALSE)
  }

  # Pull out data required to build a genovis class object.
  # Get marker allele information.
  allele <- getAllele(object = gbsr)
  allele <- strsplit(x = allele, split = ",")
  allele <- do.call("rbind", allele)

  # Confirm whether corrected genotype data is available or not.
  genotype_node <- if (exist.gdsn(node = gbsr$root, path = "annotation/format/CGT")) {
    "cor"
  } else {
    "raw"
  }

  # Organize marker information and genotype information in an output object.
  chr <- getChromosome(object = gbsr)
  pos <- getPosition(object = gbsr)
  out <- list(marker_info = data.frame(id = paste(chr, pos, sep = "_"),
                                       chr = factor(chr),
                                       pos = as.integer(pos),
                                       ref_allele = allele[, 1],
                                       alt_allele = allele[, 2]),
              sample_info = data.frame(id = getSamID(object = gbsr)),
              genotype = getGenotype(object = gbsr, node = genotype_node, parents = TRUE),
              haplotype = NULL,
              dosage = NULL
  )

  # Add haplotype information if available.
  if (exist.gdsn(node = gbsr$root, path = "annotation/format/HAP")) {
    out$haplotype <- .getHaplotypeFromGBSR(gbsr = gbsr)
  } else {
    message("Haplotype information is not available for the given dataset.")
  }

  # Add dosage information if available.
  if (exist.gdsn(node = gbsr$root, path = "annotation/format/EDS")) {
    out$dosage <- .getDosageFromGBSR(gbsr = gbsr)
  } else {
    message("Dosage information is not available for the given dataset.")
  }

  # Set the genovis class to the output.
  class(out) <- c(class(out), "genovis")

  # Return the genovis class object.
  return(out)
}

#' @importFrom gdsfmt index.gdsn read.gdsn exist.gdsn
#' @importFrom GBScleanR setParents
.setParentsFromGDS <- function(gbsr){
  parents_gdsn <- "parents/data"

  if(exist.gdsn(node = gbsr$root, path = parents_gdsn)){
    parents <- read.gdsn(node = index.gdsn(node = gbsr$root,
                                           path = parents_gdsn))
    parents <- getSamID(object = gbsr, valid = FALSE)[parents != 0]
    gbsr <- setParents(object = gbsr, parents = parents)
  }
  return(gbsr)
}

#' Get Haplotype Information from GBSR
#'
#' @param gbsr A GbsrGenotypeData object.
#' @importFrom GBScleanR getParents
#' @return A haplotype matrix with scale breaks and labels as attributes.
.getHaplotypeFromGBSR <- function(gbsr) {
  parents <- getParents(object = gbsr, verbose = FALSE)

  if (is.null(parents)) {
    hap <- getHaplotype(object = gbsr, parents = FALSE)
    max_hap <- max(hap, na.rm = TRUE)
    if (max_hap == 2) {
      scale_breaks <- seq_len(2)
      scale_labels <- paste0("Hap", seq_len(2))

    } else {
      stop("The number of haplotypes in the data is ", max_hap, ".",
           "\nWhen the number of haplotypes is more than two, ",
           "GenovisR needs parental sample information to properly process",
           " the given genotype and haplotype data.",
           "\nPlease set parental samples using setParents().")
    }

  } else {
    p_hap <- getHaplotype(object = gbsr, parents = "only")
    p_hap <- p_hap[,, 1]
    p_hap <- unique(p_hap)
    n_p_hap <- nrow(p_hap)
    n_p <- ncol(p_hap)
    n_haplotype <- length(unique(as.vector(p_hap)))
    scale_breaks <- seq_len(n_haplotype)

    p_id <- parents$sampleID[match(seq_len(n_p), parents$memberID)]
    rep_p_id <- rep(x = p_id, each = n_p_hap)
    if (n_p_hap == 1) {
      scale_labels <- rep_p_id
    } else {
      scale_labels <- paste0(rep_p_id, "_hap", seq_len(n_p_hap))
    }
  }

  out <- getHaplotype(object = gbsr, parents = "FALSE")
  attributes(out)$scale_breaks <- scale_breaks
  attributes(out)$scale_labels <- scale_labels
  return(out)
}

#' Get Dosage Information from GBSR
#'
#' @param gbsr A GbsrGenotypeData object.
#' @return A dosage matrix with scale breaks and labels as attributes.
.getDosageFromGBSR <- function(gbsr) {
  parents <- getParents(object = gbsr, verbose = FALSE)

  if(is.null(parents)){
    gbsr <- .setParentsFromGDS(gbsr = gbsr)
    parents <- getParents(object = gbsr, verbose = FALSE)
  }

  if (is.null(parents)) {
    hap <- getHaplotype(object = gbsr, parents = FALSE)
    n_ploidy <- dim(hap)[1]
    max_hap <- max(hap, na.rm = TRUE)

  } else {
    p_hap <- getHaplotype(object = gbsr, parents = "only")
    n_ploidy <- dim(p_hap)[1]
    p_hap <- p_hap[,, 1]
    max_hap <- length(unique(as.vector(p_hap)))
  }

  if (max_hap == 2) {
    out <- getGenotype(object = gbsr, node = "dosage", parents = "FALSE")
    attributes(out)$scale_breaks <- seq_len(n_ploidy)
    attributes(out)$scale_labels <- paste0("Plex", seq(0, n_ploidy))

  } else {
    message("The number of haplotypes in the data is ", max_hap, ".",
            "\nDosage information is available only when the number of ",
            "haplotypes is two.",
            "The output genovis object will has NULL for ",
            "the dosage information.")
    out <- NULL
  }

  return(out)
}
