#' Import markers, genotypes, and haplotypes information from a GDS file
#'
#' Chromosome IDs, physical positions, and alleles of markers stored in an input
#' GDS file are retrieved as well as genotypes and haplotypes information at the
#' markers. The obtained information is then organized as a genovis class
#' object that can be handled by the functions in the genovisr package.
#'
#' @param gds_fn A string to specify a path to the input GDS file.
#' @param gbscleanr_filter A logical value to indicate whether the filtering
#' information stored in the input GDS file should be loaded to apply the
#' filtering on markers and samples. See [GBScleanR::loadGDS()] and
#' [GBScleanR::closeGDS()].
#' @param ploidy Not implemented yet.
#' @param verbose A logical value to indicate whether to show messages from
#' the function.
#'
#' @importClassesFrom GBScleanR GbsrGenotypeData
#' @importMethodsFrom GBScleanR loadGDS getGenotype getHaplotype getRead
#' @importFrom gdsfmt exist.gdsn
#'
#' @examples
#' # Not run
#' # genovis <- gdscleanr2genovis(gds_fn = "path/to/gds/file.gds",
#' #                              gbscleanr_filter = TRUE)
#'
#' @export
#'
gbscleanr2genovis <- function(gds_fn,
                              gbscleanr_filter = TRUE,
                              ploidy = 2,
                              verbose = TRUE) {
  if(verbose){ message('Loading GDS file.') }

  # Check if the input is a string indicating a file path
  if(!is.character(gds_fn)){
    stop("A path to a GDS file should be specified to 'gds_fn'.")
  }

  # Load the input GDS file using the loadGDS() function
  # of the GBScleanR package.
  gds <- loadGDS(x = gds_fn,
                 load_filter = gbscleanr_filter,
                 ploidy = ploidy,
                 verbose = verbose)

  # Pull out data required to build a genovis class object.
  # Get marker allele information.
  allele <- getAllele(object = gds)
  allele <- strsplit(x = allele, split = ",")
  allele <- do.call("rbind", allele)

  # Confirm whether corrected genotype data is available or not.
  if(exist.gdsn(node = gds$root, path = "annotation/format/CGT")){
    genotype_node <- "cor"

  } else {
    genotype_node <- "raw"
  }

  # Organize marker information and genotype information in an output object.
  chr <- getChromosome(object = gds)
  pos <- getPosition(object = gds)
  out <- list(marker_info = data.frame(id = paste(chr, pos, sep = "_"),
                                       chr = factor(chr),
                                       pos = as.integer(pos),
                                       ref_allele = allele[, 1],
                                       alt_allele = allele[, 2]),
              sample_info = data.frame(id = getSamID(object = gds)),
              genotype = getGenotype(object = gds,
                                     node = genotype_node,
                                     parents = TRUE),
              haplotype = NULL,
              dosage = NULL)

  # Add haplotype information if available.
  if(exist.gdsn(node = gds$root, path = "annotation/format/HAP")){
    out$haplotype  <- getHaplotype(object = gds, parents = TRUE)

  }

  if(exist.gdsn(node = gds$root, path = "annotation/format/DS")){
    out$haplotype  <- getGenotype(object = gds, node = "dosage", parents = TRUE)
  }

  # Set the genovis class to the output.
  class(out) <- c(class(out), "genovis")

  # return the genovis class object.
  return(out)
}

