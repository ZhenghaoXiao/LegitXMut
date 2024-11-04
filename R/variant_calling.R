#' Calls Variants in a Filtered BAM File and Outputs VCF File
#'
#' This function takes a filtered BAM file and a reference genome in FASTA format to perform
#' variant calling. Detected variants, including SNPs, insertions, and deletions, are saved to a VCF file.
#' It utilizes VariantAnnotation and VariantTools packages for processing.
#'
#' @param bamfilteredPath A character string specifying the path to the filtered BAM file.
#' @param referencePath A character string specifying the path to the reference genome in FASTA format.
#' @param outputVCF A character string for the output VCF file name. Default is "output_variants.vcf".
#'
#' @return The path to the generated VCF file containing the identified variants.
#'
#' @examples
#' # Example: Perform variant calling on a filtered BAM file
#' variantResult <- variant_calling(
#'   bamPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/filtered_output.bam",
#'   refPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/sequence_demo_ref.fasta",
#'   outputVCF = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/output_variants.vcf"
#' )
#'
#' @references
#' Obenchain V, Lawrence M, Carey V, Gogarten S, Shannon P, Morgan M (2014).
#' “VariantAnnotation: a Bioconductor package for exploration and annotation of
#' genetic variants.” Bioinformatics, 30(14), 2076-2078. doi:10.1093/bioinformatics/btu168.
#'
#' Lawrence M, Degenhardt J, Gentleman R (2024). VariantTools: Tools for
#' Exploratory Analysis of Variant Calls. R package version 1.46.0.
#'
#' @export
#' @import VariantAnnotation
#' @import VariantTools
#' @import Rsamtools
variant_calling <- function(bamPath, refPath, outputVCF = "output_variants.vcf") {

  ref_genome <- FaFile(refPath)
  indexFa(ref_genome)

  bam_file <- BamFile(bamPath)
  param <- TallyVariantsParam(ref_genome, high_base_quality = 20L)

  filters <- VariantCallingFilters(read.count = 2L, p.error = 1/1000)

  message("Performing variant calling on BAM file...")
  variants <- callVariants(x = bam_file, param = param, calling.filters = filters)

  message("Writing identified variants to VCF file: ", outputVCF)
  vcf <- asVCF(variants)
  writeVcf(vcf, filename = outputVCF)

  message("Variant calling complete. VCF file saved as: ", outputVCF)

  return(outputVCF)
}

# [END] written by Zhenghao Xiao
