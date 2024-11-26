#' Update Chromosome Names to VCF File Based on Reference FASTA If Performing Single Chromosome Anlysis
#'
#' This function reads a reference FASTA file to extract chromosome mappings and updates
#' the chromosome names in a VCF file accordingly. The updated VCF file is saved to a specified output path.
#'
#' @param fastaPath A character string specifying the path to the reference FASTA file.
#' @param vcfPath A character string specifying the path to the input VCF file that needs chromosome name updating.
#' @param outputVcfPath A character string specifying the path where the updated VCF file will be saved.
#'
#' @return NULL. A message confirms that the VCF file has been updated and saved.
#'
#' @examples
#' # Example: Update chromosome names in a VCF file using a reference FASTA file/FNA file
#' # Replace with actual test file path
#' fastaPath <- system.file("extdata", "yeast.fna", package = "LegitXMut")
#' vcfPath <- system.file("extdata", "aligned_output.bam.indel.vcf", package = "LegitXMut")
#' outputVcfPath <- file.path(system.file("extdata", package = "LegitXMut"), "updated.vcf")
#'
#' # Update the VCF file with specified inputs
#' updatingvcf <- update_vcf(
#'   fastaPath = fastaPath,
#'   vcfPath = vcfPath,
#'   outputVcfPath = outputVcfPath
#' )
#'
#' @references
#' Morgan M, Pages H, Obenchain V, and Hayden N (2020).
#' “VariantAnnotation: annotation of genetic variants.”
#' Bioconductor. doi:10.18129/B9.bioc.VariantAnnotation. Available at: https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html.
#'
#' Wickham, H. (2019). stringr: Simple, Consistent Wrappers for Common String Operations.
#' R package version 1.4.0. Available at: https://CRAN.R-project.org/package=stringr
#'
#' Pagès, H., Aboyoun, P., DebRoy, S., and Lawrence, M. (2023).
#' GenomeInfoDb: Utilities for manipulating chromosome names, including modifying the names to follow a particular convention.
#' R package version 1.36.0, Bioconductor, https://bioconductor.org/packages/GenomeInfoDb.
#'
#' OpenAI. ChatGPT: Assistance with R function development for bioinformatics applications,
#' "Assignment". https://chat.openai.com (2023, accessed 5 November 2024 for debugging).
#'
#' National Center for Biotechnology Information (NCBI). Saccharomyces cerevisiae S288C Genome Assembly (GCF_000146045.2).
#' NCBI Datasets, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/. Accessed 4 Nov. 2024.
#'
#' National Center for Biotechnology Information (NCBI). Sequence Read Archive (SRA) Run: ERR12205202.
#' NCBI SRA, https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR12205202&display=download. Accessed 4 Nov. 2024.
#'
#' @export
#' @importFrom VariantAnnotation readVcf writeVcf
#' @import stringr
#' @importFrom GenomeInfoDb seqlevels
update_vcf <- function(fastaPath, vcfPath, outputVcfPath) {

  chrom_mapping <- list()
  fasta_con <- file(fastaPath, "r")

  while (TRUE) {
    line <- readLines(fasta_con, n = 1)
    if (length(line) == 0) break
    if (startsWith(line, ">")) {
      match <- stringr::str_match(line, "^>(\\S+) .*chromosome ([^, ]+)")
      if (!is.na(match[1, 2]) && !is.na(match[1, 3])) {
        ref_name <- match[1, 2]
        chrom_name <- paste0("chr", match[1, 3])
        chrom_mapping[[ref_name]] <- chrom_name
      }
    }
  }
  close(fasta_con)

  vcf <- VariantAnnotation::readVcf(vcfPath)

  current_seqlevels <- GenomeInfoDb::seqlevels(vcf)
  for (old_name in names(chrom_mapping)) {
    new_name <- chrom_mapping[[old_name]]
    if (old_name %in% current_seqlevels) {
      GenomeInfoDb::seqlevels(vcf)[GenomeInfoDb::seqlevels(vcf) == old_name] <- new_name
    }
  }

  VariantAnnotation::writeVcf(vcf, outputVcfPath)

  message("VCF chromosome names updated and saved to: ", outputVcfPath)
}

# [END] written by Zhenghao Xiao
