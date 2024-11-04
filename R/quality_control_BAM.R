#' Filters a BAM File Based on Quality Control Metrics
#'
#' This function filters a BAM file based on specified quality control criteria, such as
#' minimum read length and minimum mapping quality. It excludes reads that do not meet
#' these requirements, ensuring that only high-quality reads are retained for downstream
#' analyses, such as variant calling.
#'
#' @param bamPath A character string specifying the path to the input BAM file.
#' @param min_length An integer indicating the minimum read length to retain.
#'    Reads shorter than this length will be excluded. Default is 50.
#' @param min_depth An integer specifying the minimum mapping quality (depth) required for each read.
#'    Reads with a depth lower than this value will be excluded. Default is 10.
#' @param outputBAM A character string for the name of the output filtered BAM file. Default is "filtered_output.bam".
#'
#' @return A filtered BAM file that contains only high-quality reads meeting the specified criteria.
#'    This output BAM file can be used for further analyses, such as variant calling and error rate calculation.
#'
#' @examples
#' # Example 1: Filter a BAM file based on minimum length and depth requirements, replace it with local path.
#' filteredResult <- quality_control_BAM(
#'   bamPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/aligned_output.bam",
#'   min_length = 150,
#'   min_depth = 1,
#'   outputBAM = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/filtered_output.bam"
#' )
#'
#' @references
#' Morgan M, Pagès H, Obenchain V, Hayden N (2024). Rsamtools: Binary
#' alignment (BAM), FASTA, variant call (BCF), and tabix file import. R package
#' version 2.20.0, https://bioconductor.org/packages/Rsamtools.
#'
#' NCBI. TP53 gene (Gene ID: 7157), Homo sapiens (human). NCBI Gene Database.
#' National Center for Biotechnology Information (NCBI).
#' Available at: https://www.ncbi.nlm.nih.gov/gene/7157. Accessed on: October 31, 2024.
#'
#' NCBI. Trace Archive, Run Browser for SRR4000542. National Center for Biotechnology Information (NCBI).
#' Available at: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR4000542&display=metadata. Accessed on: October 31, 2024.
#'
#' OpenAI. ChatGPT: Assistance with R function development for bioinformatics applications,
#' "Assignment". https://chat.openai.com (2023, accessed 5 November 2024).
#'
#' @export
#' @import Rsamtools
quality_control_BAM <- function(bamPath, min_length = 50, min_depth = 10, outputBAM = "filtered_output.bam") {

  bamFile <- Rsamtools::BamFile(bamPath)
  message("Opening BAM file: ", bamPath)

  # Retrieve read lengths and mapping quality from BAM file
  bam_reads <- Rsamtools::scanBam(bamFile, param = Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE),
    what = c("qwidth", "mapq")
  ))[[1]]

  # Select indices of reads that meet filtering criteria
  selected_indices <- which(bam_reads$qwidth >= min_length & bam_reads$mapq >= min_depth)

  # Proceed if there are reads that meet the criteria
  if (length(selected_indices) > 0) {
    # Filter BAM by copying only selected reads based on criteria
    filter_param <- Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
    )
    Rsamtools::filterBam(bamFile, destination = outputBAM, param = filter_param)
    message("Quality control complete. Filtered BAM file saved as: ", outputBAM)
  } else {
    message("No reads met the specified quality criteria. No output BAM file created.")
  }

  return(invisible(outputBAM))
}

# [END] written by Zhenghao Xiao
