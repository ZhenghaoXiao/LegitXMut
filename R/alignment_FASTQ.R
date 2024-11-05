#' Aligns FASTQ File to Reference Genome and Produces a BAM File, a VCF File and a Summary
#'
#' A function to align a FASTQ file to a reference genome, generating a BAM file.
#' This function is suitable for handling large sequencing datasets and supports
#' alignment of Nanopore reads and quality control.
#'
#' @param fastqPath A character string for the path to the FASTQ file.
#' @param referencePath A character string for the path to the reference genome in FASTA format.
#' @param indels An integer specifying the maximum number of insertions or deletions allowed in the alignment. Default is 10.
#' @param maxMismatches An integer specifying the maximum number of mismatches allowed in the alignment. Default is 1000.
#' @param outputBAM A character string for the name of the output BAM file. Default is "output.bam".
#'
#' @return An object containing alignment information from the BAM file and information of variant calling from VCF file,
#'    which can be used for downstream analysis, including variant calling and visualization.
#'
#' @examples
#' # Example 1: Align a FASTQ file to a reference genome with custom parameters.
#' mutateddemoalignment <- alignment_FASTQ(fastqPath =  "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/SRR29917898.fastq",
#' referencePath =  "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/yeast.fna",
#' outputBAM =  "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/aligned_output.bam")
#'
#' @references
#' Morgan M, Lawrence M, Aboyoun P, Gentleman R, and Pages H (2019).
#' “GenomicAlignments: Efficient alignments processing in R for NGS data.”
#' Bioconductor. doi:10.18129/B9.bioc.GenomicAlignments. Available at: https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html.
#'
#' Liao Y, Smyth GK, Shi W (2019). “The R package Rsubread is easier, faster, cheaper
#' and better for alignment and quantification of RNA sequencing reads.”
#' Nucleic Acids Research, 47, e47. doi:10.1093/nar/gkz114.
#' Available at: https://academic.oup.com/nar/article/47/8/e47/5371636.
#'
#' National Center for Biotechnology Information (NCBI). Saccharomyces cerevisiae S288C Genome Assembly (GCF_000146045.2).
#' NCBI Datasets, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/. Accessed 4 Nov. 2024.
#'
#' National Center for Biotechnology Information (NCBI). Sequence Read Archive (SRA) Run: SRR29917898.
#' NCBI SRA, https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29917898&display=download. Accessed 4 Nov. 2024.
#'
#' OpenAI. ChatGPT: Assistance with R function development for bioinformatics applications,
#' "Assignment". https://chat.openai.com (2023, accessed 5 November 2024 for debugging).
#'
#' @export
#' @import Rsubread
#' @importFrom GenomicAlignments readGAlignments
alignment_FASTQ <- function(fastqPath, referencePath, indels = 10, maxMismatches = 1000, outputBAM = "output.bam") {

  indexBase <- gsub("\\.fasta$", "", referencePath)

  if (!file.exists(paste0(indexBase, ".00.b.array"))) {
    message("Building index for the reference genome...")
    Rsubread::buildindex(basename = indexBase, reference = referencePath)
  }

  message("Aligning reads to the reference genome...")
  Rsubread::align(index = indexBase,
                  readfile1 = fastqPath,
                  input_format = "FASTQ",
                  output_file = outputBAM,
                  nthreads = 1,
                  indels = indels,
                  maxMismatches = maxMismatches)

  message("Alignment complete. BAM file saved as: ", outputBAM)

  bamData <- GenomicAlignments::readGAlignments(outputBAM)
  return(bamData)
}


# [END] written by Zhenghao Xiao
