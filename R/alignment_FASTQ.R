#' Aligns FASTQ File to Reference Genome and Produces a BAM File, a VCF File and a Summary
#'
#' A function to align a FASTQ file to a reference genome, generating a BAM and a VCF file.
#' This function is suitable for handling large sequencing datasets and supports
#' alignment of Nanopore reads and quality control.
#'
#' @param fastqPath A character string for the path to the FASTQ file.
#' @param referencePath A character string for the path to the reference genome in FASTA format.
#' @param indels An integer specifying the maximum number of insertions or deletions allowed in the alignment. Default is 10.
#' @param maxMismatches An integer specifying the maximum number of mismatches allowed in the alignment. Default is 1000.
#' @param outputBAM A character string for the name and path of the output BAM and VCF file. Default is "output.bam". VCF file follows this name
#'
#' @return An object containing alignment information from the BAM file and information of variant calling from VCF file,
#'    which can be used for downstream analysis, including variant calling and visualization.
#'
#' @details
#' The LegitXMut package includes external demo data stored in the extdata directory,
#'    and these external files are provided for demonstration purposes and support reproducibility of the analysis.
#'
#' @examples
#' # Example 1: Align a FASTQ file to a reference genome with custom parameters.
#' # Replace with actual test file path
#' fastqPath <- system.file("extdata", "ERR12205202.fastq", package = "LegitXMut")
#' referencePath <- system.file("extdata", "yeast.fna", package = "LegitXMut")
#' # Locate the extdata directory for outputs
#' outputDir <- system.file("extdata", package = "LegitXMut")
#' outputBAM <- file.path(outputDir, "aligned_output.bam")
#'
#' mutateddemoalignment <- alignment_FASTQ(fastqPath =  fastqPath,
#' referencePath =  referencePath,
#' indels = 10,
#' maxMismatches = 1000,
#' outputBAM =  outputBAM)
#'
#' @references
#' Liao Y, Smyth GK, Shi W (2019)."The R package Rsubread is easier, faster, cheaper
#' and better for alignment and quantification of RNA sequencing reads."
#' Nucleic Acids Research, 47, e47. doi:10.1093/nar/gkz114.
#' Available at: https://academic.oup.com/nar/article/47/8/e47/5371636.
#'
#' Morgan M, Lawrence M, Aboyoun P, Gentleman R, and Pages H (2019).
#' "GenomicAlignments: Efficient alignments processing in R for NGS data"
#' Bioconductor.doi:10.18129/B9.bioc.GenomicAlignments.Available at: https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html.
#'
#' National Center for Biotechnology Information (NCBI). Saccharomyces cerevisiae S288C Genome Assembly (GCF_000146045.2).
#' NCBI Datasets, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/. Accessed 4 Nov. 2024.
#'
#' National Center for Biotechnology Information (NCBI). Sequence Read Archive (SRA) Run: ERR12205202.
#' NCBI SRA, https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR12205202&display=download. Accessed 4 Nov. 2024.
#'
#' OpenAI. ChatGPT: Assistance with R function development for bioinformatics applications,
#' "Assignment".https://chat.openai.com (2023, accessed 5 November 2024 for debugging).
#'
#' @export
#' @importFrom Rsubread buildindex align
#' @importFrom GenomicAlignments readGAlignments
alignment_FASTQ <- function(fastqPath, referencePath, indels = 10, maxMismatches = 1000, outputBAM = "output.bam") {
  # --- INPUT VALIDATION ---
  if (!file.exists(fastqPath)) {
    stop("FASTQ file does not exist: ", fastqPath)# FASTQ validation
  }
  if (!file.exists(referencePath)) {
    stop("Reference genome file does not exist: ", referencePath)# Reference Genome validation
  }
  if (!is.numeric(indels) || indels < 0) {
    stop("Invalid 'indels' value. Must be a non-negative integer.")# Indel value validation
  }
  if (!is.numeric(maxMismatches) || maxMismatches < 0) {
    stop("Invalid 'maxMismatches' value. Must be a non-negative integer.")# Max mismatches validation
  }
  # --- PREPARATION---
  # Remove the ".fasta" extension from the reference genome path to extract and use the index base names
  indexBase <- gsub("\\.fasta$", "", referencePath)
  # Check if the indexs are already there,If not, build the index using Rsubread
  if (!file.exists(paste0(indexBase, ".00.b.array"))) {
    message("Building index for the reference genome...")
    Rsubread::buildindex(basename = indexBase, reference = referencePath)
  }

  # --- ALIGNMENT---
  #Align the FASTQ file to the reference genome
  message("Aligning reads to the reference genome...")
  Rsubread::align(index = indexBase,# Reference genome index base
                  readfile1 = fastqPath,# Input FASTQ file
                  input_format = "FASTQ",# Specify input file format
                  output_file = outputBAM,# Output BAM file path
                  nthreads = 1,# Number of threads for the alignment
                  indels = indels,# Maximum number of insertions/deletions allowed
                  maxMismatches = maxMismatches)# Maximum number of mismatches allowed
  # Confirmation of alignment completion
  message("Alignment complete. BAM file saved as: ", outputBAM)
  message("Alignment complete. VCF file has also been saved to this directory")

  bamData <- GenomicAlignments::readGAlignments(outputBAM)
  # --- RETURN---
  return(bamData)
}

# [END] written by Zhenghao Xiao
