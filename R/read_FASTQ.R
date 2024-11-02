#' Reads FASTQ Files with Processing in Parts and Optional Sampling
#'
#' A function to read FASTQ files, either by processing in memory-efficient parts
#' or by sampling a specified number of records. This function is suitable for
#' handling large sequencing datasets in bioinformatics.
#'
#' @param filePath A character string for the path to the FASTQ file.
#' @param Size An integer indicating the number of records per region for processing.
#'    Typically useful while processing large data such as genome. Default is 1,000,000.
#' @param sampleSize An integer specifying the number of random records to sample from the file.
#'    If NULL, the function will process the entire file. Default is NULL.
#'
#' @return An object containing either the sampled or fully processed FASTQ data.
#'   Including read sequences, quality scores, and read IDs.
#'
#' @examples
#' # Example 1: Process entire FASTQ file in parts, replace it with local path.
#' fullData <- read_FASTQ(
#'   filePath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/data/sequence_demo.fastq",
#'   Size = 1e6
#' )
#'
#' # Example 2: Sample 20,000 random records from the FASTQ file, replace it with local path.
#' sampledData <- read_FASTQ(
#'   filePath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/data/sequence_demo.fastq",
#'   sampleSize = 20000
#' )
#'
#' @references
#' Morgan M, Anders S, Lawrence M, Aboyoun P, Pagès H, Gentleman R (2009).
#' “ShortRead: a Bioconductor package for input, quality assessment and exploration
#' of high-throughput sequence data.” Bioinformatics, 25, 2607-2608.
#' doi:10.1093/bioinformatics/btp450, \href{http://dx.doi.org10.1093/bioinformatics/btp450}{Link}.
#'
#' NCBI. TP53 gene (Gene ID: 7157), Homo sapiens (human). NCBI Gene Database.
#' National Center for Biotechnology Information (NCBI).
#' Available at: https://www.ncbi.nlm.nih.gov/gene/7157. Accessed on: October 31, 2024.
#'
#' OpenAI. ChatGPT: Debugging assistance and suggestions in use of R package functions,
#' "Assignment". https://chat.openai.com (2023, accessed 5 November 2024).
#'
#' @export
#' @import ShortRead
read_FASTQ <- function(filePath, Size = 1e6, sampleSize = NULL) {
  if (!is.null(sampleSize)) {
    sample <- ShortRead::FastqSampler(filePath, sampleSize)
    sampledData <- ShortRead::yield(sample)
    return(sampledData)
  }

  stream <- ShortRead::FastqStreamer(filePath, Size)
  on.exit(close(stream))

  fastqparts <- list()

  repeat {
    fqpart <- ShortRead::yield(stream)
    if (length(fqpart) == 0) break
    fastqparts <- c(fastqparts, list(fqpart))
  }

  Fastq <- do.call(c, fastqparts)
  return(Fastq)
}

# [END] written by Zhenghao Xiao
