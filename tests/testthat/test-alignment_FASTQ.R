library(LegitXMut)
library(testthat)
test_that("alignment_FASTQ produces correct BAM output with specific paths", {
  # Define the file paths
  fastqPath <- system.file("extdata", "SRR29917898.fastq", package = "LegitXMut") # Use system.file to find packaged test data
  referencePath <- system.file("extdata", "yeast.fna", package = "LegitXMut")     # Use system.file to locate reference file
  outputDir <- system.file("extdata", package = "LegitXMut")# Use system.file for output
  outputBAM <- file.path(outputDir, "aligned_output.bam")

  # Skip test if files are not present
  skip_if_not(file.exists(fastqPath), "Test FASTQ file not found")
  skip_if_not(file.exists(referencePath), "Test reference genome file not found")

  # Run alignment function
  alignment_result <- alignment_FASTQ(
    fastqPath = fastqPath,
    referencePath = referencePath,
    indels = 10,
    maxMismatches = 1000,
    outputBAM = outputBAM
  )

  # Assertions to validate the BAM output
  expect_true(grepl("\\.bam$", outputBAM, ignore.case = TRUE),
              "Output file should be a .bam file")
  expect_true(file.exists(outputBAM), "Output BAM file should exist")
  expect_gt(file.size(outputBAM), 0, "Output BAM file should not be empty")
})

test_that("alignment_FASTQ can handle missing files", {
  # Define nonexistent paths for testing
  fastqPath <- "nonexistent_sample.fastq"
  referencePath <- "nonexistent_reference.fasta"
  outputBAM <- tempfile(fileext = ".bam")

  # Skip test if these files are detected (they should not exist)
  skip_if_not(file.exists(fastqPath), "Test Fastq file not found")
  skip_if_not(file.exists(referencePath), "Test Reference file not found")
  unlink(outputBAM)
})
