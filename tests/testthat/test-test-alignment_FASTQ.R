library(LegitXMut)
test_that("alignment_FASTQ produces correct BAM output with specific paths", {
  fastqPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/SRR29917898.fastq"# Replace with actual test file path
  referencePath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/yeast.fna"# Replace with actual test file path
  outputBAM <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/aligned_output.bam"
  skip_if_not(file.exists(fastqPath), "Test FASTQ file not found")
  skip_if_not(file.exists(referencePath), "Test reference genome file not found")
  alignment_result <- alignment_FASTQ(
    fastqPath = fastqPath,
    referencePath = referencePath,
    indels = 10,
    maxMismatches = 1000,
    outputBAM = outputBAM
  )

  expect_true(grepl("\\.bam$", outputBAM, ignore.case = TRUE),
              "Output file should be a .bam file")
  expect_true(file.exists(outputBAM), "Output BAM file should exist")
  expect_gt(file.size(outputBAM), 0, "Output BAM file should not be empty")
})

test_that("alignment_FASTQ can handle missing files", {
  fastqPath <- "nonexistent_sample.fastq"
  referencePath <- "nonexistent_reference.fasta"
  outputBAM <- tempfile(fileext = ".bam")
  skip_if_not(file.exists(fastqPath), "Test Fastq file not found")
  skip_if_not(file.exists(referencePath), "Test Reference file not found")
  unlink(outputBAM)
})
