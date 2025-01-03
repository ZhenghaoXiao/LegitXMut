library(LegitXMut)
test_that("update_vcf correctly updates chromosome names in VCF file", {
  # Use system.file to locate the packaged test data
  fastaPath <- system.file("extdata", "yeast.fna", package = "LegitXMut") # Reference genome path
  vcfPath <- system.file("extdata", "aligned_output.bam.indel.vcf", package = "LegitXMut") # Test VCF file path
  outputVcfPath <- tempfile(fileext = ".vcf") # Temporary path for updated VCF file

  # Skip test if files are not present
  skip_if_not(file.exists(fastaPath), "Test FASTA file not found")
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  update_vcf(fastaPath, vcfPath, outputVcfPath)# Run the update_vcf function

  expect_true(file.exists(outputVcfPath), info = "Output VCF file should be created.")# Validate that the output VCF file is created

  # Read the updated VCF file and check sequence levels
  output_vcf <- VariantAnnotation::readVcf(outputVcfPath)
  output_seqlevels <- GenomeInfoDb::seqlevels(output_vcf)

  # Validate chromosome names
  expect_true("chrI" %in% output_seqlevels, info = "Chromosome 'chrI' should be present in output.")
  expect_true("chrII" %in% output_seqlevels, info = "Chromosome 'chrII' should be present in output.")

  unlink(outputVcfPath)
})

test_that("update_vcf handles cases with no matching chromosomes gracefully", {

  # Use an alternate FASTA file with non-matching chromosome names for testing
  fastaPath <- system.file("extdata", "non_matching.fna", package = "LegitXMut") # Use a FASTA file with different chromosome names
  vcfPath <- system.file("extdata", "aligned_output.bam.indel.vcf", package = "LegitXMut") # Test VCF file path
  outputVcfPath <- tempfile(fileext = ".vcf") # Temporary path for updated VCF file

  # Skip test if files are not present
  skip_if_not(file.exists(fastaPath), "Test Fasta file not found")
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  unlink(outputVcfPath)
})

test_that("update_vcf throws error for incorrect file types", {
  # Define paths to incorrect file types
  incorrectFastaPath <- system.file("extdata", "not_a_fasta.txt", package = "LegitXMut") # Invalid file type for FASTA
  vcfPath <- system.file("extdata", "aligned_output.bam.indel.vcf", package = "LegitXMut") # Valid VCF file path
  incorrectVcfPath <- system.file("extdata", "not_a_vcf.txt", package = "LegitXMut") # Invalid file type for VCF
  outputVcfPath <- tempfile(fileext = ".vcf") # Temporary path for updated VCF file

  # Skip test if files are not present
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  # Test with incorrect FASTA file type
  expect_error(
    update_vcf(incorrectFastaPath, vcfPath, outputVcfPath),
    "The specified FASTA file does not exist"
  )

  # Test with incorrect VCF file type
  expect_error(
    update_vcf(fastaPath = system.file("extdata", "yeast.fna", package = "LegitXMut"),
               vcfPath = incorrectVcfPath,
               outputVcfPath = outputVcfPath),
    "The specified VCF file does not exist"
  )

  unlink(outputVcfPath)
})
