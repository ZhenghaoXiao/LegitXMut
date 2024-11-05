library(LegitXMut)
test_that("update_vcf correctly updates chromosome names in VCF file", {

  fastaPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/yeast.fna" # Replace with actual test file path
  vcfPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/aligned_output.bam.indel.vcf" # Replace with actual test VCF file path
  outputVcfPath <- tempfile(fileext = ".vcf")

  update_vcf(fastaPath, vcfPath, outputVcfPath)

  expect_true(file.exists(outputVcfPath), info = "Output VCF file should be created.")

  output_vcf <- VariantAnnotation::readVcf(outputVcfPath)
  output_seqlevels <- GenomeInfoDb::seqlevels(output_vcf)

  # Assuming "chrI" and "chrII" are example chromosome names in the reference genome
  expect_true("chrI" %in% output_seqlevels, info = "Chromosome 'chrI' should be present in output.")
  expect_true("chrII" %in% output_seqlevels, info = "Chromosome 'chrII' should be present in output.")

  unlink(outputVcfPath)
})

test_that("update_vcf handles cases with no matching chromosomes gracefully", {

  fastaPath <- "path/to/test_reference.fasta"   # Use a FASTA file with chromosome names not in the VCF
  vcfPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/aligned_output.bam.indel.vcf" # Replace with actual test VCF file path
  outputVcfPath <- tempfile(fileext = ".vcf")

  skip_if_not(file.exists(fastaPath), "Test Fasta file not found")
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  unlink(outputVcfPath)
})
