library(LegitXMut)
test_that("plot_vcf_mutation_data generates a heatmap without errors", {
  vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut") # Path to VCF
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  expect_no_error(plot_vcf_mutation_data(
    vcfPath = vcfPath,
    plotType = "heatmap",
    title = "Heatmap Test",
    font_size = 12,
  ))
})

test_that("plot_vcf_mutation_data generates a rainfall plot without errors", {
  # Locate the packaged test data
  vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")  # Path to VCF
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  expect_no_error(plot_vcf_mutation_data(
    vcfPath = vcfPath,
    plotType = "rainfall",
    title = "Rainfall Plot Test",
    color_scheme = c("C>A" = "red", "C>G" = "orange", "C>T" = "green",
                     "T>A" = "yellow", "T>C" = "blue", "T>G" = "purple",
                     "indel" = "grey"),
    alpha = 0.7,
    font_size = 10,
    xlab = "Genomic Position",
    ylab = "Inter-Variant Distance (log10)",
    legend_position = "right"
  ))
})

test_that("plot_vcf_mutation_data generates a Manhattan plot without errors", {
  # Locate the packaged test data
  vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")  # Path to VCF
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  expect_no_error(plot_vcf_mutation_data(
    vcfPath = vcfPath,
    plotType = "manhattan",
    title = "Manhattan Plot Test",
    ylab = "-log10(Supporting Reads)",
    xlab = "Genomic Position",
    font_size = 11,
    alpha = 0.6,
    legend_position = "top"
  ))
})

test_that("plot_vcf_mutation_data handles invalid plot type", {
  # Locate the packaged test data
  vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")  # Path to VCF
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  expect_error(
    plot_vcf_mutation_data(vcfPath = vcfPath, plotType = "invalid_type"),
    "Invalid plot type"
  )
})

test_that("plot_vcf_mutation_data allows custom color schemes for heatmap", {
  # Locate the packaged test data
  vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")  # Path to VCF
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  custom_colors <- circlize::colorRamp2(breaks = c(0, 10), colors = c("blue", "yellow"))
  expect_no_error(plot_vcf_mutation_data(
    vcfPath = vcfPath,
    plotType = "heatmap",
    color_scheme = custom_colors,
    title = "Custom Heatmap Test",
    font_size = 10
  ))
})

test_that("plot_vcf_mutation_data allows custom color schemes for rainfall plot", {
  # Locate the packaged test data
  vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")  # Path to VCF
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  custom_colors <- c("C>A" = "purple", "C>G" = "cyan", "C>T" = "pink",
                     "T>A" = "brown", "T>C" = "green", "T>G" = "orange",
                     "indel" = "grey")

  expect_no_error(plot_vcf_mutation_data(
    vcfPath = vcfPath,
    plotType = "rainfall",
    color_scheme = custom_colors,
    title = "Custom Rainfall Plot Test",
    alpha = 0.7,
    font_size = 10,
    xlab = "Genomic Position",
    ylab = "Inter-Variant Distance (log10)",
    legend_position = "right"
  ))
})

test_that("plot_vcf_mutation_data handles invalid alpha parameter", {
  # Locate the packaged test data
  vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  # Invalid alpha: negative value
  expect_error(
    plot_vcf_mutation_data(
      vcfPath = vcfPath,
      plotType = "rainfall",
      title = "Negative Alpha Test",
      alpha = -0.5,  # Invalid alpha
      font_size = 10
    ),
    "Alpha must be a numeric value between 0 and 1."  # Error message match
  )

  # Invalid alpha: greater than 1
  expect_error(
    plot_vcf_mutation_data(
      vcfPath = vcfPath,
      plotType = "rainfall",
      title = "Alpha Greater Than 1 Test",
      alpha = 1.5,  # Invalid alpha
      font_size = 10
    ),
    "Alpha must be a numeric value between 0 and 1."  # Error message match
  )
})

test_that("plot_vcf_mutation_data throws error for incorrect file type", {
  # Use a file path with an invalid file type (not a VCF file)
  incorrectFilePath <- system.file("extdata", "not_a_vcf.txt", package = "LegitXMut") # Invalid file
  skip_if_not(file.exists(incorrectFilePath), "Incorrect file type not found")

  expect_error(
    plot_vcf_mutation_data(
      vcfPath = incorrectFilePath,
      plotType = "manhattan",
      title = "Invalid File Type Test",
      alpha = 0.6
    ),
    "The specified VCF file does not exist or is invalid"
  )
})
