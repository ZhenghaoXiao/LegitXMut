library(LegitXMut)
test_that("plot_vcf_mutation_data generates a heatmap without errors", {
  vcfPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf"  # Replace with actual test file path
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  expect_silent(plot_vcf_mutation_data(
    vcfPath = vcfPath,
    plotType = "heatmap",
    title = "Heatmap Test",
    font_size = 12,
    xlab = "Chromosome",
    ylab = "Mutation Frequency",
    legend_position = "bottom"
  ))
})

test_that("plot_vcf_mutation_data generates a rainfall plot without errors", {
  vcfPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf"  # Replace with actual test file path
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  expect_silent(plot_vcf_mutation_data(
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

test_that("plot_vcf_mutation_data generates a manhattan plot without errors", {
  vcfPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf"  # Replace with actual test file path
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  expect_silent(plot_vcf_mutation_data(
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
  vcfPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf"  # Replace with actual test file path
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  expect_error(
    plot_vcf_mutation_data(vcfPath = vcfPath, plotType = "invalid_type"),
    "Invalid plot type"
  )
})

test_that("plot_vcf_mutation_data allows custom color schemes for heatmap", {
  vcfPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf"  # Replace with actual test file path
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  custom_colors <- circlize::colorRamp2(breaks = c(0, 10), colors = c("blue", "yellow"))
  expect_silent(plot_vcf_mutation_data(
    vcfPath = vcfPath,
    plotType = "heatmap",
    color_scheme = custom_colors,
    title = "Custom Heatmap Test",
    font_size = 10
  ))
})

test_that("plot_vcf_mutation_data allows custom color schemes for rainfall plot", {
  vcfPath <- "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf"  # Replace with actual test file path
  skip_if_not(file.exists(vcfPath), "Test VCF file not found")

  custom_colors <- c("C>A" = "purple", "C>G" = "cyan", "C>T" = "pink",
                     "T>A" = "brown", "T>C" = "green", "T>G" = "orange",
                     "indel" = "grey")

  expect_silent(plot_vcf_mutation_data(
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
