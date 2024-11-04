#' Visualize Mutation Data from VCF File in Various Plot Styles
#'
#' This function takes a VCF file and generates different mutation visualization plots,
#' including karyograms, density plots, lollipop plots, and circos plots. It provides
#' genome-specific and mutation-focused visualizations directly from the VCF file.
#'
#' @param vcfPath A character string specifying the path to the VCF file.
#' @param plotType A character string specifying the type of plot to generate.
#'    Options are "karyogram" (karyogram of variants), "density" (variant density along chromosomes),
#'    "lollipop" (lollipop plot for specific genes), or "circos" (circular genome plot).
#' @param gene A character string specifying the gene name for the lollipop plot.
#'    Required if `plotType` is "lollipop".
#'
#' @return A plot visualizing mutation information based on the selected plot type.
#'
#' @examples
#' # Example 1: Generate a karyogram of variants
#' plot_vcf_mutation_data(
#'   vcfPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf",
#'   plotType = "karyogram"
#' )
#'
#' # Example 2: Generate a density plot of variants across chromosomes
#' plot_vcf_mutation_data(
#'   vcfPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf",
#'   plotType = "density"
#' )
#'
#' # Example 3: Generate a lollipop plot for mutations in the TP53 gene
#' plot_vcf_mutation_data(
#'   vcfPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf",
#'   plotType = "lollipop",
#'   gene = "TP53"
#'   title = "Lollipop Plot of Mutations in TP53 Gene"
#' )
#'
#' # Example 4: Generate a circular plot of whole-genome mutations
#' plot_vcf_mutation_data(
#'   vcfPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf",
#'   plotType = "circos"
#'   title = "Circos Plot of Whole-Genome Variants"
#' )
#'
#' @references
#' Morgan M, Pages H, Obenchain V, and Hayden N (2020).
#' “VariantAnnotation: annotation of genetic variants.”
#' Bioconductor. doi:10.18129/B9.bioc.VariantAnnotation. Available at: https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html.
#'
#' Yin T, Cook D, Lawrence M (2012). “ggbio: an R package for extending
#'  visualization capabilities of the popular ggplot2 package to genomic data.”
#' Genome Biology, 13(8):R77. doi:10.1186/gb-2012-13-8-r77.
#'
#' Wang, L. et al. (2017). “trackViewer: a Bioconductor package for interactive
#' and integrative visualization of multi-omics data.”
#' BMC Bioinformatics 18, 124. doi:10.1186/s12859-017-1542-4.
#'
#' Gu, Z. et al. (2014). “circlize implements and enhances circular visualization in R.”
#' Bioinformatics, 30(19), 2811–2812. doi:10.1093/bioinformatics/btu393.
#'
#' OpenAI. ChatGPT: Assistance with R function development for bioinformatics applications,
#' "Assignment". https://chat.openai.com (2023, accessed 5 November 2024).
#'
#' @export
#' @import VariantAnnotation
#' @import ggbio
#' @import trackViewer
#' @import circlize
#' @import GenomicRanges
plot_vcf_mutation_data <- function(vcfPath, plotType = "karyogram", gene = NULL, title = "Plot") {
  # Load the VCF file
  vcf <- VariantAnnotation::readVcf(vcfPath, genome = "hg38")

  # Extract row ranges as individual components
  row_ranges <- SummarizedExperiment::rowRanges(vcf)

  # Convert seqnames to character and ensure it matches length of ranges
  seqnames <- as.character(GenomicRanges::seqnames(row_ranges))
  ranges <- GenomicRanges::ranges(row_ranges)

  # Check and adjust seqnames length if necessary
  if (length(seqnames) == 1 && length(ranges) > 1) {
    seqnames <- rep(seqnames, length(ranges))
  } else if (length(seqnames) != length(ranges)) {
    stop("Length of 'seqnames' does not match the number of ranges.")
  }

  # Construct GRanges using validated seqnames and ranges
  gr <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = ranges,
    strand = "*",
    REF = mcols(row_ranges)$REF,
    ALT = mcols(row_ranges)$ALT
  )

  # Generate plot based on plotType
  if (plotType == "karyogram") {
    p <- ggbio::autoplot(gr, layout = "karyogram", aes(color = seqnames)) +
      labs(title = "Karyogram of Variants")
    print(p)

  } else if (plotType == "density") {
    p <- ggbio::autoplot(gr, stat = "coverage", geom = "line") +
      labs(title = "Variant Density Across Chromosomes")
    print(p)

  } else if (plotType == "lollipop") {
    if (is.null(gene)) {
      stop("Please specify a gene name for the lollipop plot.")
    }

    gr_gene <- subset(gr, seqnames == gene)
    trackViewer::lolliplot(gr_gene, gr, ylab = "Mutation Frequency")

  } else if (plotType == "circos") {
    chrom_data <- as.data.frame(gr)

    circlize::circos.clear()

    circlize::circos.initializeWithIdeogram(
      species = "hg38",    # Adding human ideograms if the genome is hg38
      plotType = c("ideogram", "labels")   # Include ideograms and chromosome labels
    )

    circlize::circos.track(
      factors = chrom_data$seqnames,
      x = chrom_data$start,
      y = chrom_data$start,
      panel.fun = function(x, y) {
        circos.points(x, y, col = "blue", pch = 20, cex = 0.5)
      },
      track.height = 0.1
    )

    circlize::circos.track(
      factors = chrom_data$seqnames,
      x = chrom_data$start,
      y = chrom_data$start,
      panel.fun = function(x, y) {
        circos.rect(xleft = x, xright = x + 1, ybottom = 0, ytop = 1, col = "red", border = NA)
      },
      track.height = 0.1,
      bg.col = "gray90"   # Light background for better visibility
    )

    title(title)

  } else {
    stop("Invalid plot type. Choose from 'karyogram', 'density', 'lollipop', or 'circos'.")
  }
}

# [END] written by Zhenghao Xiao
