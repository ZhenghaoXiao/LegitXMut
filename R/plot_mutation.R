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
#' # Example 2: Generate a heatmap of variants' frequency across chromosomes
#' plot_vcf_mutation_data(
#'   vcfPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf",
#'   plotType = "density"
#' )
#'
#' # Example 3: Generate a manhattan plot for mutations in the TP53 gene
#' plot_vcf_mutation_data(
#'   vcfPath = "C:/Users/rjay1/Desktop/BCB410/LegitXMut/inst/extdata/updated.vcf",
#'   plotType = "manhattan",
#'   title = "Manhattan Plot of Mutations in TP53 Gene on chromosome 17"
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
#' @import circlize
#' @import GenomicRanges
#' @import ggplot2
#' @import ComplexHeatmap
plot_vcf_mutation_data <- function(vcfPath, plotType = "karyogram", species = "hg38", title = "Plot") {
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

  if (plotType == "heatmap") {
    chrom_counts <- as.data.frame(table(seqnames(gr)))
    colnames(chrom_counts) <- c("Chromosome", "Mutation_Count")

    mutation_matrix <- as.matrix(chrom_counts["Mutation_Count"])
    rownames(mutation_matrix) <- chrom_counts$Chromosome

    if (length(unique(mutation_matrix)) == 1) {
      warning("Only one unique mutation count present. Adjusting color scale for single value.")
      mutation_matrix <- cbind(mutation_matrix, mutation_matrix + 0.01)
      colnames(mutation_matrix) <- c("Mutation_Count", "Adjusted_Value")
    }

    annotation_row <- data.frame(
      Group = ifelse(mutation_matrix[, 1] > median(mutation_matrix[, 1]), "High Mutation", "Low Mutation")
    )
    rownames(annotation_row) <- rownames(mutation_matrix)

    annotation_colors <- list(
      Group = c("High Mutation" = "darkgreen", "Low Mutation" = "purple")
    )

    col_fun <- circlize::colorRamp2(
      breaks = c(min(mutation_matrix), max(mutation_matrix)),
      colors = c("white", "red")
    )

    row_ha <- ComplexHeatmap::rowAnnotation(
      Group = annotation_row$Group,
      col = annotation_colors
    )

    ComplexHeatmap::Heatmap(
      mutation_matrix,
      name = "Mutation Count",
      col = col_fun,
      cluster_rows = nrow(mutation_matrix) > 1,
      cluster_columns = FALSE,
      row_names_side = "left",
      column_names_side = "top",
      row_title = "Chromosomes",
      column_title = title,
      show_row_dend = nrow(mutation_matrix) > 1,
      row_names_gp = grid::gpar(fontsize = 10),
      column_names_gp = grid::gpar(fontsize = 10, rot = 45),
      left_annotation = row_ha,
      heatmap_legend_param = list(
        title = "Mutation Count",
        at = c(min(mutation_matrix), max(mutation_matrix)),
        labels = c("Low", "High")
      )
    )
  }
  else if (plotType == "manhattan") {

    vcf_data <- VariantAnnotation::readVcf(vcfPath, genome = "hg38")

    info_data <- info(vcf_data)

    if ("SR" %in% colnames(info_data)) {
      SR_values <- as.numeric(info_data$SR)
    } else {
      stop("The 'SR' field is not present in the INFO column. Ensure the VCF has SR values.")
    }

    SR_values[is.na(SR_values)] <- 1

    if (length(SR_values) != length(rowRanges(vcf_data))) {
      stop("Mismatch in SR values and number of VCF records.")
    }

    chrom_data <- data.frame(
      Chromosome = as.character(seqnames(rowRanges(vcf_data))),
      Position = start(rowRanges(vcf_data)),
      SR = SR_values
    )

    p <- ggplot2::ggplot(chrom_data, aes(x = Position, y = -log10(SR), color = Chromosome)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::facet_wrap(~ Chromosome, scales = "free_x") +
      ggplot2::labs(title = title, x = "Position", y = "-log10") +
      ggplot2::theme_minimal()
    print(p)
  }
  else if (plotType == "circos") {
    chrom_data <- as.data.frame(gr)

    circlize::circos.clear()

    circlize::circos.initializeWithIdeogram(
      species = species,
      plotType = c("ideogram", "labels")
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
      bg.col = "gray90"
    )

    title(title)

  } else {
    stop("Invalid plot type. Choose from 'karyogram', 'heatmap', 'manhattan', or 'circos'.")
  }
}

# [END] written by Zhenghao Xiao
