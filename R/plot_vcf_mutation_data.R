#' Visualize Mutation Data from VCF File in Various Plot Styles
#'
#' This function takes a VCF file and generates different mutation visualization plots,
#' including heatmap, rainfall plot, and manhattan plot. It provides
#' genome-specific and, across chromosomes mutation visualizations directly from the VCF file.
#'
#' @param vcfPath A character string specifying the path to the updated VCF file.
#' @param plotType A character string specifying the type of plot to generate.
#'    Options are "heatmap" (heatmap of variants across chromosomes), "manhattan" (variant distribution),
#'    or "rainfall" (distance between variants).
#' @param title A character string specifying the title of the plot.
#' @param color_scheme A named vector of colors for different mutation types.
#'    Colors are fixed for "manhattan", but users can customize for "heatmap" or "rainfall".
#' @param alpha A numeric value between 0 and 1 to set the transparency of points in plots (default is 0.6).
#' @param font_size A numeric value for the font size of plot labels and titles (default is 10).
#' @param xlab A character string for the x-axis label (default is "Position").
#' @param ylab A character string for the y-axis label, which varies by plot type.
#' @param legend_position A character string specifying the legend position ("top", "bottom", "left", "right").
#'
#' @return A plot visualizing mutation information based on the selected plot type and parameters
#'
#' @examples
#' # Example 1: Generate a heatmap of variants' frequency across chromosomes
#' #xlab, ylab, and legend position are fixed for heatmap
#' vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")
#' plot_vcf_mutation_data(
#'   vcfPath = vcfPath,
#'   plotType = "heatmap",
#'   title = "Heatmap of Mutations Density",
#'   font_size = 12,
#' )
#'
#' # Example 2: Generate a rainfall plot of mutations across chromosomes with customized colors
#' # Demo data only have insertion-deletion variants of more than one nucleotides changes
#' # Only indel will be visualized because the this alignment method only output indel into VCF
#' vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")
#' plot_vcf_mutation_data(
#'   vcfPath = vcfPath,
#'   plotType = "rainfall",
#'   title = "Rainfall Plot of Mutations",
#'   color_scheme = c("C>A" = "red", "C>G" = "orange", "C>T" = "green",
#'                    "T>A" = "yellow", "T>C" = "blue", "T>G" = "purple",
#'                    "indel" = "grey"),
#'   alpha = 0.7,
#'   font_size = 10,
#'   xlab = "Genomic Position",
#'   ylab = "Inter-Variant Distance (log10)",
#'   legend_position = "right"
#' )
#'
#' # Example 3: Generate a manhattan plot for mutations with modified y-axis label
#' # Manhattan plots do not allow color modifications
#' # Since the potentially high number of chromosomes displayed
#' vcfPath <- system.file("extdata", "updated.vcf", package = "LegitXMut")
#' plot_vcf_mutation_data(
#'   vcfPath = vcfPath,
#'   plotType = "manhattan",
#'   title = "Manhattan Plot of Variants",
#'   ylab = "-log10(Supporting Reads)",
#'   xlab = "Genomic Position",
#'   font_size = 11,
#'   alpha = 0.6,
#'   legend_position = "top"
#' )
#'
#' @references
#' Gu, Z., Eils, R., and Schlesner, M. (2016).
#' "Complex heatmaps reveal patterns and correlations in multidimensional genomic data."
#' Bioinformatics, 32(18), 2847-2849. doi:10.1093/bioinformatics/btw313.
#'
#' Gu, Z., Gu, L., Eils, R., Schlesner, M., and Brors, B. (2014).
#' "circlize implements and enhances circular visualization in R."
#' Bioinformatics, 30(19), 2811-2812. doi:10.1093/bioinformatics/btu393.
#'
#' Lawrence, Michael, et al. "Software for computing and annotating genomic ranges."
#' PLoS Computational Biology, vol. 9, no. 8, 2013, e1003118. https://doi.org/10.1371/journal.pcbi.1003118.
#'
#' Morgan M, Pages H, Obenchain V, and Hayden N (2020).
#' "VariantAnnotation:annotation of genetic variants."
#' Bioconductor. doi:10.18129/B9.bioc.VariantAnnotation. Available at: https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html.
#'
#' Morgan, Martin, Vincent Obenchain, James Hester, and Hervé Pagès.
#' SummarizedExperiment: Summarized Experiment Container. R package version 1.28.0, Bioconductor, 2022, https://bioconductor.org/packages/SummarizedExperiment.
#'
#' National Center for Biotechnology Information (NCBI). Saccharomyces cerevisiae S288C Genome Assembly (GCF_000146045.2).
#' NCBI Datasets, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/. Accessed 4 Nov. 2024.
#'
#' National Center for Biotechnology Information (NCBI). Sequence Read Archive (SRA) Run: ERR12205202.
#' NCBI SRA, https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR12205202&display=download. Accessed 4 Nov. 2024.
#'
#' OpenAI. ChatGPT: Assistance with R function development for bioinformatics applications,
#' "Assignment". https://chat.openai.com (2023, accessed 5 November 2024 for debugging).
#'
#' R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria, 2024, https://www.R-project.org/.
#'
#' Wickham, Hadley, et al. dplyr: A Grammar of Data Manipulation.
#' R package version 1.1.2, 2023, https://CRAN.R-project.org/package=dplyr.
#'
#' Wickham, Hadley. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, 2016. https://ggplot2.tidyverse.org.
#'
#' @export
#' @importFrom VariantAnnotation readVcf ref alt info
#' @importFrom circlize colorRamp2
#' @importFrom GenomicRanges seqnames ranges start GRanges mcols
#' @importFrom ggplot2 ggplot aes geom_point scale_y_log10 scale_color_manual labs facet_wrap element_text theme_minimal theme guides guide_legend
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom dplyr lag
#' @importFrom stats median
plot_vcf_mutation_data <- function(vcfPath, plotType = "manhattan", title = "Plot",
                                   color_scheme = NULL, alpha = 0.6, font_size = 10,
                                   xlab = "Position", ylab = "Value", legend_position = "top") {

  # --- INPUT VALIDATION ---
  if (!file.exists(vcfPath)) {
    stop("The specified VCF file does not exist: ", vcfPath)
  }

  # --- LOAD VCF FILE ---
  # Read the VCF file to extract mutation data
  message("Reading VCF file: ", vcfPath)
  vcf <- VariantAnnotation::readVcf(vcfPath)
  row_ranges <- SummarizedExperiment::rowRanges(vcf)
  seqnames <- as.character(GenomicRanges::seqnames(row_ranges))
  ranges <- GenomicRanges::ranges(row_ranges)
  positions <- GenomicRanges::start(row_ranges)

  # Ensure seqnames match the number of ranges for consistent data
  if (length(seqnames) == 1 && length(ranges) > 1) {
    seqnames <- rep(seqnames, length(ranges))
  } else if (length(seqnames) != length(ranges)) {
    stop("Length of 'seqnames' does not match the number of ranges.")
  }

  # Validate alpha
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("Alpha must be a numeric value between 0 and 1.")
  }

  # --- CREATE GRanges OBJECT ---
  # Store information
  gr <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = ranges,
    strand = "*",
    REF = GenomicRanges::mcols(row_ranges)$REF,
    ALT = GenomicRanges::mcols(row_ranges)$ALT
  )

  # --- HEATMAP ---
  if (plotType == "heatmap") {
    # Aggregate mutation counts per chromosome
    chrom_counts <- as.data.frame(table(seqnames(gr)))
    colnames(chrom_counts) <- c("Chromosome", "Mutation_Count")
    chrom_counts <- chrom_counts[order(-chrom_counts$Mutation_Count), ]

    mutation_matrix <- as.matrix(chrom_counts["Mutation_Count"])# Create a matrix of chrom_counts for heatmap visualization
    rownames(mutation_matrix) <- chrom_counts$Chromosome

    # Define color mapping
    col_fun <- if (is.null(color_scheme)) {
      circlize::colorRamp2(breaks = c(min(mutation_matrix), max(mutation_matrix)),
                           colors = c("white", "darkred"))
    } else {
      color_scheme
    }

    annotation_row <- data.frame(Group = ifelse(mutation_matrix[, 1] > stats::median(mutation_matrix[, 1]),# Create annotations for rows
                                                "High Mutation", "Low Mutation"))
    rownames(annotation_row) <- rownames(mutation_matrix)

    annotation_colors <- list(Group = c("High Mutation" = "darkgreen", "Low Mutation" = "purple"))

    # Plot
    message("Generating heatmap...")
    ComplexHeatmap::Heatmap(mutation_matrix, name = "Mutation Count", col = col_fun,
                            cluster_rows = FALSE, cluster_columns = FALSE,
                            show_row_names = TRUE, show_column_names = FALSE,
                            row_title = "Chromosomes (High to Low Density)",
                            column_title = title, row_names_gp = grid::gpar(fontsize = font_size),
                            left_annotation = ComplexHeatmap::rowAnnotation(Group = annotation_row$Group, col = annotation_colors),
                            heatmap_legend_param = list(title = "Mutation Count", at = c(min(mutation_matrix), max(mutation_matrix)),
                            labels = c("Low", "High"))
    )
  }
  else if (plotType == "manhattan") {
    # --- INPUT VALIDATION ---
    info_data <- VariantAnnotation::info(vcf)
    if (!"SR" %in% colnames(info_data)) {
      stop("The 'SR' field is not present in the INFO column of the VCF.")
    }
    # Convert SR values to numeric and handle missing data
    SR_values <- as.numeric(info_data$SR)# Supporting Reads(SR) values are used in -log10(SR) for visualization
    SR_values[is.na(SR_values)] <- 1  # Replace missing values with 1 to avoid errors in calculation and visualization

    # Ensure SR values match the number of VCF records
    if (length(SR_values) != length(SummarizedExperiment::rowRanges(vcf))) {
      stop("Mismatch in SR values and number of VCF records.")
    }

    # --- DATA PREPARATION ---
    # Create a data frame of SR values
    chrom_data <- data.frame(
      Chromosome = as.character(seqnames(SummarizedExperiment::rowRanges(vcf))),
      Position = start(SummarizedExperiment::rowRanges(vcf)),
      SR = SR_values
    )

    # --- PLOT GENERATION: MANHATTAN ---
    # Create a Manhattan plot
    p <- ggplot2::ggplot(chrom_data, ggplot2::aes(x = Position, y = -log10(SR), color = Chromosome)) +
      ggplot2::geom_point(alpha = alpha) +
      ggplot2::facet_wrap(~ Chromosome, scales = "free_x") +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = legend_position, text = element_text(size = font_size))# Place the legend and set front size
    print(p)# Display the plot
  }
  # --- PLOT GENERATION: RAINFALL ---
  else if (plotType == "rainfall") {
    # --- INPUT VALIDATION ---
    # Extract REF and ALT alleles from the VCF
    ref_alleles <- as.character(VariantAnnotation::ref(vcf))
    alt_alleles <- as.character(unlist(VariantAnnotation::alt(vcf)))

    # Ensure REF and ALT lengths match the data
    if (length(ref_alleles) != length(alt_alleles)) {
      stop("Mismatch in REF and ALT allele lengths.")
    }

    # --- DATA PREPARATION ---
    # Create the data frame
    variant_data <- data.frame(
      Chromosome = as.character(seqnames(row_ranges)),
      Position = positions,
      REF = ref_alleles,
      ALT = alt_alleles)

    variant_data$Mutation_Type <- paste(variant_data$REF, ">", variant_data$ALT, sep = "")# Calculate mutation type
    variant_data <- variant_data[order(variant_data$Chromosome, variant_data$Position), ]# Sort data by chromosome and position
    variant_data$Distance <- c(NA, diff(variant_data$Position))# Calculate inter-variant distances
    variant_data$Distance[variant_data$Chromosome != dplyr::lag(variant_data$Chromosome)] <- NA
    variant_data <- variant_data[!is.na(variant_data$Distance) & variant_data$Distance > 0, ]# Remove invalid or NA distances

    # Define default colors
    mutation_colors <- if (is.null(color_scheme)) {
      c("C>A" = "red", "C>G" = "orange", "C>T" = "green", "T>A" = "yellow",
        "T>C" = "blue", "T>G" = "purple", "indel" = "grey")
    } else {
      color_scheme
    }
    # Handle indels by reclassifying mutation types
    variant_data$Mutation_Type <- ifelse(
      nchar(variant_data$REF) > 1 | nchar(variant_data$ALT) > 1,
      "indel",
      variant_data$Mutation_Type
    )
    # Create a Rainfall plot
    p <- ggplot2::ggplot(variant_data, aes(x = Position, y = Distance, color = Mutation_Type)) +
      ggplot2::geom_point(alpha = alpha) +
      ggplot2::scale_y_log10() + # Log scale for distance
      ggplot2::scale_color_manual(values = mutation_colors) + # Assign colors
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::facet_wrap(~ Chromosome, scales = "free_x") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = legend_position,
                     text = ggplot2::element_text(size = font_size)) +
      ggplot2::guides(color = ggplot2::guide_legend(title = "Mutation Type"))# Add legend title

    print(p)
  } else {
    stop("Invalid plot type. Choose from 'heatmap', 'manhattan', or 'rainfall'.")
  }
}

# [END] written by Zhenghao Xiao
