
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LegitXMut

<!-- badges: start -->
<!-- badges: end -->

A User-Friendly R Package for Nanopore Sequencing Mutation Analysis.
Offering A Efficient, Accessible Workflow from FASTQ File Alignment to
Visualizations.

# Description

`LegitXMut` is an R package designed for the analysis of genetic
variants, specifically designed for workflows that handle Nanopore
sequencing data. The workflow of the package involves aligning
sequencing reads, updating chromosome names in VCF files, and
visualizing mutation data, which enhance the efficiency and
accessibility of mutation analysis in bioinformatics pipelines.
`LegitXMut` integrates with widely used bioinformatics tools like
`Rsubread`, `VariantAnnotation`, `GenomicAlignments`, and others, to
help users to perform end-to-end genomic variant analysis within
`ONE PACKAGE` and `THREE FUNCTIONS`. Users have to change the
`file paths` to `local paths` for `tests` and `usage`.

This package fills a gap in bioinformatics workflows by allowing direct
handling of `FASTQ`, and `VCF` files for visualization of mutations,
which can aid in exploratory data analysis and the interpretation of
mutation patterns across various chromosomes. Unlike many existing
tools, `LegitXMut` is uniquely focused on a complete solution from
alignment to visualization in `R language`, offering functions for
generating compelling graphical outputs such as `heatmaps`,
`rainfall plots`, and `manhattan plots` to illustrate mutation
distributions and result confidences. The `LegitXMut` package was
developed using `R version 4.4.1 (June, 2024)`, Platform:
`x86_64-w64-mingw32/x64 (64-bit)`, running under `Windows 11`.

# Installation

``` r
# Install devtools if needed
install.packages("devtools")
library(devtools)

# Install LegitXMut from GitHub with vignettes
devtools::install_github("ZhenghaoXiao/LegitXMut", build_vignettes = TRUE)
library(LegitXMut)

# Rtools is a mandatory requirement for this package
# If errors occur mentioning any dependencies was failed to install
# Please manually install them on Bioconductor (Likely to be VariantAnnotation and GenomicAlignments)
# For Example:
#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("GenomicAlignments")
```

To run the Shiny app:

``` r
LegitXMut::runLegitXMut()
```

# Overview

The overview of the package:

``` r
ls("package:LegitXMut")
data(package = "LegitXMut") # empty
# Demo data of this package was downloaded from NCBI, as listed in the reference
# Users can access them in inst/exdata/, as mentioned in the function examples and vignettes
browseVignettes("LegitXMut")
```

`LegitXMut` provides the following functions for now: -
`alignment_FASTQ`:This function creates `BAM`and `VCF` file by aligning
a `FASTQ` file to a provided reference genome in `FASTA` format and
calls variant. It allows the setting of maximal insertions and
deletions(indel) mutations and mismatches criteria for quality control
and optimization. The result `.VCF` file should be input into
`update_vcf` and `plot_vcf_mutation_data` later on.

\-`update_vcf`:This tool updates the chromosome names in a `VCF` file by
extracting chromosomal mappings from a reference genome in `FASTA`
format. It ensures the workflows that demand multiple chromosomes and
whole genome analysis by producing a modified `VCF` file with
standardized chromosome names.

\-`plot_vcf_mutation_data`:This function takes the VCF file generated
from update_vcf or users can input external files and use them to create
different mutation visualizations. The three graphic types that users
can select from—“heatmap,” “manhattan,” and “rainfall”—each offers
distinct viewpoints on mutation data, including mutation density,
confidences of the result, and inter-variant distances across
chromosomes.

![](inst/extdata/workflow.png)

# Contributions

The author of `LegitXMut` is `Zhenghao Xiao`. The package provides three
functions to analyze genetic variation by allowing users to align
sequence reads to reference genome, standardize chromosome names, and
visualize mutation patterns in one package.

`alignment_FASTQ()`: `Rsubread` and `GenomicAlignments` was employed to
align a FASTQ file to a reference genome, generating `BAM` and `VCF`
files while allowing quality control through configurable parameters for
`indels` and `mismatches`. OpenAI was used for debugging suggestions and
potential function suggestions only, and the function was fully
constructed and written by the author.

`update_vcf()`: This function utilizes `VariantAnnotation`, `stringr`
and `GenomeInfoDb` to update chromosome names in VCF files according to
a reference genome `FASTA` file, the updated file will used for
visualizations. OpenAI was used for debugging suggestions and potential
function suggestions only, and the function was fully constructed and
written by the author.

`plot_vcf_mutation_data()`: This function requires `VariantAnnotation`,
`GenomicRanges`, `ggplot2`, `ComplexHeatmap`, `SummarizedExperiment`,
`dplyr`, `stats`, and `circlize`. Using a variety of plot patterns
(“heatmap,” “manhattan,” and “rainfall”), this tool visualizes mutation
data from VCF files to help with understanding chromosome-to-chromosome
inter-variant distances, mutation frequency density, and
distribution.OpenAI was used for debugging suggestions and potential
function suggestions only, and the function was fully constructed and
written by the author.

# References

1.BioRender. (2024). Image created by Zhenghao Xiao. Retrieved November
5, 2024, from <https://app.biorender.com/>

2.Chang, W., J. Cheng, J. Allaire, C. Sievert, B. Schloerke, Y. Xie, J.
Allen, J. McPherson, A. Dipert, B. Borges (2024). shiny: Web Application
Framework for R. R package version 1.9.1,
<https://CRAN.R-project.org/package=shiny>

3.Gu, Z. et al. “circlize Implements and Enhances Circular Visualization
in R.” Bioinformatics, vol. 30, no. 19, 2014, pp. 2811–2812.
<doi:10.1093/bioinformatics/btu393>.

4.Gu, Zuguang, Roland Eils, and Matthias Schlesner. “Complex Heatmaps
Reveal Patterns and Correlations in Multidimensional Genomic Data.”
Bioinformatics, vol. 32, no. 18, 2016, pp. 2847–2849.
<https://doi.org/10.1093/bioinformatics/btw313>.

5.Lawrence, Michael, et al. “Software for Computing and Annotating
Genomic Ranges.” PLoS Computational Biology, vol. 9, no. 8, 2013,
e1003118. <https://doi.org/10.1371/journal.pcbi.1003118>.

6.Liao, Y., Gordon K. Smyth, and Wei Shi. “The R Package Rsubread is
Easier, Faster, Cheaper and Better for Alignment and Quantification of
RNA Sequencing Reads.” Nucleic Acids Research, vol. 47, no. 8, 2019,
e47. <doi:10.1093/nar/gkz114>. Available at:
<https://academic.oup.com/nar/article/47/8/e47/5371636>.

7.Morgan, M., P. Aboyoun, R. Gentleman, M. Lawrence, and H. Pages.
“GenomicAlignments: Efficient Alignments Processing in R for NGS Data.”
Bioconductor, 2019, <doi:10.18129/B9.bioc.GenomicAlignments>. Available
at:
<https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html>.

8.Morgan, Martin, Vincent Obenchain, James Hester, and Hervé Pagès.
SummarizedExperiment: Summarized Experiment Container. R package version
1.28.0, Bioconductor, 2022,
<https://bioconductor.org/packages/SummarizedExperiment>.

9.National Center for Biotechnology Information (NCBI). Saccharomyces
cerevisiae S288C Genome Assembly (GCF_000146045.2). NCBI Datasets,
<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/>.
Accessed 4 Nov. 2024.

10.National Center for Biotechnology Information (NCBI). Sequence Read
Archive (SRA) Run: ERR12205202. NCBI SRA,
<https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR12205202&display=download>.
Accessed 4 Nov. 2024.

11.National Center for Biotechnology Information (NCBI). Sequence Read
Archive (SRA) Run: SRR29917898. NCBI SRA,
<https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29917898&display=download>.
Accessed 4 Nov. 2024.

12.OpenAI. ChatGPT: Assistance with R Function Development for
Bioinformatics Applications, “LegitxMut”. <https://chat.openai.com>.
Accessed 5 Nov. 2024.

13.Pagès, Hervé, Patrick Aboyoun, S. DebRoy, and Michael Lawrence.
GenomeInfoDb: Utilities for Manipulating Chromosome Names, Including
Modifying the Names to Follow a Particular Convention. R package version
1.36.0, Bioconductor, 2023,
<https://bioconductor.org/packages/GenomeInfoDb>.

14.R Core Team. R: A Language and Environment for Statistical Computing.
R Foundation for Statistical Computing, 2024, Vienna, Austria,
<https://www.R-project.org/>.

15.Silva, Anjali. TestingPackage. GitHub,
<https://github.com/anjalisilva/TestingPackage>. Accessed 5 Nov. 2024.

16.Wickham, Hadley. ggplot2: Elegant Graphics for Data Analysis.
Springer-Verlag, 2016, <https://ggplot2.tidyverse.org>.

17.Wickham, Hadley, et al. dplyr: A Grammar of Data Manipulation. R
package version 1.1.2, 2023, <https://CRAN.R-project.org/package=dplyr>.

18.Wickham, H. stringr: Simple, Consistent Wrappers for Common String
Operations. R package version 1.4.0, 2019,
<https://CRAN.R-project.org/package=stringr>.

# Acknowledgements

This package was developed as part of an assessment for 2024 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `LegitXMut`welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues.

# Notes

The package may take a while to download, because it requires multiple
dependencies. The workflow of this version of the package only supports
the display of insertions and deletions(indel) mutation patterns due to
the current implementation of the alignment function. However, the
Manhattan plot visualization in this version can display single
nucleotide polymorphisms(SNPs) when provided with appropriate input
data.

# Demonstration

The following figures were generating using the data set SRR29917898
![](inst/extdata/heatmap.png)

![](inst/extdata/manhattan.png)

![](inst/extdata/rainfall.png)

The package tree structure is provided below.

    - LegitXMut
      |- LegitXMut.Rproj
      |- DESCRIPTION
      |- NAMESPACE
      |- LICENSE
      |- README.md
      |- inst
        CITATION
        |- extdata
          |- ERR12205202.fastq
          |- aligned_output.bam
          |- aligned_output.bam.indel.vcf
          |- ERR12205202.fastq
          |- heatmap.png
          |- manhattan.png
          |- rainfall.png
          |- updated.vcf
          |- workflow.png
          |- yeast.fna
        |- shiny-scripts
          |- app.R
      |- man
        |- alignment_FASTQ.Rd
        |- update_vcf.Rd
        |- plot_vcf_mutation_data.Rd
        |- runLegitXMut.Rd
      |- R
        |- alignment_FASTQ.R
        |- update_vcf.R
        |- plot_vcf_mutation_data.R
        |- runLegitXMut.R
      |- vignettes
        |- heatmap.png
        |- manhattan.png
        |- rainfall.png
        |- Introduction_Of_LegitXMut.Rmd
      |- tests
        |- testthat.R
        |- testthat
          |- test-alignment_FASTQ.R
          |- test-update_vcf.R
          |- test-plot_vcf_mutation_data.R
