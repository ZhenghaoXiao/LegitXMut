# This example is adapted from
# Chang, W., J. Cheng, J. Allaire, C. Sievert, B. Schloerke, Y. Xie, J. Allen, J. McPherson, A. Dipert, B. Borges (2024).
# shiny: Web Application Framework for R. R package version 1.9.1, https://CRAN.R-project.org/package=shiny

# Silva, Anjali. TestingPackage. GitHub, https://github.com/anjalisilva/TestingPackage. Accessed 5 Nov. 2024.

library(shiny)
library(LegitXMut)

# Set maximum file upload size to handle large genomic data
options(shiny.maxRequestSize = 100000 * 1024^2)

# Define the User Interface (UI)
ui <- fluidPage(
  # App title
  titlePanel("Alignment and Visualization in LegitXMut"),

  # Introduction and instructions section
  fluidRow(
    column(12,
           tags$p("Welcome to the LegitXMut Shiny Application! This app provides a streamlined workflow for LegitXMut,
              focusing on alignment, VCF updates, and mutation patterns visualizations."),
           tags$h3("Description:"),
           tags$p("LegitXMut is an R package tailored for mutation analyze for Nanopore sequencing data. With this app,
              you can align FASTQ files to a reference genome, update VCF files, and create visually engaging plots
              to explore mutation patterns. Please fully read through the following instructions for the best user experience"),
           tags$h3("Instructions:"),
           tags$ol(
             tags$li("In 'Step 1', upload your sequencing data and reference genome to generate BAM and VCF files."),
             tags$li("In 'Step 2', update chromosome mappings in your VCF file."),
             tags$li("In 'Step 3', customize and generate mutation data visualizations from the uploaded VCF file.")
           ),
           tags$p("Navigate through the tabs on the top right to view logs and plots. Customize settings in each step for better control
              of the analysis."),
           tags$h3("Demo Data:"),
           tags$p("Demo files are available in the package's 'extdata' directory. Use the following R code to access them:"),
           tags$pre(
             "demo_fastq <- system.file(\"extdata\", \"ERR12205202.fastq\", package = \"LegitXMut\")\n",
             "demo_fasta <- system.file(\"extdata\", \"yeast.fna\", package = \"LegitXMut\")\n",
             "demo_vcf <- system.file(\"extdata\", \"aligned_output.bam.indel.vcf\", package = \"LegitXMut\")\n",
             "print(demo_fastq)\n",
             "print(demo_fasta)\n",
             "print(demo_vcf)"
           ),
           tags$h4("For visualizations, xlab, ylab, and legend position are fixed for heatmap;
                  Demo data only have insertion-deletion variants of more than one nucleotides
                  changes and only indel will be visualized at the moment because the this alignment method only output indel into VCF;
                  Remember to reupload the updated.vcf to the shiny app before visualization;
                  Manhattan plots do not allow color modifications since the potentially high number of chromosomes displayed")
  ),

  # Sidebar layout with inputs and main outputs
  sidebarLayout(
    sidebarPanel(
      # Alignment inputs
      h4("Step 1: Align FASTQ to Reference"),
      fileInput("fastqFile", "Upload FASTQ File:", accept = c(".fastq")),# Input for FASTQ file
      fileInput("referenceFile", "Upload Reference Genome (FASTA):", accept = c(".fasta", ".fna")),# Input for reference genome
      numericInput("indels", "Maximum Indels Allowed:", value = 10, min = 0, step = 1), # indels allowed
      numericInput("maxMismatches", "Maximum Mismatches Allowed:", value = 1000, min = 0, step = 1), # Number of mismatches allowed
      textInput("outputBAM", "Output BAM File PATH:", value = "output.bam"), # Output BAM file path
      actionButton("alignBtn", "Align FASTQ"), # Button to trigger alignment
      hr(),

      # VCF update inputs
      h4("Step 2: Update VCF with Reference"),
      fileInput("vcfFile", "Upload VCF File:", accept = c(".vcf")), # Input for VCF file
      textInput("outputVCF", "Output VCF File PATH:", value = "updated.vcf"), # Output VCF file path
      actionButton("updateBtn", "Update VCF"), # Button to trigger VCF update
      hr(),

      # Visualization inputs
      h4("Step 3: Visualize Mutation Data"),
      selectInput("plotType", "Choose Plot Type:", # Select plot type
                  choices = c("Heatmap" = "heatmap",
                              "Rainfall" = "rainfall",
                              "Manhattan" = "manhattan")),
      sliderInput("fontSize", "Font Size:", min = 8, max = 20, value = 12), # Font size for plots
      textInput("plotTitle", "Plot Title:", value = "Mutation Data Plot"), # Plot title
      textInput("xlab", "X-axis Label:", value = "Genomic Position"), # X-axis label
      textInput("ylab", "Y-axis Label:", value = "Mutation Frequency"), # Y-axis label
      selectInput("legendPosition", "Legend Position:", # Legend position
                  choices = c("top", "bottom", "left", "right"), selected = "top"),
      textInput("colorScheme", "Custom Color Scheme (Comma-separated):",
                placeholder = "e.g., red,blue,green"), # Custom color scheme
      actionButton("plotBtn", "Generate Plot") # Button to trigger plot generation
    ),

    mainPanel(
      # Tabbed outputs for alignment, VCF update, and visualizations
      tabsetPanel(
        tabPanel("Alignment",
                 tags$h5("Alignment Logs"),
                 tags$p("Logs for the FASTQ alignment process will be displayed here."),
                 verbatimTextOutput("alignmentLog")),
        tabPanel("VCF Update",
                 tags$h5("VCF Update Logs"),
                 tags$p("Logs for the VCF update process will be displayed here."),
                 verbatimTextOutput("vcfLog")),
        tabPanel("Visualization",
                 tags$h5("Visualization Output"),
                 tags$p("The plot generated from the VCF file will appear below."),
                 plotOutput("mutationPlot"))
      )
    )
  )
)
)

# Define the server
server <- function(input, output, session) {
  # Dynamically set the output BAM path based on the uploaded FASTQ file
  observeEvent(input$fastqFile, {
    if (!is.null(input$fastqFile)) {
      fastqDir <- dirname(input$fastqFile$datapath)
      updateTextInput(session, "outputBAM", value = file.path(fastqDir, "output.bam")) # Update output BAM path
    }
  })

  # Dynamically set the output VCF path based on the uploaded VCF file
  observeEvent(input$vcfFile, {
    if (!is.null(input$vcfFile)) {
      vcfDir <- dirname(input$vcfFile$datapath)
      updateTextInput(session, "outputVCF", value = file.path(vcfDir, "updated.vcf")) # Update output VCF path
    }
  })

  # Alignment
  output$alignmentLog <- renderText({
    req(input$alignBtn)
    validate(
      need(input$fastqFile$datapath, "Please upload a FASTQ file."),
      need(input$referenceFile$datapath, "Please upload a reference genome.")
    )
    isolate({
      tryCatch({
        # Perform alignment
        LegitXMut::alignment_FASTQ(
          fastqPath = input$fastqFile$datapath,
          referencePath = input$referenceFile$datapath,
          outputBAM = input$outputBAM,
          indels = input$indels,
          maxMismatches = input$maxMismatches
        )

        # Message to notify about the BAM and VCF files
        paste("Alignment completed.",
              "\nBAM file saved at:", input$outputBAM,
              "\nThe VCF file has been output to the same directory as the BAM file and ends with '.indel.vcf'")
      }, error = function(e) {
        paste("Error during alignment:", e$message)
      })
    })
  })

  # VCF Update
  output$vcfLog <- renderText({
    req(input$updateBtn)
    validate(
      need(input$vcfFile$datapath, "Please upload a VCF file."),
      need(input$referenceFile$datapath, "Please upload a reference genome.")
    )
    isolate({
      tryCatch({
        LegitXMut::update_vcf(
          fastaPath = input$referenceFile$datapath,
          vcfPath = input$vcfFile$datapath,
          outputVcfPath = input$outputVCF
        )
        paste("VCF updated successfully. File saved at:", input$outputVCF)
      }, error = function(e) {
        paste("Error during VCF update:", e$message)
      })
    })
  })

  # Visualization
  output$mutationPlot <- renderPlot({
    req(input$plotBtn)
    validate(
      need(input$vcfFile$datapath, "Please upload a VCF file.")
    )
    isolate({
      tryCatch({
        # Handle custom color scheme
        color_scheme <- if (!is.null(input$colorScheme) && nchar(input$colorScheme) > 0) {
          strsplit(input$colorScheme, ",")[[1]]
        } else {
          NULL
        }

        LegitXMut::plot_vcf_mutation_data(
          vcfPath = input$vcfFile$datapath,
          plotType = input$plotType,
          title = input$plotTitle,
          color_scheme = color_scheme,
          font_size = input$fontSize,
          xlab = input$xlab,
          ylab = input$ylab,
          legend_position = input$legendPosition
        )
      }, error = function(e) {
        showNotification(paste("Error generating plot:", e$message), type = "error")
      })
    })
  })
}

# Run the Shiny App
shiny::shinyApp(ui = ui, server = server)
