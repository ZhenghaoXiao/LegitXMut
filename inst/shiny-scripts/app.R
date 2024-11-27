# This example is adapted from
# Chang, W., J. Cheng, J. Allaire, C. Sievert, B. Schloerke, Y. Xie, J. Allen, J. McPherson, A. Dipert, B. Borges (2024).
# shiny: Web Application Framework for R. R package version 1.9.1, https://CRAN.R-project.org/package=shiny
# Silva, Anjali. TestingPackage. GitHub, https://github.com/anjalisilva/TestingPackage. Accessed 5 Nov. 2024.

library(shiny)
library(LegitXMut)

# Increase maximum file upload size
options(shiny.maxRequestSize = 100000 * 1024^2)

# Define the UI
ui <- fluidPage(
  titlePanel("Genomic Data Processing and Visualization in LegitXMut"),

  sidebarLayout(
    sidebarPanel(
      h4("Step 1: Align FASTQ to Reference"),
      fileInput("fastqFile", "Upload FASTQ File:", accept = c(".fastq")),
      fileInput("referenceFile", "Upload Reference Genome (FASTA):", accept = c(".fasta", ".fna")),
      numericInput("indels", "Maximum Indels Allowed:", value = 10, min = 0, step = 1),
      numericInput("maxMismatches", "Maximum Mismatches Allowed:", value = 1000, min = 0, step = 1),
      textInput("outputBAM", "Output BAM File PATH:", value = "output.bam"),
      actionButton("alignBtn", "Align FASTQ"),
      hr(),

      h4("Step 2: Update VCF with Reference"),
      fileInput("vcfFile", "Upload VCF File:", accept = c(".vcf")),
      textInput("outputVCF", "Output VCF File PATH:", value = "updated.vcf"),
      actionButton("updateBtn", "Update VCF"),
      hr(),

      h4("Step 3: Visualize Mutation Data"),
      selectInput("plotType", "Choose Plot Type:",
                  choices = c("Heatmap" = "heatmap",
                              "Rainfall" = "rainfall",
                              "Manhattan" = "manhattan")),
      sliderInput("fontSize", "Font Size:", min = 8, max = 20, value = 12),
      textInput("plotTitle", "Plot Title:", value = "Mutation Data Plot"),
      textInput("xlab", "X-axis Label:", value = "Genomic Position"),
      textInput("ylab", "Y-axis Label:", value = "Mutation Frequency"),
      selectInput("legendPosition", "Legend Position:",
                  choices = c("top", "bottom", "left", "right"), selected = "top"),
      textInput("colorScheme", "Custom Color Scheme (Comma-separated):",
                placeholder = "e.g., red,blue,green"),
      actionButton("plotBtn", "Generate Plot")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Alignment", verbatimTextOutput("alignmentLog")),
        tabPanel("VCF Update", verbatimTextOutput("vcfLog")),
        tabPanel("Visualization", plotOutput("mutationPlot"))
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
      updateTextInput(session, "outputBAM", value = file.path(fastqDir, "output.bam"))
    }
  })

  # Dynamically set the output VCF path based on the uploaded VCF file
  observeEvent(input$vcfFile, {
    if (!is.null(input$vcfFile)) {
      vcfDir <- dirname(input$vcfFile$datapath)
      updateTextInput(session, "outputVCF", value = file.path(vcfDir, "updated.vcf"))
    }
  })

  # Alignment Logic
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

  # VCF Update Logic
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

  # Visualization Logic
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
