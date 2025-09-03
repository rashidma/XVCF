
# Script to make a shiny app for VCFR and Maftool packages 
# optional for cancer genoume  
#  10/11/2024 
# 2ANNOVAR_output_for_Shiny_app.R
# ANNOVAR output (tabular format) + clinical data
#getwd()
#setwd("C:/Users/ghaid/Documents/training scripts/App-2/")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
 
BiocManager::install("maftools")
#install.packages("DT")
#install.packages("officer")
#install.packages("flextable")
#install.packages("shinyjs")


library(shiny)
library(vcfR)
library(reshape2)
library(ggplot2)
library(shinydashboard)
library(maftools)
library(tidyr)
library(readxl)
library(DT)
library(shinyjs)  # For enabling/disabling components dynamically

#library(R.utils)
#library(officer)
#library(flextable)
 
 

# Increase maximum file size limit to 200MB (adjust as needed)
options(shiny.maxRequestSize = 2000 * 1024^2)

# Define a reactiveVal to store VCF data
vcf_data <- reactiveVal(NULL)
ANNOVAR_data <- reactiveVal(NULL)
clindata <- reactiveVal(NULL)

# UI for the app
ui <- dashboardPage(
  dashboardHeader(title = "XVCF"),
  dashboardSidebar(
    fileInput("vcfFile", "Choose VCF File", multiple = TRUE, accept = ".vcf"),
     tags$hr(),
    sidebarMenu(
     # id = "tabsVCF",  # This ID is referenced in `updateTabItems`
      menuItem("Home", tabName = "home" , icon = icon("home"), startExpanded = TRUE,  #,icon = icon("chart-bar")
           id = "tabsVCF",  # This ID is referenced in `updateTabItems`
               
                   menuSubItem("VCF Summary ",
                           tabName = "home"),
               
               menuSubItem("Read Depth Plot as a Violin plot",
                           tabName = "page2"),
               menuSubItem("Read Depth as a Box Plot",
                           tabName = "page3"),
               menuSubItem("Heatmap of Read Depth Plot",
                           tabName = "page6"),
               menuSubItem("Bar Plot of Read Depth",
                           tabName = "page7"),
               menuSubItem("Genotype Plot", tabName = "page4"),
               menuSubItem("Genotype Quality Plot", tabName = "page5"),
               menuSubItem("Quality Metrics Overview", tabName = "page8"),
               menuSubItem("Histogram Summary of Quality Metrics", tabName = "page9"),
               menuSubItem("Allele Frequencies", tabName = "page10"),
               menuSubItem("Allele Frequencies For Samples", tabName = "page11")
             #  menuSubItem("Minor Allele Frequency", tabName = "page12")
               ),
      tags$hr(),
      
      ###
      
      useShinyjs(),  # Initialize shinyjs
      
      # Input for the number of samples to extract
      fluidRow(
        numericInput(
          inputId = "num_samples",
          label = "Enter The Number of Samples In Your ANNOVAR File To upload The File:",
          value = NULL,    # Start with no value
          min = 1,
          max = 200
        )
      ),
      
      # Placeholder for fileInput
      fluidRow(
        uiOutput("fileInputUI")
      ),
      
      # Add a message for better user experience
      fluidRow(
        textOutput("validationMessage")
      ),
      
     
      selectInput("AF_Threshold", 
                  label = "AF Enter a value for filtering on AF to enrich somatic mutation:", 
                  choices = c("None" = NA,  "0.1" = 0.1, "0.2" = 0.2,"0.3" = 0.3,"0.4" = 0.4, "0.5" = 0.5), 
                  selected = NA)
      ,
      fileInput("xlsxFile", "Choose Clinical Data File", multiple = TRUE, accept = ".Xlsx"),
      tags$hr(),
      id = "tabs",  # This ID is referenced in `updateTabItems`
      menuItem("Maftool Summary", tabName = "page13", icon =  icon("table")),
      menuItem("Maftool Summary plot", tabName = "page14"),
      menuItem("oncoplot for top ten mutated genes", tabName = "page15"),
      menuItem("plot titv summary", tabName = "page16"),
      menuItem("lollipop plot", tabName = "page17"),
      
      menuItem("rainfall plot", tabName = "page18"),
      
      menuItem("Compare with TCGA cohorts", tabName = "page19"),
      
      menuItem("Oncogenic Pathways plot", tabName = "page20"),
      tags$hr(),
      menuItem("Contact", tabName = "contact")
    ) 
  ),
  dashboardBody(
    
    tabItems(
      tabItem(tabName = "home",
              h3("Welcome to XVCF"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("An Extensible Data Visualization and Quality Assurance Platform for Integrated Variant Analysis. An interactive Rshiny tool designed to support the evaluation of mutation calls from sequencing data. The tool takes as input a single variant call format (VCF) file and enables researchers to explore the impacts of analytical choices on the mutant allele frequency spectrum, on mutational signatures and on annotated somatic variants in genes of interest. It allows variants that have failed caller filters to be re-examined in order to improve sensitivity or guide strategies for optimal sample / sequencing preparation. It is extensible allowing other algorithms to take advantage of its VCF preprocessing capabilities.")
                  ) , width = 12 ,background = "light-blue"  #"navy"   # "black"  #"light-blue"
                ),
                
                h4 ( "  VCF File validation "),
                fluidRow(
                  verbatimTextOutput("validation")
                  
                ),
                
                
                box(
                  radioButtons("disp", "Below is a Part of  the FORMAT column", choices = c(Head = "head"),  selected = "head" ), # downloadFORMAT
                  
                  style = "overflow-x: auto;",
                  tableOutput("contents") , width = "auto",  autowidth = TRUE  ,status = "primary", solidHeader = FALSE , collapsible = TRUE)
              ),      #   downloadButton( "downloadFORMAT", "Download Table"),
              
              # downloadsummary
              h4("A summary of the VCF file"),
              
              fluidRow(
                #           box(
                
                box(tableOutput("summary"), width = 12 ,status = "primary", solidHeader = FALSE , collapsible = TRUE )
                
              ), 
              downloadButton( "downloadsummary", "Download Table"),
              
             
      ),
      tabItem(tabName = "page2",
              h2("Read Depth Plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("Violin plot of read depth (DP) for the samples in the VCF file. A numeric matrix was produced from the variant call format (VCF) file.")
                  ) , width = 12
                ),
                box(plotOutput("violinPlot"), width = 12),            downloadButton("downloadViolinPlot", "Download Plot"),
                
              )
      ),
      tabItem(tabName = "page3",
              h2("Read Depth as a Box Plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("use the sequence depth (DP) from the VCF FORMAT column. We can isolate matrices of sequence depth.")
                  ) , width = 12
                ),
                box(plotOutput("DPPlot"), width = 12),
                downloadButton("downloadDPPlot", "Download Plot")
              )
      ),
      tabItem(tabName = "page4",
              h2("Genotype Plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("This Plot shows the genotype of all samples in the VCF with their unique variant, where the x-axis presents the samples and the Y-axis shows the GT.This plot was obtained from VCF genotype data using the extract.gt:vcfR function element = 'GT' to extract the genotype from the extracted matrix")
                  ) , width = 12
                ),
                box(plotOutput("genotypePlot"), width = 12),
                downloadButton("downloadGenotypePlot", "Download Plot")
              )
      ),
      tabItem(tabName = "page5",
              h2("Genotype Quality Plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("This plot visualizing genotype quality (GQ) scores to assess the reliability of genotype calls. This plot was obtained from VCF genotype data using the extract.gt function, from the arguments we select: element = 'GQ' to extract genotype quality from the matrix.")
                  ) , width = 12
                ),
                box(plotOutput("GQPlot"), width = 12),
                downloadButton("downloadGQPlot", "Download Plot")
              )
      ),
      
      tabItem(tabName = "page6",
              h2("Heatmap of Read Depth Plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("Heat map of read depth (DP) for the samples in your VCF file. Each column is a sample, and each row is a variant. The colour of each cell corresponds to each variant’s read depth (DP). Cells in white contain missing or zero data. Marginal bar plots summarize row and column sums and aid in identifying variants or samples of low depth.")
                  ) , width = 12
                ),
                box(plotOutput("HeatmapPlot"), width = 12),
                downloadButton("downloadHeatmapPlot", "Download Plot")
              )
      ),
      tabItem(tabName = "page7",
              h2("Bar Plot of the actual values of Read Depth"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("Bar plot of read depth (DP) for the samples in your VCF file. where X axes is a sample, and Y axes is a variant. Here we see that the mean number of high-quality bases per variant is shown for the samples. Some sample appears to have an abundance of missing data in the heatmap. In the barplot, we’ve validated that this sample lacks information. While some samples might have a high average coverage in the heatmap, and in the barplot, we can see the actual value.")
                  ) , width = 12
                ),
                box(plotOutput("barPlot"), width = 12),
                downloadButton("downloadBarPlot", "Download Plot")
              )
      ),
      
      
      #########################
      
      tabItem(tabName = "page8",
              h2("Multi-Layered Summary of Quality Metrics"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("This plot summarizes data from your VCF file. This track summarizes the number of variants per window. Above that, we see dot plots for quality, mapping quality, and read depth.")
                  ) , width = 12
                ),
                box(plotOutput("chromoplot"), width = 12),
                downloadButton("downloadChromoPlot", "Download Plot")
              )
      ),
      tabItem(tabName = "page9",
              h2("Histogram Summary of Quality Metrics"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("Here we see the distribution of some quality metrics. The raw read depth and mapping quality have been extracted from the INFO column of the VCF data. The quality is from the QUAL column of the VCF data. And the variant count per window was summarized during the windowing analysis performed by proc.chromR().")
                  ) , width = 12
                ),
                box(plotOutput("chromplot"), width = 12),
                downloadButton("downloadChromPlot", "Download Plot")
              )
      ),
      tabItem(tabName = "page10",
              h2("Allele Frequencies"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("This heatmap exhibits allele frequencies across samples and variants. It serves as a tool for representing allele frequency distribution and identifying potential rare or common variants. allele frequency positions are stored in rows each color represents a different value and individual samples are stored in columns. The plots utilize the AD_frequency function, from the vcfR package to generate allele frequencies based on matrices of allelic depths (AD).")
                  ) , width = 12
                ),
                box(plotOutput("AlleleFreqplot"), width = 12),
                downloadButton("downloadAlleleFreqPlot", "Download Plot")
              )
      ),
      tabItem(tabName = "page11",
              h2("Allele Frequencies For Samples"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("This box plot displays allele frequencies offering a depiction of how allele frequencies are distributed in the dataset, however, it will represent the mean value of AF in the samples. The plots utilize the AD_frequency function, from the vcfR package to generate allele frequencies based on matrices of allelic depths (AD).")
                  ) , width = 12
                ),
                box(plotOutput("AlleleFreqSamplesplot"), width = 12 ),
                downloadButton("downloadAlleleFreqSamplesplot", "Download Plot")
              )
      ),
     
      # # MAF  #
      tabItem(tabName = "page13",
              h2("Maftool Summary"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("This feature gives researchers a snapshot of gene, sample, and metadata details using the function maftool:: getSampleSummary for sample summary, getGeneSummary to Show gene summary, and getClinicalData to display clinical data associated with samples. This information is derived from the ANNOVAR output offering a dataset summary.")
                  ) , width =  "auto"
                ),
                h4("Below is the Sample Summary"),
                # DTOutput
                box(
                  
                  style = "overflow-x: auto;",
                  DTOutput ("getSampleSummary"), width = "auto" ,status = "primary", solidHeader = FALSE , collapsible = TRUE ),
                downloadButton( "downloadSampleSummary", "Download Sample Summary Table"),
                
                #### downloadgeneSummary
                h4("Gene Summary for the top 10 genes"),
                box(
                  style = "overflow-x: auto;",
                  DTOutput  ("geneSummary"), width = "auto" ,status = "primary", solidHeader = FALSE , collapsible = TRUE ),
                downloadButton( "downloadgeneSummary", "Download Gene Summary Table"),
                # 
                h4("Clinical Summary"),
                box(
                  style = "overflow-x: auto;",
                  DTOutput  ("clinicalSummary"), width = "auto" ,status = "primary", solidHeader = FALSE , collapsible = TRUE ),
                downloadButton( "downloadclinicalSummary", "Download Clinical Summary Table"),
                

              )
              
      ),
      
      tabItem(tabName = "page14",
              h2("Maftool Summary Plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li(" This detailed plot created with the function plotmafSummary:: maftools visually represents mutation annotation data providing insights into the mutation landscape within the dataset.")
                  ) , width = 12
                ),
                box(plotOutput("ANNOVARsummaryplot")
                    , width = 12) ,
                downloadButton( "downloadANNOVARsummary", "Download Plot")
                
                
                
              )
      ),
      # oncoplot  clindata2
      tabItem(tabName = "page15",
              h2("Oncoplot for top ten mutated genes"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("This heatmap was generated using oncoplot function to display the ten frequently mutated genes highlighting significant genetic changes in the dataset.")
                  ) , width = 12
                ),
                #######
                ### menuItem   clindata
                selectInput(inputId = "clindata2",
                            label = "Choose a column name of your clinical data:",
                            choices = NULL) ,
                selectInput(inputId = "clindata1",
                            label = "Choose a column name of your clinical data:",
                            choices = NULL) ,
                #####
                box(plotOutput("oncoplot")
                    , width = 12)
                , downloadButton("downloadOncoplot", "Download Plot")
                
              )
      ),
      
      tabItem(tabName = "page16",
              h2("plot titv summary"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li(" This plot illustrates the Ti/Tv ratio, a metric that indicates the balance between transition mutations (purine to purine or pyrimidine to pyrimidine) and transversion mutations (purine to pyrimidine or vice versa). This plot was done using the function plotTiTv")
                  ) , width = 12
                ),
                box(plotOutput("plotTiTv")
                    , width = 12),
                downloadButton("downloadTiTv", "Download Plot")
                
              )
      ),
      tabItem(tabName = "page17",
              h2("Lollipop plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li(" The lollipop plot offers an in-depth view of mutations, in specific genes of interest using the lollipopPlot function, by showing their position and frequency along the gene structure.")
                  ) , width = 12
                ),
                selectInput(inputId = "genename",
                            label = "Choose Gene Name:",
                            choices = NULL) ,
                
                box(plotOutput("lollipopPlot")
                    , width = 12),
                downloadButton("downloadLollipopPlot", "Download Plot")
                
              )
      ),
      tabItem(tabName = "page18",
              h2("Rainfall Plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("The rainfall chart shows how mutations are spread across the genome, helping to pinpoint areas with mutation rates or clusters of mutations. This plot was generated using rainfallPlot:: maftools function.")
                  ) , width = 12
                ),
                selectInput(inputId = "sample_names",
                            label = "Choose Sample Name:",
                            choices = NULL) ,
                
                box(plotOutput("rainfallPlot")
                    , width = 12),
                downloadButton("downloadRainfallPlot", "Download Plot")
                
              )
      ),
      tabItem(tabName = "page19",
              h2("Compare Mutational Load Against TCGA Cohorts"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li(" This chart is generated using the tcgaCompare::maftools function to compare the amount of mutations in the dataset with The Cancer Genome Atlas (TCGA) groups providing a reference point for the mutation levels seen in the study.")
                  ) , width = 12
                ), #
                textInput(inputId = "cohortName",
                          label = "Enter the name of your data:",
                          value = ""), #
                box(plotOutput("Comparemutation")
                    , width = 12),
                downloadButton("downloadComparemutation", "Download Plot")
                
              )
      ),
      tabItem(tabName = "page20",
              h2("Oncogenic Pathways plot"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li(" This visual representation indicates which cancer-causing pathways are influenced by the observed mutations giving insights into the effects of these genetic changes. We utilized the OncogenicPathways::maftools function to create this plot.")
                  ) , width = 12
                ),
                box(plotOutput("Pathways")
                    , width = 12),
                downloadButton("downloadPathways", "Download Plot")
                
                
              )
      ),
      tabItem(tabName = "contact",
              h2("Contact Us"),
              fluidRow(
                box(
                  tags$ul(
                    tags$li("Have questions or need to report an issue with the service? We've got you covered."),
                    tags$li("Contact Us at this email : XVCF@7OMICS.com")
                  ),
                  width = 12
                )
              )
      )
    )
  )
)


# Server logic for the app
server <- function(input, output, session) {
  observeEvent(input$vcfFile, {
    req(input$vcfFile)
    vcf_data(read.vcfR(input$vcfFile$datapath, verbose = FALSE))
    
    # Direct the user to the "VCF Summary" page (home)
    updateTabItems(session, "tabsVCF", "home")
  })
  
 
##########################################################################
  ## 1. Validation and Dynamic UI for ANNOVAR File Upload ##
  
  # Monitor `num_samples` and conditionally show/hide fileInput

  output$fileInputUI <- renderUI({
    # Validate input$num_samples to ensure it is not NULL and within range
    if (!is.null(input$num_samples) && is.numeric(input$num_samples) &&
        input$num_samples >= 1 && input$num_samples <= 200) {
      fileInput(
        inputId = "txtFile",
        label = "Choose ANNOVAR File (Tabular Format)",
        multiple = TRUE,
        accept = ".txt"
      )
    } else {
      NULL  # Do not show fileInput if num_samples is invalid
    }
  })
 
  ## Please enter the number of samples first to upload the file.##
  
  output$validationMessage <- renderText({
    if (is.null(input$num_samples) || !is.numeric(input$num_samples) ||
        input$num_samples < 1 || input$num_samples > 200)  { "....    Between 1-200 "
    } else {
      ""  # Clear the message if input is valid
    }
  })
  ## 2. Reactive Values and Observers ##
  
  # Observe ANNOVAR file upload and process data
  observeEvent(input$txtFile, {
    req(input$txtFile)
    ANNOVAR_data(read.table(
      input$txtFile$datapath,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  })
  ###
  observeEvent(input$txtFile, {
    req(input$num_samples)  # Ensure num_samples is provided
    req(input$num_samples >= 1 && input$num_samples <= 200)  # Validate range
    req(input$txtFile)  # Ensure the file is uploaded
    
    # Process ANNOVAR data
    ANNOVAR_data(read.table(
      input$txtFile$datapath,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
    # Direct the user to the "ANNOVAR Summary" page (page13)
    updateTabItems(session, "tabs", "page13")
  })
  ## 3. Reactive Data Processing ##
  
  # Reactive value for sample column names
  sample_column_names <- reactive({
    req(ANNOVAR_data())  # Ensure the ANNOVAR file is uploaded
    num_samples <- input$num_samples
    tail(names(ANNOVAR_data()), num_samples)  # Extract the last `num_samples` column names
  })
  # ANNOVAR summary processing
  ANNOVARsummary <- reactive({
    req(input$txtFile, sample_column_names())
    
    # Convert ANNOVAR data to MAF
    annovar_data <- annovarToMaf(
      annovar = input$txtFile$datapath,
      Center = "your data",
      tsbCol = sample_column_names(),
      table = "refGene",
      ens2hugo = TRUE,
      sep = "\t"
    )
    
    # Transform data
    ball_maf1 <- annovar_data %>%
      pivot_longer(
        cols = all_of(sample_column_names()),
        names_to = "Tumor_Sample_Barcode",
        values_to = "Genotype"
      )
    
    # Replace empty genotypes with NA
    ball_maf1$Genotype <- gsub("^\\.\\/.+?\\.$", "NA", ball_maf1$Genotype, perl = TRUE)
    
    # Filter out rows with NA genotypes
    ball_maf2 <- ball_maf1[ball_maf1$Genotype != "NA", ]
    
    # Validate and apply filters
    if (!"1000G_ALL" %in% names(ball_maf2)) {
      stop("The column `1000G_ALL` does not exist in the data.")
    }
    
    ball_maf3 <- ball_maf2[ball_maf2$`1000G_ALL` <= AF_Threshold(), ]
    return(ball_maf3)
  })
  
  ## 4. Output Rendering ##
  
  # Render the extracted sample column names
  output$sample_column_names <- renderTable({
    req(sample_column_names())  # Ensure sample names are available
    data.frame(Sample_Names = sample_column_names())  # Display as a table
  })
  
  
  ########################################################################
 
  
  # Read clinical data xlsx 
  observeEvent(input$xlsxFile, {
    req(input$xlsxFile)
    clindata(read_xlsx(input$xlsxFile$datapath,sheet = 1))
  })
  
  
  #################################################################
  
  
  # VCF validation
  
  output$validation <- renderText({
    req(vcf_data())
    vcf_header <- vcf_data()@meta
    if (!grepl("^##fileformat", vcf_header[1])) {
      stop("Error: This does not appear to be a valid VCF file.")
    }
    
    required_headers <- c("##fileformat", "##INFO", "##FORMAT", "##FILTER", "##contig")
    missing_headers <- required_headers[lengths(sapply(required_headers, grep, x = vcf_data()@meta)) == 0]
    
    if (length(missing_headers) > 0) {
      error_message <- paste("Error: The following required VCF headers are missing:", paste(missing_headers, collapse = ", "))
      stop(error_message)
    } else {
      success_message <-  "  Validation successful. This is a valid VCF file."
      cat(success_message, "\n")
      return(success_message)
    }
  })
  
  ######################### VCF summary  ##################################################
  
  output$summary <- renderTable({
    
    req(vcf_data())
    vcf_header <- vcf_data()@meta
    
    reference_line <- grep("^##reference=", vcf_header, value = TRUE)
    reference <- ifelse(length(reference_line) > 0, gsub("^##reference=", "", reference_line), NA)
    
    fileformat_line <- grep("^##fileformat=", vcf_header, value = TRUE)
    fileformat <- ifelse(length(fileformat_line) > 0, gsub("^##fileformat=", "", fileformat_line), NA)
    
    # Extract all contig lines
    chromosome_name_lines <- grep("^##contig=<ID=", vcf_header, value = TRUE)
    contigs <- gsub("^##contig=<ID=", "", chromosome_name_lines)
    # Number of sample
    num_samples = ncol(vcf_data()@gt) - 1
    # Number of Variants
    num_variants = nrow(vcf_data()@gt)
      # Is the VCF Phased or Unphased 
    num_phased <- sum(grepl(pattern="\\|",x= vcf_data()@gt[1,])) # sum(TRUE, FALSE, ..)
  
    phased <- if (num_phased == 0) {
      "The variant is Unphased."
    } else {
      "The variant is Phased."
    }
    
    # SNP INS DEL
    ref_column <- vcf_data()@fix[, "REF", drop = FALSE]
    alt_column <- vcf_data()@fix[, "ALT", drop = FALSE]
    
    num_snps <- sum(nchar(ref_column) == 1 & nchar(alt_column) == 1)
    num_insertions <- sum(nchar(ref_column) < nchar(alt_column))
    num_deletions <- sum(nchar(ref_column) > nchar(alt_column))
    
    # Check for annotation
    annotation <- vcf_data()@meta[grep("^##annotationSources=", vcf_data()@meta)]
    annotation_status <- if (length(annotation) > 0) {
      "The VCF file is annotated."
    } else {
      "The VCF file is not annotated."
    }
    # Extract sample names
    samples <- paste(colnames(vcf_data()@gt)[-1], collapse = ", ")
    
    
    # Create a data frame with column names and values
    metadata <- data.frame(
      Type.Of.Data = c("Reference", "File Format", "Contig", "Number of SNPs", "Number of Insertions", "Number of Deletions", "Number of Samples", "Number of Variants", "Phased", "Is The VCF File Annotated?", "Samples Name"),
      Description = c(reference, fileformat, paste(contigs, collapse = ", ") , num_snps, num_insertions, num_deletions, num_samples, num_variants, phased, annotation_status, samples),
      stringsAsFactors = FALSE
    )
    
    return(metadata)
  })
  
  output$downloadsummary <- downloadHandler(
    filename = function() {
      "ANNOVAR_VCF_summary_table.csv"
    },
    content = function(file) {
      vcf_header <- vcf_data()@meta
      
      reference_line <- grep("^##reference=", vcf_header, value = TRUE)
      reference <- ifelse(length(reference_line) > 0, gsub("^##reference=", "", reference_line), NA)
      
      fileformat_line <- grep("^##fileformat=", vcf_header, value = TRUE)
      fileformat <- ifelse(length(fileformat_line) > 0, gsub("^##fileformat=", "", fileformat_line), NA)
      
      # Extract all contig lines
      chromosome_name_lines <- grep("^##contig=<ID=", vcf_header, value = TRUE)
      contigs <- gsub("^##contig=<ID=", "", chromosome_name_lines)
      
      num_samples = ncol(vcf_data()@gt) - 1
      num_variants = nrow(vcf_data()@gt)
      
      num_phased <- sum(grepl(pattern="\\|",x= vcf_data()@gt[1,])) # sum(TRUE, FALSE, ..)
      phased <- if (num_phased == 0) { "The variant is Unphased."  } else {  "The variant is Phased."  }
      
      ref_column <- vcf_data()@fix[, "REF", drop = FALSE]
      alt_column <- vcf_data()@fix[, "ALT", drop = FALSE]
      
      num_snps <- sum(nchar(ref_column) == 1 & nchar(alt_column) == 1)
      num_insertions <- sum(nchar(ref_column) < nchar(alt_column))
      num_deletions <- sum(nchar(ref_column) > nchar(alt_column))
      
      #annotation 
    
      # Check for annotation
      annotation <- vcf_data()@meta[grep("^##annotationSources=", vcf_data()@meta)]
      annotation_status <- if (length(annotation) > 0) {  "The VCF file is annotated."  } else {  "The VCF file is not annotated." }
      
      # Extract sample names
      samples <- paste(colnames(vcf_data()@gt)[-1], collapse = ", ")
      
      # Create a data frame with column names and values
      metadata <- data.frame(
        Variable = c("Reference", "File Format", "Contig", "Number of SNPs", "Number of Insertions", "Number of Deletions", "Number of Samples", "Number of Variants", "Phased", "Is The VCF file annotated?", "Samples Name"),
        Value = c(reference, fileformat, paste(contigs, collapse = ", "), num_snps, num_insertions, num_deletions, num_samples, num_variants, phased, annotation_status, samples),
        stringsAsFactors = FALSE
      )
      
      # return(metadata)
      # Save the data frame to a CSV file
      write.csv(metadata, file)
    }
  )
  
  ##################################################################################
 
  # FORMAT
  output$contents <- renderTable({
    req(vcf_data())
    df <- as.data.frame(vcf_data()@gt)
    
    if (input$disp == "head") {
      return(head(df))
    } else {
      return(df)
    }
  })
  
  #  download FORMAT
  output$downloadFORMAT <- downloadHandler(
    filename = function() {
      "ANNOVAR_FORMAT_Head_table.csv"
    },
    content = function(file) {
      df <- as.data.frame(vcf_data()@gt)
      
      if (input$disp == "head") {
        return(head(df))
      } else {
        return(df)
      }
      # Save the data frame to a CSV file
      write.csv(df, file)
    }
  )
  
  
  
  ####### violinPlot
  dpf_data <- reactive({
    req(vcf_data())
    dp <- extract.gt(vcf_data(), element = 'DP', as.numeric = TRUE)
    dpf <- melt(dp, varnames = c('Index', 'Sample'), value.name = 'Depth', na.rm = TRUE)
    dpf[dpf$Depth > 0, ]
  })
  
  output$violinPlot <- renderPlot({
    p <- ggplot(dpf_data(), aes(x = Sample, y = Depth)) +
      geom_violin(fill = '#C0C0C0', adjust = 1.0, scale = 'count', trim = TRUE) +
      theme_bw() +
      ylab('Read Depth (DP)') +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) +
      stat_summary(fun.data = mean_sdl, geom = 'pointrange', color = 'black') +
      scale_y_continuous(trans = scales::log2_trans(), breaks = c(1, 10, 100, 1000))
    
    print(p)
  })
  
  output$downloadViolinPlot <- downloadHandler(
    filename = function() {
      "violin_plot_of_Read_depth.png"  
    },
    content = function(file) {
      png(file)
      p <- ggplot(dpf_data(), aes(x = Sample, y = Depth)) +
        geom_violin(fill = '#C0C0C0', adjust = 1.0, scale = 'count', trim = TRUE) +
        theme_bw() +
        ylab('Read Depth (DP)') +
        theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) +
        stat_summary(fun.data = mean_sdl, geom = 'pointrange', color = 'black') +
        scale_y_continuous(trans = scales::log2_trans(), breaks = c(1, 10, 100, 1000))
      print(p)
      dev.off()
    }
  )
  
  
  # Sequence Depth (DP) box plot
  dp_data <- reactive({
    req(vcf_data())
    dp <- extract.gt(vcf_data(), element = "DP", as.numeric = TRUE)
    dp2 <- dp
    dp2[dp2 == 0] <- NA
    dp2
  })
  
  output$DPPlot <- renderPlot({
    boxplot(dp_data(), las = 2, col = 2:5, main = "Sequence Depth (DP)", log = "y")
    abline(h = 10^c(0:4), lty = 3, col = "#808080")
  })
  
  output$downloadDPPlot <- downloadHandler(
    filename = function() {
      "boxplot_depth.png"  
    },
    content = function(file) {
      png(file)
      boxplot(dp_data(), las = 2, col = 2:5, main = "Sequence Depth (DP)", log = "y")
      abline(h = 10^c(0:4), lty = 3, col = "#808080")
      dev.off()
    }
  )
  
  
  
  # genotype Plot 
  
  geno_data <- reactive({
    req(vcf_data())
    geno_matrix <- extract.gt(vcf_data())
    # Convert Genotypes to a data frame for ggplot
    df_Geno <- as.data.frame(geno_matrix)
    df_Geno$Variant <- rownames(df_Geno)
    # Melt the data frame for plotting
    df_genotypes <- melt(df_Geno, id.vars = "Variant", variable.name = "Samples", value.name = "Genotypes")
    df_genotypes
  })
  
  output$genotypePlot <- renderPlot({
    ggplot(geno_data(), aes(x = Samples, fill = Genotypes)) +
      geom_bar(position = "stack") +
      labs(title = "Genotype Plot of All Samples", x = "Samples", y = "Variant") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_blank())
  })
  
  output$downloadGenotypePlot <- downloadHandler(
    filename = function() {
      "genotype_plot.png"  
    },
    content = function(file) {
      png(file)
      ggplot(geno_data(), aes(x = Samples, fill = Genotypes)) +
        geom_bar(position = "stack") +
        labs(title = "Genotype Plot of All Samples", x = "Samples", y = "Variant") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_blank())
      dev.off()
    }
  )
  
  
  # Genotype Quality plot
  gq_data <- reactive({
    req(vcf_data())
    gq <- extract.gt(vcf_data(), element = "GQ", as.numeric = TRUE)
    gq
  })
  
  output$GQPlot <- renderPlot({
    par(mar = c(8, 4, 4, 2))
    boxplot(gq_data(), las = 2, col = 2:5, main = "Genotype Quality (GQ)")
  })
  
  output$downloadGQPlot <- downloadHandler(
    filename = function() {
      "Genotype_Quality_plot.png"  
    },
    content = function(file) {
      png(file)
      par(mar = c(8, 4, 4, 2))
      boxplot(gq_data(), las = 2, col = 2:5, main = "Genotype Quality (GQ)")
      dev.off()
    }
  )
  
  
  # DP heatmap plot
  heatmap_data <- reactive({
    req(vcf_data())
    dp <- extract.gt(vcf_data(), element = "DP", as.numeric = TRUE)
    dp
  })
  
  output$HeatmapPlot <- renderPlot({
    heatmap.bp(heatmap_data()[1001:1500, ])
  })
  
  output$downloadHeatmapPlot <- downloadHandler(
    filename = function() {
      "RD_heatmap_plot.png"  
    },
    content = function(file) {
      png(file)
      heatmap.bp(heatmap_data()[1001:1500, ])
      dev.off()
    }
  )
  
  
  #
  barplot_data <- reactive({
    req(vcf_data())
    dp <- extract.gt(vcf_data(), element = "DP", as.numeric = TRUE)
    #rownames(dp) <- 1:nrow(dp)
    dp
  })
  
  output$barPlot <- renderPlot({
    par(mar = c(8, 4, 4, 2))
    barplot(apply(barplot_data(), MARGIN = 2, mean, na.rm = TRUE), las = 3)
  })
  
  output$downloadBarPlot <- downloadHandler(
    filename = function() {
      "Read_depth_bar_plot.png"
    },
    content = function(file) {
      png(file)
      par(mar = c(8, 4, 4, 2))
      barplot(apply(barplot_data(), MARGIN = 2, mean, na.rm = TRUE), las = 3)
      dev.off()
    }
  )
  #
  
  
  chromo_data <- reactive({
    req(vcf_data())
    chrom <- create.chromR(name = "Chromosome 1", vcf = vcf_data(), seq = NULL, ann = NULL, verbose = TRUE)
    chrom <- proc.chromR(chrom, verbose = TRUE)
    chromoqc(chrom, dp.alpha = 22)
    chrom
  })
  
 
  output$chromoplot <- renderPlot({
    chromo_data()
  })
  
  output$downloadChromoPlot <- downloadHandler(
    filename = function() {
      "Summarized_plot_of_Ch1.png"  
    },
    content = function(file) {
      png(file)
      chromo_data()
      dev.off()
    }
  )
  
  #
  chrom_data <- reactive({
    req(vcf_data())
    chrom <- create.chromR(name = "Chromosome 1", vcf = vcf_data(), seq = NULL, ann = NULL, verbose = TRUE)
    # chrom <- masker(chrom, min_QUAL = 0, min_DP = 350, max_DP = 650, min_MQ = 59.5, max_MQ = 60.5)
    chrom <- proc.chromR(chrom, verbose = TRUE)
    chromoqc(chrom, dp.alpha = 22)
    #chrom <- proc.chromR(chrom, verbose = FALSE, win.size = 1e3)
    chrom
  })
  
  output$chromplot <- renderPlot({
    plot(chrom_data())
  })
  
  output$downloadChromPlot <- downloadHandler(
    filename = function() {
      "qualit_ metrics_of_Ch1_plot.png"  
    },
    content = function(file) {
      png(file)
      plot(chrom_data())
      dev.off()
    }
  )
  
  
  # Allele Freq plot
  AlleleFreqplot <- reactive({
    req(vcf_data())
    
    
    # Extract the AD (Allelic Depth) matrix
    ad_matrix_AO <- extract.gt(vcf_data(), element = "AO")
    ad_matrix_AD <- extract.gt(vcf_data(), element = "AD")
    
    # Check if "AO" is present
    if (!is.null(ad_matrix_AO) && !all(is.na(ad_matrix_AO))) {
      # Calculate allele frequencies using AO
      allele_frequencies <- AD_frequency(ad_matrix_AO, delim = ",", allele = 1L, sum_type = 0L, decreasing = 1L)
    } else if (!is.null(ad_matrix_AD) && !all(is.na(ad_matrix_AD))) {
      # Calculate allele frequencies using AD if AO is not available
      allele_frequencies <- AD_frequency(ad_matrix_AD, delim = ",", allele = 1L, sum_type = 0L, decreasing = 1L)
    } else {
      # Handle the case where neither "AO" nor "AD" is present
      stop("Error: Both 'AO' and 'AD' elements are missing or contain only 'NA' values.")
    }
    
    
    # Convert allele frequencies to a data frame for ggplot
    df <- as.data.frame(allele_frequencies)
    df$Variant <- rownames(df)
    # Melt the data frame for plotting
    df_melted <- melt(df, id.vars = "Variant", variable.name = "Samples", value.name = "Frequency")
    # Plotting AD_frequency as a heatmap
    ggplot(df_melted, aes(x = Variant, y = Samples, fill = Frequency)) +
      geom_tile() +  # Using geom_tile for a heatmap
      labs(title = "Allele Frequencies Heatmap", x = "Variant", y = "Samples") +
      theme(axis.text.x = element_blank())
    
  })
  output$AlleleFreqplot <- renderPlot({
    plot(AlleleFreqplot())
  })
  output$downloadAlleleFreqPlot <- downloadHandler(
    filename = function() {
      "allele_frequency_plot.png"  
    },
    content = function(file) {
      png(file)
      plot(AlleleFreqplot())
      dev.off()
    }
  )
  
  # Allele Freq Samples plot
  
  AlleleFreqSamplesplot <- reactive({
    req(vcf_data())
    
    
    # Extract the AD (Allelic Depth) matrix
    ad_matrix_AO <- extract.gt(vcf_data(), element = "AO")
    ad_matrix_AD <- extract.gt(vcf_data(), element = "AD")
    
    # Check if "AO" is present
    if (!is.null(ad_matrix_AO) && !all(is.na(ad_matrix_AO))) {
      # Calculate allele frequencies using AO
      allele_frequencies <- AD_frequency(ad_matrix_AO, delim = ",", allele = 1L, sum_type = 0L, decreasing = 1L)
    } else if (!is.null(ad_matrix_AD) && !all(is.na(ad_matrix_AD))) {
      # Calculate allele frequencies using AD if AO is not available
      allele_frequencies <- AD_frequency(ad_matrix_AD, delim = ",", allele = 1L, sum_type = 0L, decreasing = 1L)
    } else {
      # Handle the case where neither "AO" nor "AD" is present
      stop("Error: Both 'AO' and 'AD' elements are missing.")
    }
    
    # Convert allele frequencies to a data frame for ggplot
    df <- as.data.frame(allele_frequencies)
    df$Variant <- rownames(df)
    # Melt the data frame for plotting
    df_melted <- melt(df, id.vars = "Variant", variable.name = "Samples", value.name = "Frequency")
    # Plotting AD_frequency as a boxplot, X = Samples
    ggplot(df_melted, aes(x = Samples, y = Frequency)) +
      geom_boxplot() +  # Using geom_boxplot for a box plot
      labs(title = "Allele Frequencies Boxplot", x = "Samples", y = "Frequency") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  })
  output$AlleleFreqSamplesplot <- renderPlot({
    plot(AlleleFreqSamplesplot())
  })
  
  # Add a download button for the plot
  output$downloadAlleleFreqSamplesplot <- downloadHandler(
    filename = function() {
      "allele_freq_samples_plot.png"
    },
    content = function(file) {
      png(file)
      plot(AlleleFreqSamplesplot())
      dev.off()
    }
  )

  
  #########################################################################
  #                             ANNOVAR                                   ##
  ##########################################################################
  
  ############################### Annovar Summary  #########################################
  ## renderDT  renderDataTable
  output$getSampleSummary <- renderDT({
   # req(ANNOVAR_data())
    req(ANNOVARsummary())
    
    
    ## reading MAF file with clinical data
    ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    
    # Get the sample summary
    sampleSummary <- getSampleSummary(ball)
     datatable(
      sampleSummary,
      
      options = list(
        scrollX = TRUE,
        scrollY = "250px",
        
        pageLength = 5,
        lengthMenu = c(5, 10, 15),
        dom = 'tp',
        columnDefs = list(
          list(width = '100px', targets = '_all')
        )
      )
    )
    return(sampleSummary)
  })
  
  output$downloadSampleSummary <- downloadHandler(
    filename = function() {
      "ANNOVAR_sample_Summary_table.csv"
    },
    content = function(file) {
      ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      # Get the sample summary
      sampleSummary <- getSampleSummary(ball)
      
     
      # Save the data frame to a CSV file
      write.csv(sampleSummary, file)
    }
  )
  
  
  #####  gene summary.########
  output$geneSummary <- renderDT({
   # req(ANNOVAR_data())
    req(ANNOVARsummary())
    
    
    ## reading MAF file with clinical data
    ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    
    
    #Shows gene summary.
    GeneSummary <- getGeneSummary(ball)
    
    
    # Select and display the top 10 genes
    topGenes <- head(GeneSummary[order(GeneSummary$total, decreasing = TRUE), ], 10)
    
    datatable(
      topGenes,
      
      options = list(
        scrollX = TRUE,
        scrollY = "250px",
        
        pageLength = 5,
        lengthMenu = c(5, 10, 15),
        dom = 'tp',
        columnDefs = list(
          list(width = '100px', targets = '_all')
        )
      )
    )
    return(topGenes)
  })
  
  
  
  ##
  
  output$downloadgeneSummary <- downloadHandler(
    filename = function() {
      "ANNOVAR_gene_Summary_table.csv"
    },
    content = function(file) {
      ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      # Shows gene summary.
      GeneSummary <- getGeneSummary(ball)
      # Select and display the top 10 genes
      topGenes <- head(GeneSummary[order(GeneSummary$total, decreasing = TRUE), ], 10)
      
      # Save the data frame to a CSV file
      write.csv(topGenes, file)
    }
  )
  

  
  #####  clinical summary.########
  output$clinicalSummary <- renderDT({
   # req(ANNOVAR_data())
    req(ANNOVARsummary())
    
    
    ## reading MAF file with clinical data
    ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    
    #shows clinical data associated with samples
    getClinicalData(ball)
    
    ClinicalData <- getClinicalData(ball)
    
    datatable(
      ClinicalData,
      
      options = list(
        scrollX = TRUE,
        scrollY = "250px",
        
        pageLength = 5,
        lengthMenu = c(5, 10, 15),
        dom = 'tp',
        columnDefs = list(
          list(width = '100px', targets = '_all')
        )
      )
    )
    
  })
  
  output$downloadclinicalSummary <- downloadHandler(
    filename = function() {
      "ANNOVAR_clinical_Summary_table.csv"
    },
    content = function(file) {
      ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      #shows clinical data associated with samples
      getClinicalData(ball)
      
      ClinicalData <- getClinicalData(ball)
      
      # Save the data frame to a CSV file
      write.csv(ClinicalData, file)
    }
  )

   ###############################
  
  # Reactive AF_Threshold
  AF_Threshold <- reactive({
    input$AF_Threshold
  })
  
  output$AF_Threshold <- renderText({
    req(AF_Threshold())
    if (is.na(AF_Threshold())) {
      "No threshold selected"
    } else {
      paste("AF Threshold:", AF_Threshold())
    }
  })
  
  
  ########################## plot 1 ########################################
  
  
  output$ANNOVARsummaryplot <- renderPlot({
    
    ## reading MAF file with clinical data
    ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    
    plotmafSummary(maf = ball, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, showBarcodes = TRUE, titvRaw = FALSE)
    
  })
  
  output$downloadANNOVARsummary <- downloadHandler(
    filename = function() {
      "ANNOVAR_summary_plot.png"  
    },
    content = function(file) {
      png(file)
      ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      plotmafSummary(maf = ball, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, showBarcodes = TRUE, titvRaw = FALSE)
      dev.off()
    }
  )
  
  ############################################################################
  
  # Dynamically update clindata names choices
  observe({
    req(clindata())
    
    clindata1 <- colnames(clindata())
    
    updateSelectInput(session, "clindata1", choices = clindata1)
  })
  
  observe({
    req(clindata())
    
    clindata2 <- colnames(clindata())
    
    updateSelectInput(session, "clindata2", choices = clindata2)
    
  })
  
  ##############################################################################

  output$oncoplot <- renderPlot({
    req(ANNOVARsummary()) 
    
    # Read the MAF file with or without clinical data
    clinical_data <- clindata() # May be NULL
    ball <- read.maf(maf = ANNOVARsummary(), clinicalData = clinical_data)
    
    # Get selected clinical feature columns
    selected_clndata1 <- input$clindata1
    selected_clndata2 <- input$clindata2
    
    # To Ensure selected columns exist in clinicalData (if provided)
    clinical_features <- if (!is.null(clindata()) &&
                             !is.null(selected_clndata1) &&
                             selected_clndata1 %in% colnames(clindata())) {
      if (!is.null(selected_clndata2) && selected_clndata2 %in% colnames(clindata())) {
        c(selected_clndata1, selected_clndata2)
      } else {
        selected_clndata1
      }
    } else {
      NULL
    }
    
    # Generate the plot
    if (is.null(clinical_features)) {
      # Plot without clinical annotations
      oncoplot(
        maf = ball,
        top = 20,
        showTumorSampleBarcodes = TRUE
      )
    } else {
      # Plot with clinical annotations
      oncoplot(
        maf = ball,
        top = 20,
        showTumorSampleBarcodes = TRUE,
        clinicalFeatures = clinical_features
      )
    }
  })
  
  ############
  
  output$downloadOncoplot <- downloadHandler(
    filename = function() {
      "oncoplot.png"
    },
    content = function(file) {
      png(file)
      #ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      #oncoplot(maf = ball, top = 20, showTumorSampleBarcodes = TRUE) #, clinicalFeatures = c("Progression", "Gender"))
      
      # Read the MAF file with or without clinical data
      clinical_data <- clindata() # May be NULL
      ball <- read.maf(maf = ANNOVARsummary(), clinicalData = clinical_data)
      
      # Get selected clinical feature columns
      selected_clndata1 <- input$clindata1
      selected_clndata2 <- input$clindata2
      
      # To Ensure selected columns exist in clinicalData (if provided)
      clinical_features <- if (!is.null(clindata()) &&
                               !is.null(selected_clndata1) &&
                               selected_clndata1 %in% colnames(clindata())) {
        if (!is.null(selected_clndata2) && selected_clndata2 %in% colnames(clindata())) {
          c(selected_clndata1, selected_clndata2) } else {  selected_clndata1  }  } else {   NULL  }
      
      # Generate the plot
      if (is.null(clinical_features)) {
        # Plot without clinical annotations
        oncoplot(
          maf = ball,
          top = 20,
          showTumorSampleBarcodes = TRUE
        )
      } else {
        # Plot with clinical annotations
        oncoplot(
          maf = ball,
          top = 20,
          showTumorSampleBarcodes = TRUE,
          clinicalFeatures = clinical_features
        )
      }
       dev.off()
    }
  )
  
  
  ########################## plot 3 ########################################
  output$plotTiTv <- renderPlot({
    
    ## reading MAF file with clinical data
    ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    
    #plot 3
    ball.titv = titv(maf = ball, plot = FALSE, useSyn = TRUE)
    #plot titv summary
    plotTiTv(res = ball.titv, showBarcodes = TRUE)
  })
  
  output$downloadTiTv <- downloadHandler(
    filename = function() {
      "plot_TiTv.png"
    },
    content = function(file) {
      png(file)
      ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      
      ball.titv = titv(maf = ball, plot = FALSE, useSyn = TRUE)
      #plot titv summary
      plotTiTv(res = ball.titv, showBarcodes = TRUE)
      dev.off()
    }
  )
  
  
  ########################## plot 4 ########################################
  # Dynamically update gene names choices
  observe({
   # req(ANNOVAR_data())
    req(ANNOVARsummary())
    
    ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    GeneSummary <- getGeneSummary(ball)
    genename <- GeneSummary$Hugo_Symbol
    
    
    updateSelectInput(session, "genename", choices = genename)
  })
  
  #########################################
  
  output$lollipopPlot <- renderPlot({
   # req(ANNOVAR_data())
    req(ANNOVARsummary())
    
    # reading MAF file with clinical data
    ball <- read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    
    # Get selected gene name
    selected_gene <- input$genename
    
    # Check if selected gene name is not empty
    if (!is.null(selected_gene)) {
      
      # lollipop plot
      lollipopPlot(maf = ball, gene = selected_gene , AACol = 'aaChange', showMutationRate = TRUE)
      
    } else {
      # If no gene name selected, show a message or default plot
      plot(NULL, type = "n", xlab = "", ylab = "", main = "No gene selected")
    }
  })
  
  # 
  output$downloadLollipopPlot <- downloadHandler(
    filename = function() {
      "lollipop_plot.png"
    },
    content = function(file) {
      png(file)
     
      ball <- read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      
      # Get selected gene name
      selected_gene <- input$genename
      
      # Check if selected gene name is not empty
      if (!is.null(selected_gene)) {
        
        # lollipop plot
        lollipopPlot(maf = ball, gene = selected_gene , AACol = 'aaChange', showMutationRate = TRUE)
        
      } else {
        # If no gene name selected, show a message or default plot
        plot(NULL, type = "n", xlab = "", ylab = "", main = "No gene selected")
      }
      dev.off()
    }
  )
  
  
  ########################## plot 5 ########################################
  
 # Dynamically update sample names choices
  observe({
   # req(sample_column_names())
    req(ANNOVAR_data())

    num_samples <- input$num_samples
    sample_names <- tail(names(ANNOVAR_data()), num_samples)
    updateSelectInput(session, "sample_names", choices = sample_names)
  })

  ################################
  output$rainfallPlot <- renderPlot({
   # req(ANNOVAR_data())
    req(ANNOVARsummary())


    # reading MAF file with clinical data
    ball <- read.maf(maf = ANNOVARsummary(), clinicalData = clindata())

    # Get selected sample name
    selected_sample <- input$sample_names

    # Check if selected sample name is not empty
    if (!is.null(selected_sample)) {
      # rainfall plot
      rainfallPlot(maf = ball, tsb = selected_sample, detectChangePoints = TRUE, pointSize = 0.6)
    } else {
      # If no sample name selected, show a message or default plot
      plot(NULL, type = "n", xlab = "", ylab = "", main = "No sample selected")
    }
  })

  
  ################ 
  output$downloadRainfallPlot <- downloadHandler(
    filename = function() {
      "rainfall_plot.png"
    },
    content = function(file) {
      png(file)
    
      # reading MAF file with clinical data
      ball <- read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      
      # Get selected sample name
      selected_sample <- input$sample_names
      
      # Check if selected sample name is not empty
      if (!is.null(selected_sample)) {
        # rainfall plot
        rainfallPlot(maf = ball, tsb = selected_sample, detectChangePoints = TRUE, pointSize = 0.6)
      } else {
        # If no sample name selected, show a message or default plot
        plot(NULL, type = "n", xlab = "", ylab = "", main = "No sample selected")
      }
      
      dev.off()
    }
  )
  ########################## plot 6 ########################################
  # reactive cohortName
  cohortName <- reactive({
    input$cohortName
    
  })
  
  output$cohortName <- renderText({
    req(cohortName())
    input$cohortName
    
    return(cohortName())
  })
  output$Comparemutation <- renderPlot({
    
    # ## reading MAF file with clinical data
    ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    
    ## Compare mutational load gainst TCGA cohorts
    ball.mutload = tcgaCompare(maf = ball, cohortName = cohortName())
    
  })
  
  output$downloadComparemutation <- downloadHandler(
    filename = function() {
      "Compare_mutation_plot.png"  
    },
    content = function(file) {
      png(file)
      ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      ball.mutload = tcgaCompare(maf = ball, cohortName = cohortName())
      
      # plotVaf(maf = ball.mutload, vafCol = 'i_TumorVAF_WU')
      dev.off()
    }
  )
  ########################## plot 7 ########################################
  output$Pathways <- renderPlot({
    ## reading MAF file with clinical data
    ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
    
    OncogenicPathways(maf = ball)
    
  })
  
  
  output$downloadPathways <- downloadHandler(
    filename = function() {
      "Pathways_plot.png"
    },
    content = function(file) {
      png(file)
      
      ball = read.maf(maf = ANNOVARsummary(), clinicalData = clindata())
      
      OncogenicPathways(maf = ball)
      dev.off()
    }
  )
  
}
# Run the application
shinyApp(ui = ui, server = server)



