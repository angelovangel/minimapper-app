library(shiny)
library(shinyWidgets)
library(shinyjs)
library(bslib)
library(bsicons)
library(stringr)
library(dplyr)
library(processx)
library(shinybusy)
library(digest)
library(shinyvalidate)
library(hover)

bin_on_path = function(bin) {
  exit_code = suppressWarnings(system2("command", args = c("-v", bin), stdout = FALSE))
  return(exit_code == 0)
}

sidebar <- sidebar(
  title = 'Inputs',
  
  selectInput(
    'format', 
    tags$a(
      'Select reference format',
      tooltip(bsicons::bs_icon("question-circle"),
        "Select correct format, it is not guessed from the input",
        placement = "right")
    ),
    choices = c('fasta', 'genbank', 'embl', 'snapgene'), multiple = FALSE
  ),
  fileInput(
    'upload_ref', 
    tags$a(
      'Upload reference',  
      tooltip(bsicons::bs_icon("question-circle"),
        "Upload reference file containing one sequence (fasta, genbank, embl or snapgene)",
        placement = "right")
    ), 
    multiple = F, accept = c('.fa', '.fasta', '.dna'), placeholder = 'reference'),
  
  # if more than 1, the uploaded files will be placed in a tmp folder on the server and
  # then passed to the nxf script as dir
  fileInput(
      'upload_fastq', 
    tags$a(
      'Upload fastq',  
      tooltip(bsicons::bs_icon("question-circle"),
        "Upload either one fastq file or more fastq files to be mapped to the reference",
        placement = "right")
    ), 
    multiple = T, accept = c('.fastq', '.gz', '.fq'), placeholder = 'fastq file(s)'),
  
  hover_action_button('start', 'Run mapping', icon = icon('play'), button_animation = 'icon-fade'),
  hover_action_button('reset', 'Reset inputs', icon = icon('rotate'), button_animation = 'icon-fade'),
  
)

ui <- page_navbar(
  header = 
    tags$a('Live terminal view',
    tooltip(
    bsicons::bs_icon("question-circle"),
    "Preview of the terminal, for viewing the selected parameters and monitor output",
    placement = "right")),
  
  useShinyjs(),
  use_hover(),
  
  fillable = T,
  title = 'nxf-minimapper app',
  theme = bs_theme(bootswatch = 'yeti', primary = '#196F3D'),
  sidebar = sidebar,
  #nav_panel(
    card(
      verbatimTextOutput('stdout')
    )
  #)
)

server <- function(input, output, session) {
  
  ##
  options(shiny.maxRequestSize=1000*1024^2)
  ##
  
  # check nextflow and docker is on path
  if (!bin_on_path('nextflow') | !bin_on_path('docker')) {
    notify_failure('nextflow and/or docker not found!', position = 'center-bottom')
  } else {
    notify_success('The server is ready!', position = 'center-bottom')
  }
  
  
  ref <- reactive({
    rname <- input$upload_ref$name
    rpath <- input$upload_ref$datapath
    
    fs::file_move(rpath, fs::path(fs::path_dir(rpath), rname))
    fs::path(fs::path_dir(rpath), rname)
  })
  
  fastq <- reactive({
    fnames <- input$upload_fastq$name
    fpaths <- input$upload_fastq$datapath
    
    fs::file_move(fpaths, fs::path(fs::path_dir(fpaths), fnames))
    fs::path_dir(fpaths[1])
  })
  
  output$stdout <- renderPrint({
    req(input$upload_fastq)
    req(input$upload_ref)
    
    # validations
    arguments <<- c('--ref', ref(), '--fastq', fastq(), '--format', input$format)
  })
  
  #main
  observeEvent(input$start, {
    #req(input$upload_fastq)
    
    withCallingHandlers({
      p <- processx::run(
        'nextflow',
        args = c('run', 'angelovangel/nxf-minimapper', arguments),
        echo_cmd = T,
        stderr_to_stdout = TRUE,
        error_on_status = FALSE,
        stdout_line_callback = function(line, proc) {
          message(line)
        }
      )
    }, message = function(m) {
      shinyjs::html(
        id = "stdout",
        html = m$message,
        add = TRUE)
      runjs("document.getElementById('stdout').parentElement.scrollTo({ top: 1e9, behavior: 'smooth' });")
    })
  })
  
}

shinyApp(ui, server)