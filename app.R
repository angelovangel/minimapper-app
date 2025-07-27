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
  
  # controls
  shiny::div(
  id = 'controls',
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
  hr(),
  hover_action_button('reset', 'Reset inputs', icon = icon('rotate'), button_animation = 'icon-fade'),
  
)
)

ui <- page_navbar(
  # header = 
  #   tags$a('Live terminal view',
  #   tooltip(
  #   bsicons::bs_icon("question-circle"),
  #   "Preview of the terminal, for viewing the selected parameters and monitor output",
  #   placement = "right")),
  
  useShinyjs(),
  use_hover(),
  
  fillable = F,
  title = 'nxf-minimapper app',
  theme = bs_theme(bootswatch = 'yeti', primary = '#196F3D'),
  sidebar = sidebar,
  #nav_panel(
  
  card(
    card_header(id = 'header1', class = 'bg-secondary', 'Results'),
    height = 150,
    card_body(
      uiOutput('download_ui'),
      uiOutput('report_ui')
      #downloadLink('download', 'Download results')
    )
  ),
  
  card(
    card_header(id = 'header2', class = 'bg-secondary', 'Output'),
    height = 450,
    card_body(
      verbatimTextOutput('stdout')
    )
  ),

  tags$style(HTML("
  .progress-number { color: transparent !important; }
  .progress-number { display: none !important; }
  ")),
  progressBar(id = "pb", value = 0, total = 100, status = "warning", display_pct = FALSE, title = "")
)

server <- function(input, output, session) {
  
  ##
  options(shiny.maxRequestSize=1000*1024^2)
  ##
  
  ##
  runid <- digest::digest(runif(1), algo = 'crc32')
  ##
  
  # check nextflow and docker is on path
  if (!bin_on_path('nextflow') | !bin_on_path('docker')) {
    notify_failure('nextflow and/or docker not found!', position = 'center-bottom')
  } else {
    notify_success('The server is ready!', position = 'center-bottom')
    #show_alert('OK', 'The server is ready!', type = 'success', )
  }
  
  output$stdout <- renderText({
    paste0('Run id: ', runid, '\n')
  })
  
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
  
  
  # progress bar
  progress_val <- reactiveVal(0)
  proc <- reactiveVal(NULL) # store process object per session
  
  observeEvent(input$start, {
    req(input$upload_fastq)
    req(input$upload_ref)
    
    arguments <- c('--ref', ref(), '--fastq', fastq(), '--format', input$format, '-ansi-log', 'false')
    
    shinyjs::disable('controls')
    lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-secondary')})
    lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-warning')})
    
    nsamples <- nrow(input$upload_fastq)
    inc <- 100/(2 + (6 * nsamples)) # 2 single proc and 6 proc that are per sample
    
    # Start Nextflow asynchronously
    p <- processx::process$new(
      'nextflow',
      args = c(
        'run', 
        '-w', paste0('work/', runid), # allows per session cleanup
        'angelovangel/nxf-minimapper', 
        '--outdir', file.path('www', runid), arguments),
      stdout = "|", stderr = "2>&1"
      # env = c(NXF_OFFLINE = "TRUE") # uncomment if needed
    )
    proc(p) # store process object
    
    # Poll process output
    observe({
      invalidateLater(500, session)
      p <- proc()
      if (!is.null(p)) {
        if (p$is_alive()) {
          lines <- p$read_output_lines()
          for (line in lines) {
            message(line)
            if (grepl('Submitted process >', line)) {
              progress_val(progress_val() + inc)
              updateProgressBar(session, id = "pb", value = progress_val())
            }
            shinyjs::html(
              id = "stdout",
              html = paste0(gsub("\n", "<br>", line), "<br>"),
              add = TRUE
            )
            runjs("document.getElementById('stdout').parentElement.scrollTo({ top: 1e9, behavior: 'smooth' });")
          }
        } else {
          # Process finished
          lines <- p$read_all_output_lines()
          for (line in lines) {
            shinyjs::html(
              id = "stdout",
              html = paste0(gsub("\n", "<br>", line), "<br>"),
              add = TRUE
            )
          }
          if (p$get_exit_status() == 0) {
            notify_success('Processing finished', position = 'center-bottom')
            shinyjs::enable('controls')
            lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-warning')})
            lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-success')})
            shinyjs::hide(id = 'pb')
            output$download_ui <- renderUI({
              downloadLink('download', 'Download data')
            })
            output$report_ui <- renderUI({
              pathtoreport <- paste0(runid, '/00-alignment-summary.html')
              actionLink(
                'report', 'View HTML report',
                onclick = sprintf("window.open('%s', '_blank')", pathtoreport)
              )
            })
          }
          proc(NULL) # clear process object
          return(NULL)
        }
      }
    })
  })
  
  observeEvent(input$reset, {
    session$reload()
  })
  
  output$download <- downloadHandler(
  filename = function() {
    paste0('results-', runid, '.tar')
  },
  content = function(file) {
    run_dir <- file.path("www", runid)
    files_to_tar <- list.files(run_dir, full.names = FALSE) # relative names
    old_wd <- setwd(run_dir)
    on.exit(setwd(old_wd))
    utils::tar(
      tarfile = file,
      files = files_to_tar
    )
  }
)
  
  session$onSessionEnded(function() {
    # Clean up only this user's files, e.g. by runid
    user_www <- file.path("www", runid)
    user_work <- file.path("work", runid)
    if (dir.exists(user_www)) fs::dir_delete(user_www)
    if (dir.exists(user_work)) fs::dir_delete(user_work)
  })

}

# cleanup <- function() {
#   workfiles <- list.files(path = "work", full.names = T)
#   #wwwfiles <- list.files(path = 'www', full.names = T)
#   lapply(c(workfiles), fs::dir_delete)
# }
# onStop(function() { cleanup() })

shinyApp(ui, server)