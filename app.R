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
library(hover)

bin_on_path = function(bin) {
  exit_code = suppressWarnings(system2("command", args = c("-v", bin), stdout = FALSE))
  return(exit_code == 0)
}

log_run <- function(runid, nsamples, start_time, exit_status) {
  log_file <- "run_log.csv"
  log_entry <- data.frame(
    runid = runid,
    start_time = format(start_time),
    end_time = format(Sys.time(), tz = 'UTC', usetz = FALSE), # it is called when p finishes, tz!
    #datetime = as.character(Sys.time()),
    nsamples = nsamples,
    exit_status = exit_status,
    #ip = session$request$REMOTE_ADDR, # errors on the ws
    #agent = session$request$HTTP_USER_AGENT,
    stringsAsFactors = FALSE
  )
  if (!file.exists(log_file)) {
    write.table(log_entry, log_file, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
  } else {
    write.table(log_entry, log_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

sidebar <- sidebar(
  title = '',
  
  # controls
  shiny::div(
  id = 'controls',
  checkboxInput('demo', 'Use demo data', value = FALSE),
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
        "Upload reference file containing ONE sequence (fasta, genbank, embl or snapgene).",
        placement = "right")
    ), 
    multiple = F, accept = c('.fa', '.fasta', '.dna', '.gbk', '.genbank', '.embl'), placeholder = 'reference'),
  
  # if more than 1, the uploaded files will be placed in a tmp folder on the server and
  # then passed to the nxf script as dir
  fileInput(
      'upload_fastq', 
    tags$a(
      'Upload fastq',  
      tooltip(bsicons::bs_icon("question-circle"),
        "Upload either one or more fastq files to be mapped to the reference. Max file size allowed: 1 Gb",
        placement = "right")
    ), 
    multiple = T, accept = c('.fastq', '.gz', '.fq'), placeholder = 'fastq file(s)'),
  #hr(),
  hover_action_button('start', 'Run mapping', icon = icon('play'), button_animation = 'icon-fade'),
  hr(),
  hover_action_button('reset', 'Reset input', icon = icon('rotate'), button_animation = 'icon-fade'),
  hr(),
  checkboxInput('advanced', 'Advanced settings', value = F),
  conditionalPanel(
    condition = "input.advanced",
    selectInput(
      "profile", "Run profile",
      c("Singularity" = "singularity", "Docker" = "standard"), 
      selected = "singularity"
    )
    #textInput('minimap_params', 'Minimap parameters', ''),
    #checkboxInput('include_variants', 'Include variants in report', value = F)
  ),
  uiOutput("copy_error_btn")
  )
)

ui <- page_navbar(
  
  useShinyjs(),
  use_hover(),
  
  
  fillable = F,
  title = tags$span(
  tags$span(
    "nxf-minimapper app",
    style = "font-size: 1.3rem; font-weight: normal;"
  ),
  tags$span(
    icon('align-center'),
    "Nextflow pipeline to map long reads to a plasmid/amplicon reference.", icon('align-center'),
    style = "font-size: 0.9rem; font-weight: normal; margin-left: 3em;"
  )
  ),
  theme = bs_theme(bootswatch = 'yeti', primary = '#196F3D'),
  sidebar = sidebar,
  
  card(
    card_header(
      id = 'header1', 
      class = 'bg-secondary', 
      tags$a('Results', tooltip(bsicons::bs_icon("question-circle"), 'Results are not stored on the server, please download.'))
    ),
    height = 150,
    card_body(
      uiOutput('download_ui'),
      uiOutput('report_ui')
      #downloadLink('download', 'Download results')
    )
  ),
  
  card(
    card_header(
      id = 'header2', 
      class = 'bg-secondary', 
      tags$a('Output', tooltip(bsicons::bs_icon("question-circle"), 'Output from the nxf-minimapper pipeline'))
  ),
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
    
    # Remove whitespace and non-ASCII characters from reference name
    rname <- gsub("\\s+", "_", rname)                # Replace whitespace with underscore
    rname <- iconv(rname, to = "ASCII//TRANSLIT")    # Remove/convert non-ASCII
    rname <- gsub("[^A-Za-z0-9._-]", "", rname)      # Remove any remaining unwanted chars
    
    # Auto-select correct format
    
    # Check file extension matches selected format
    ext <- tools::file_ext(rname)
    # Auto-select correct format based on extension
    ext <- tolower(tools::file_ext(rname))

    valid_ext <- switch(
      input$format,
      fasta = c("fa", "fasta"),
      genbank = c("gbk", "genbank"),
      embl = c("embl"),
      snapgene = c("dna"),
      character(0)
    )
    if (!(tolower(ext) %in% valid_ext)) {
      notify_failure('Invalid reference file, please select correct format', position = 'center-bottom') 
      validate("Invalid reference file format.")
    }
    
    fs::file_move(rpath, fs::path(fs::path_dir(rpath), rname))
    fs::path(fs::path_dir(rpath), rname)
  })

  observeEvent(input$upload_ref, {
    req(input$upload_ref)
    rname <- input$upload_ref$name
    ext <- tolower(tools::file_ext(rname))
    auto_format <- dplyr::case_when(
      ext %in% c("fa", "fasta") ~ "fasta",
      ext %in% c("gbk", "genbank") ~ "genbank",
      ext %in% c("embl") ~ "embl",
      ext %in% c("dna") ~ "snapgene",
      TRUE ~ "fasta"
    )
    updateSelectInput(session, "format", selected = auto_format)
  })
  
  fastq <- reactive({
    fnames <- input$upload_fastq$name
    fpaths <- input$upload_fastq$datapath
    
    # Remove whitespace and non-ASCII characters from reference name
    fnames <- gsub("\\s+", "_", fnames)                # Replace whitespace with underscore
    fnames <- iconv(fnames, to = "ASCII//TRANSLIT")    # Remove/convert non-ASCII
    fnames <- gsub("[^A-Za-z0-9._-]", "", fnames)      # Remove any remaining unwanted chars
    fnames <- make.unique(fnames, sep = "_")           # make sure sample names are unique
    
    fs::file_move(fpaths, fs::path(fs::path_dir(fpaths), fnames))
    fs::path_dir(fpaths[1])
  })
  
  observeEvent(input$demo, {
    if (input$demo) {
      lapply(c('format', 'upload_ref', 'upload_fastq'), function(x) {shinyjs::hide(x)})
      
    } else {
      lapply(c('format', 'upload_ref', 'upload_fastq'), function(x) {shinyjs::show(x)})
    }
    
  })
  
  # progress bar
  progress_val <- reactiveVal(0)
  proc <- reactiveVal(NULL) # store process object per session
  
  observeEvent(input$start, {
    if ((is.null(input$upload_fastq) || is.null(input$upload_ref)) & !input$demo) {
      notify_failure('Please upload both reference and fastq files before starting.', position = 'center-bottom')
      return()
    }
    
    #req(input$upload_fastq)
    #req(input$upload_ref)

    nsamples <- ifelse(input$demo, 3, nrow(input$upload_fastq))
    #log_run(runid, nsamples) 
    inc <- 100/(2 + (6 * nsamples)) # 2 single proc and 6 proc that are per sample

     
      if (input$demo) {
        arguments <- c('--ref', 'www/demo/reference.gbk', '--fastq', 'www/demo/samples', '--format', 'genbank', '-ansi-log', 'false')
      } else {
        arguments <- c('--ref', ref(), '--fastq', fastq(), '--format', input$format, '-ansi-log', 'false')
      }
      
      
    shinyjs::disable('controls')
    lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-secondary')})
    lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-warning')})
    notify_info('Run started, please wait...', position = 'center-bottom')

    # Start Nextflow asynchronously
    p <- processx::process$new(
      'nextflow',
      args = c(
        'run', 
        '-w', paste0('work/', runid), # allows per session cleanup
        'angelovangel/nxf-minimapper', 
        '--outdir', file.path('www', runid), arguments, '-profile', input$profile),
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
              html = paste0("<br>", line),
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
            log_run(runid = runid, start_time = p$get_start_time(), nsamples = nsamples, exit_status = 0)
            notify_success('Processing finished, please download results!', position = 'center-bottom')
            shinyjs::enable('controls')
            lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-warning')})
            lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-success')})
            shinyjs::hide(id = 'pb')
            output$copy_error_btn <- renderUI({ NULL })
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
          } else {
            log_run(runid = runid, start_time = p$get_start_time(), nsamples = nsamples, exit_status = p$get_exit_status())
            notify_failure('Processing failed', position = 'center-bottom')
            shinyjs::enable('controls')
            lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-warning')})
            lapply(c('header1', 'header2'), function(x) {shinyjs::toggleClass(id = x, class = 'bg-danger')})
            shinyjs::hide(id = 'pb')

            output$copy_error_btn <- renderUI({
              hover_action_button("copy_error", "Copy error to clipboard", icon = icon("clipboard"), button_animation = 'icon-fade', style = "color: red;")
            })

            output$download_ui <- renderUI({ NULL })
            output$report_ui <- renderUI({ NULL })
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

  observeEvent(input$copy_error, {
  runjs(sprintf(
    "navigator.clipboard.writeText(document.getElementById('stdout').innerText);"
  ))
  notify_success("Error output copied to clipboard!", position = "center-bottom")
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