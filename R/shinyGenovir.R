#' @export
#'
#' @import shiny ggplot2 plotly shinyjs
shinyGenovir <- function() {
  shinyApp(ui = .app_ui, server = .app_server)
}


# Define the UI with custom font for the title
.app_ui <- navbarPage(
  title = div(
    style = "font-family: 'Dancing Script', cursive; font-weight: bold; color: #4CAF50; font-size: 36px;",
    "GenovisR"
  ),

  tabPanel("File Loading",
           fluidPage(
             useShinyjs(),
             tags$head(
               tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Dancing+Script:wght@700&display=swap")
             ),
             sidebarLayout(
               sidebarPanel(
                 numericInput("max_file_size", "Max File Size (MB)", value = 1024, min = 1, max = 2048, step = 1),
                 fileInput("file", "Choose GDS File", accept = ".gds"),
                 checkboxInput("load_filter", "Load Filter", value = TRUE)
               ),
               mainPanel(
                 h5("Please select the file."),
                 hr(),
                 h4("Select Parent Samples"),
                 actionButton("set_parents", "Set Parents"),
                 hr(),
                 h4("Rank Samples"),
                 fluidRow(
                   uiOutput("rank_samples_ui")
                 ),
                 hr(),
                 h4("Sample Selection"),
                 uiOutput("parent_toggle_buttons")
               )
             )
           )
  ),

  tabPanel("Plotting",
           fluidPage(
             sidebarLayout(
               sidebarPanel(
                 h4("Tweaking histograms"),
                 selectInput("data_type_hist", "Data Type for Histogram", choices = c("haplotype", "dosage")),
                 sliderInput("binwidth", "Binwidth", min = 1e5, max = 1e7, value = 1e6),
                 textInput("fill_color", "Fill Color", value = "darkgreen"),
                 checkboxInput("chrwise", "Chromosome-wise", value = FALSE),
                 checkboxInput("samplewise", "Sample-wise", value = FALSE),
                 numericInput("ncol_hist", "Number of Columns", value = 3, min = 1, step = 1),
                 actionButton("select_samples_hist", "Select Samples for Histogram"),
                 hr(),
                 h4("Graph Geno Parameters"),
                 selectInput("data_type_graph", "Data Type for Graph Geno", choices = c("haplotype", "dosage")),
                 selectInput("direction", "Direction", choices = c("h", "v"), selected = "h"),
                 numericInput("width", "Width", value = 0.1, min = 0.1, max = 0.9, step = 0.01),
                 actionButton("select_samples_graph", "Select Samples for Graph Geno"),
                 hr(),
                 h4("Stats Parameters"),
                 selectInput("data_type_stats", "Data Type for Stats", choices = c("haplotype", "dosage", "recomb")),
                 selectInput("group", "Group By", choices = c("hap", "sample", "chr", "genome")),
                 selectInput("value", "Value", choices = c("class", "segment"))
               ),
               mainPanel(
                 plotlyOutput("hist_plot"),
                 plotlyOutput("graph_geno_plot"),
                 plotlyOutput("stats_plot"),
                 textOutput("no_samples_selected")
               )
             )
           )
  )
)

# Define the server logic
.app_server <- function(input, output, session) {
  # Initialize prev_samples as a reactive value
  prev_samples <- reactiveVal(NULL)
  selected_samples_hist <- reactiveVal(NULL)
  selected_samples_graph <- reactiveVal(NULL)

  # Observe changes to the max_file_size input and update shiny.maxRequestSize accordingly
  observe({
    shiny.maxRequestSize <- input$max_file_size * 1024^2  # Convert MB to bytes
    options(shiny.maxRequestSize = shiny.maxRequestSize)
  })

  # Reactive expression to read the uploaded GDS file and extract necessary information
  gbsr_obj <- reactiveVal()
  genovis_obj <- reactiveVal()

  observeEvent(input$file, {
    req(input$file)
    tryCatch({
      gbsr <- loadGDS(input$file$datapath, load_filter = input$load_filter)
      gbsr_obj(gbsr)

      # Generate checkboxes for sample selection
      sample_names <- getSamID(gbsr, valid = FALSE)
      output$parent_toggle_buttons <- renderUI({
        fluidRow(
          column(12,
                 checkboxGroupInput(
                   "parent_samples", "Parent Samples", choices = sample_names, selected = NULL,
                   inline = TRUE,
                   width = '100%'
                 )
          )
        )
      })

      # Generate numeric inputs for ranking samples
      output$rank_samples_ui <- renderUI({
        req(input$parent_samples)
        sample_ranks <- seq_along(input$parent_samples)
        lapply(seq_along(input$parent_samples), function(i) {
          column(3, numericInput(paste0("rank_", input$parent_samples[i]), input$parent_samples[i], value = sample_ranks[i], min = 1, max = length(input$parent_samples), step = 1))
        })
      })

      # Show the 'Set Parents' button
      output$file_uploaded <- reactive(TRUE)
      outputOptions(output, 'file_uploaded', suspendWhenHidden = FALSE)
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        "Failed to load the GDS file. Please check the file format.",
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })

  observeEvent(input$set_parents, {
    req(gbsr_obj(), input$parent_samples)

    # Check if same numbers were assigned to multiple samples
    selected_samples <- input$parent_samples
    sample_ranks <- sapply(selected_samples, function(sample) {
      input[[paste0("rank_", sample)]]
    })

    if (any(duplicated(sample_ranks))) {
      showModal(modalDialog(
        title = "Error",
        "Duplicate ranks detected. Please assign unique ranks to each sample.",
        easyClose = TRUE,
        footer = NULL
      ))
      return()
    }

    # Confirm if users want to rerun setParents
    if (!is.null(prev_samples()) && !identical(prev_samples(), selected_samples)) {
      showModal(modalDialog(
        title = "Confirm",
        "Are you sure you want to rerun setParents?",
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_set_parents", "Yes")
        )
      ))
      return()
    }

    shinyjs::disable("set_parents")

    ranked_samples <- selected_samples[order(sample_ranks)]
    progress <- Progress$new(session, min = 0, max = 3)
    progress$set(message = "Processing...", value = 0)

    on.exit({
      progress$close()
      shinyjs::enable("set_parents")
    })

    progress$inc(1, detail = "Running setParents")
    gbsr <- setParents(gbsr_obj(), parents = ranked_samples)
    gbsr_obj(gbsr)
    prev_samples(selected_samples)

    progress$inc(1, detail = "Running gbscleanr2genovis")
    genovis <- gbscleanr2genovis(gbsr)

    progress$inc(1, detail = "Running evalSegments and statsGeno")
    genovis <- evalSegments(genovis)
    genovis <- statsGeno(genovis)
    genovis_obj(genovis)

    # Get sample names for default selections
    sample_names <- getSamID(gbsr)
    selected_samples_hist(sample_names)
    selected_samples_graph(sample_names[1:10])

    # Render initial plots
    output$hist_plot <- renderPlotly({
      req(selected_samples_hist())
      if (input$samplewise && length(selected_samples_hist()) > 10) {
        showModal(modalDialog(
          title = "Too many samples selected",
          "Please select 10 or fewer samples for sample-wise histogram.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else {
        ggplotly(plotHist(genovis_obj(), data = input$data_type_hist, sample = selected_samples_hist(), binwidth = input$binwidth, fill = input$fill_color, chrwise = input$chrwise, samplewise = input$samplewise, ncol = input$ncol_hist))
      }
    })

    output$graph_geno_plot <- renderPlotly({
      req(selected_samples_graph())
      if (length(selected_samples_graph()) > 10) {
        showModal(modalDialog(
          title = "Too many samples selected",
          "Please select 10 or fewer samples for Graph Geno.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else {
        ggplotly(plotGraphGeno(genovis_obj(), data = input$data_type_graph, sample = selected_samples_graph(), direction = input$direction, width = input$width))
      }
    })

    output$stats_plot <- renderPlotly({
      ggplotly(plotStats(genovis_obj(), data = input$data_type_stats, group = input$group, value = input$value))
    })

    output$no_samples_selected <- renderText({
      if (is.null(selected_samples_hist()) && is.null(selected_samples_graph())) {
        "Please select samples to be shown"
      } else {
        NULL
      }
    })
  })

  observeEvent(input$confirm_set_parents, {
    removeModal()
    selected_samples <- input$parent_samples
    sample_ranks <- sapply(selected_samples, function(sample) {
      input[[paste0("rank_", sample)]]
    })

    if (any(duplicated(sample_ranks))) {
      showModal(modalDialog(
        title = "Error",
        "Duplicate ranks detected. Please assign unique ranks to each sample.",
        easyClose = TRUE,
        footer = NULL
      ))
      return()
    }

    shinyjs::disable("set_parents")
    ranked_samples <- selected_samples[order(sample_ranks)]
    progress <- Progress$new(session, min = 0, max = 3)
    progress$set(message = "Processing...", value = 0)

    on.exit({
      progress$close()
      shinyjs::enable("set_parents")
    })

    progress$inc(1, detail = "Running setParents")
    gbsr <- setParents(gbsr_obj(), parents = ranked_samples)
    gbsr_obj(gbsr)
    prev_samples(selected_samples)

    progress$inc(1, detail = "Running gbscleanr2genovis")
    genovis <- gbscleanr2genovis(gbsr)

    progress$inc(1, detail = "Running evalSegments and statsGeno")
    genovis <- evalSegments(genovis)
    genovis <- statsGeno(genovis)
    genovis_obj(genovis)

    # Get sample names for default selections
    sample_names <- getSamID(gbsr)
    selected_samples_hist(sample_names)
    selected_samples_graph(sample_names[1:10])

    # Render initial plots
    output$hist_plot <- renderPlotly({
      req(selected_samples_hist())
      if (input$samplewise && length(selected_samples_hist()) > 10) {
        showModal(modalDialog(
          title = "Too many samples selected",
          "Please select 10 or fewer samples for sample-wise histogram.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else {
        ggplotly(plotHist(genovis_obj(), data = input$data_type_hist, sample = selected_samples_hist(), binwidth = input$binwidth, fill = input$fill_color, chrwise = input$chrwise, samplewise = input$samplewise, ncol = input$ncol_hist))
      }
    })

    output$graph_geno_plot <- renderPlotly({
      req(selected_samples_graph())
      if (length(selected_samples_graph()) > 10) {
        showModal(modalDialog(
          title = "Too many samples selected",
          "Please select 10 or fewer samples for Graph Geno.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else {
        ggplotly(plotGraphGeno(genovis_obj(), data = input$data_type_graph, sample = selected_samples_graph(), direction = input$direction, width = input$width))
      }
    })

    output$stats_plot <- renderPlotly({
      ggplotly(plotStats(genovis_obj(), data = input$data_type_stats, group = input$group, value = input$value))
    })

    output$no_samples_selected <- renderText({
      if (is.null(selected_samples_hist()) && is.null(selected_samples_graph())) {
        "Please select samples to be shown"
      } else {
        NULL
      }
    })
  })

  # Auto update plots when parameters are changed
  observe({
    req(genovis_obj())
    output$hist_plot <- renderPlotly({
      req(selected_samples_hist())
      if (input$samplewise && length(selected_samples_hist()) > 10) {
        showModal(modalDialog(
          title = "Too many samples selected",
          "Please select 10 or fewer samples for sample-wise histogram.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else {
        ggplotly(plotHist(genovis_obj(), data = input$data_type_hist, sample = selected_samples_hist(), binwidth = input$binwidth, fill = input$fill_color, chrwise = input$chrwise, samplewise = input$samplewise, ncol = input$ncol_hist))
      }
    })
  })

  observe({
    req(genovis_obj())
    output$graph_geno_plot <- renderPlotly({
      req(selected_samples_graph())
      if (length(selected_samples_graph()) > 10) {
        showModal(modalDialog(
          title = "Too many samples selected",
          "Please select 10 or fewer samples for Graph Geno.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else {
        ggplotly(plotGraphGeno(genovis_obj(), data = input$data_type_graph, sample = selected_samples_graph(), direction = input$direction, width = input$width))
      }
    })
  })

  observe({
    req(genovis_obj())
    output$stats_plot <- renderPlotly({
      ggplotly(plotStats(genovis_obj(), data = input$data_type_stats, group = input$group, value = input$value))
    })
  })

  observeEvent(input$select_samples_hist, {
    req(gbsr_obj())
    sample_names <- getSamID(gbsr_obj())
    showModal(modalDialog(
      title = "Select Samples for Histogram",
      checkboxGroupInput("hist_samples", "Samples", choices = sample_names, selected = selected_samples_hist(),
                         inline = FALSE, width = '100%'),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_hist_samples", "OK")
      )
    ))
  })

  observeEvent(input$confirm_hist_samples, {
    selected_samples_hist(input$hist_samples)
    removeModal()
  })

  observeEvent(input$select_samples_graph, {
    req(gbsr_obj())
    sample_names <- getSamID(gbsr_obj())
    showModal(modalDialog(
      title = "Select Samples for Graph Geno",
      checkboxGroupInput("graph_samples", "Samples", choices = sample_names, selected = selected_samples_graph(),
                         inline = FALSE, width = '100%'),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_graph_samples", "OK")
      )
    ))
  })

  observeEvent(input$confirm_graph_samples, {
    selected_samples_graph(input$graph_samples)
    removeModal()
  })
}
