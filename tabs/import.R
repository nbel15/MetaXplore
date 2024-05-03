import <- tabItem(tabName = "import_data", 
                  fluidRow(
                    box(
                      width = 12,
                      title = "Load Data", status = "primary", solidHeader = TRUE,
                      collapsible = TRUE,
                      column(width = 3, 
                             fileInput('otumat_input', 'Load OTUs file',
                                       accept=c('text/tsv', 
                                                'text/comma-separated-values,text/plain', '.tsv')),
                             selectInput("normalize", "Select normalization", 
                                         choices = c("None"="None", "Rarefy"="Rarefy", 
                                                     "Total sum scaling (TSS)" = "TSS",
                                                     "Centered log-ratio normalization (CLR)" = "CLR",
                                                     "log(X+1)" = "log",
                                                     "Square root" = "sqrt",
                                                     "OTUs Norm"= "otunorm"),
                                         selected =  "None", multiple=F, selectize=TRUE)
                      ),
                      column(width = 3, 
                             fileInput('taxamat_input', 'Load Taxa file',
                                       accept=c('text/tsv',
                                                'text/comma-separated-values,text/plain', '.tsv')),
                             numericInput("min.freq", "Minimum frequency (%)",min = 0, max = 90, value = 0, step = 0.001)
                             ),
                      column(width = 3, 
                             fileInput('metamat_input', 'Load Metadata file',
                                       accept=c('text/tsv',
                                                'text/comma-separated-values,text/plain', '.tsv'))
                      ),
                      column(width = 3, 
                             fileInput('tree_file_input', 'Load tree file (optional)',
                                       accept=c(".nwk",'.NWK','.tre', '.tree')),
                             div(class = "pull-right", actionButton("read_input", "Read Input", icon = icon("circle-down")))
                             
                      )
                    )
                  ),
                  fluidRow(
                    tags$head(tags$style(HTML('.info-box {min-height: 50px;} .info-box-icon {height: 50px; line-height: 50px;} .info-box-content {padding-top: 0px; padding-bottom: 0px;}'))),
                    infoBoxOutput("nsamples"),
                    infoBoxOutput("ntaxa"),
                    infoBoxOutput("info_input")
                  ),
                  fluidRow(
                    tabBox(width = 12, id = "tabview",
                           title = "Visualize Input",
                           tabPanel("OTUs table",
                                    div(style = 'overflow-x: scroll', DT::dataTableOutput("otu_tab_view",width = "100%"))),
                           tabPanel("Taxonomy table", 
                                    div(style = 'overflow-x: scroll', DT::dataTableOutput("taxa_tab_view",width = "100%"))),
                           tabPanel("Metatable", 
                                    div(style = 'overflow-x: scroll', DT::dataTableOutput("meta_tab_view",width = "100%")))
                    )
                  )
)