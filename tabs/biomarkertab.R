biomarkertab <- tabItem(tabName = "biomarker_tab",
                   box(
                     width = 12, title = "Biomarker discovery parameters", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                     fluidRow(
                       column(width = 3, 
                              selectInput("biomarker.factor", "Factor", choices = "choose")
                       ),
                       column(width = 3,
                              selectInput("biomarker.level.taxa", "Taxa Level", 
                                          choices = c("Phylum"="Phylum", "Class"="Class", "Order"= "Order", "Family"="Family", "Genus"="Genus", "OTU"=NULL),
                                          selected =  "Genus", multiple=F, selectize=TRUE)
                       ),
                       column(width = 3,
                              sliderInput("biomarker.lda.threshold", "LDA threshold",min = 0, max = 5, value = 2, step = 0.5)

                       ),
                       column(width = 3,
                              checkboxInput("biomarker.strict", "Strict (one-against-one)", value = TRUE, width = NULL),
                              checkboxInput("biomarker.rb.factor", "Relative abundance by sample", value = FALSE, width = NULL),
                              div(style="display: inline-block;vertical-align:top; width: 50px;", disabled(actionButton("biomarker.go", "Go!")))
                       )
                     )
                   ),
                   box(width = 12, #height = "800px", 
                       title = "Discriminant  taxa ", collapsible = T, collapsed = F, solidHeader = T, status = "primary",
                       plotOutput("biomarker.plot",height = "1000px"),
                       uiOutput("biomarker.plot.download")
                       #div(style = 'overflow-x: scroll',DT::dataTableOutput("diff.biomarker.tab",width = "100%"))
                       
                   )
                   )
