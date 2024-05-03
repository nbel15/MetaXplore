rbtab <- tabItem(tabName = "rb_tab", 
                 
                 fluidRow(
                   #tags$head(tags$style(HTML('.box {margin: 1px;}'))),
                   column(width = 12, offset = 0, style='padding:0px;',
                          box(
                            width = 12, title = "Stacked barchart", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                            fluidRow(
                              column(width = 4,
                                     selectInput("st.barchart.factor", "Factor", choices = "choose"),
                                     checkboxInput("st.barchart.change.order", label = "Change X-axis order", value = F),
                                     uiOutput("st.barchart.order")
                              ),
                              column(width = 4,
                                     selectInput("st.barchart.level.taxa", "Taxa Level(s)", 
                                                 choices = c("Kingdom"="Kingdom", "Phylum"="Phylum", "Class"="Class", "Order"= "Order", "Family"="Family", "Genus"="Genus", "Species"="Species"),
                                                 selected =  "Class", multiple=TRUE, selectize=TRUE)
                              ), 
                              column(width = 4,
                                     sliderInput("st.barchart.threshold", "Threshold (%)",min = 0, max = 10, value = 1, step = 0.5),
                                     div(style="display: inline-block;vertical-align:top; width: 50px;", actionButton("st.barchart.go", "Go!")),
                                     div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("st.barchart.show.tab"))
                              )
                            ),
                            plotlyOutput("st.barchart", height = "500px"),
                            uiOutput("param.st.barchart")
                          )
                   )),
                 fluidRow(
                   column(width = 12, offset = 0, style='padding:0px;',
                          box(width = 12,
                              title = "Heatmap", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                              fluidRow(
                                column(width = 4,
                                       selectInput("heatmap.factor", "Factor", choices = "choose"),
                                       checkboxInput("heatmap.change.order", label = "Change X-axis order", value = F),
                                       uiOutput("heatmap.order")
                                ),
                                column(width = 4,
                                       selectInput("heatmap.level.taxa", "Taxa Level(s)", 
                                                   choices = c("Kingdom"="Kingdom", "Phylum"="Phylum", "Class"="Class", "Order"= "Order", "Family"="Family", "Genus"="Genus", "Species"="Species"),
                                                   selected =  "Genus", multiple=TRUE, selectize=TRUE),
                                       checkboxInput("heatmap.rot", label = "Rotate heatmap", value = F)
                                ),
                                column(width = 4,
                                       sliderInput("heatmap.threshold", "Threshold (%)",min = 0, max = 10, value = 1, step = 0.1),
                                       div(style="display: inline-block;vertical-align:top; width: 50px;", actionButton("heatmap.go", "Go!")),
                                       div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("heatmap.show.tab"))
                                )
                              ),
                              plotOutput("heatmap", brush = "ht_brush", click = "ht_click", height = "500px"),
                              uiOutput("param.heatmap")
                          )
                   )
                 ),
                 fluidRow(
                   column(width = 12,  offset = 0, style='padding:0px;',
                          box(width = 12, 
                              title = "Barchart", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                              fluidRow(
                                column(width = 4,
                                       selectInput("barchart.factor", "Factor", choices = "choose"),
                                       checkboxInput("barchart.change.order", label = "Change X-axis order", value = F),
                                       uiOutput("barchart.order")
                                ),
                                column(width = 4,
                                       selectInput("barchart.level.taxa", "Taxa Level(s)", 
                                                   choices = c("Kingdom"="Kingdom", "Phylum"="Phylum", "Class"="Class", "Order"= "Order", "Family"="Family", "Genus"="Genus"),
                                                   selected =  "Genus", multiple=TRUE, selectize=TRUE)
                                ),
                                column(width = 4,
                                       sliderInput("barchart.num.OTUs", "Number of OTUs",min = 0, max = 12, value = 3, step = 1),
                                       div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton("barchart.go", "Go!")),
                                       div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("barchart.show.tab"))
                                       
                                )
                              ),
                              plotlyOutput("barchart"),
                              uiOutput("param.barchart")
                          )
                   )
                 )
)