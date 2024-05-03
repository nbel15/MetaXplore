diffabundtab <- tabItem(tabName = "diff_abund_tab",
                   box(
                     width = 12, title = "Differential abundance parameters", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                     fluidRow(
                       column(width = 3, 
                              selectInput("diff.abund.factor", "Factor", choices = "choose")
                       ),
                       column(width = 3,
                              selectInput("diff.abund.level.taxa", "Taxa Level", 
                                          choices = c("Phylum"="Phylum", "Class"="Class", "Order"= "Order", "Family"="Family", "Genus"="Genus"),
                                          selected =  "Genus", multiple=F, selectize=TRUE)
                       ),
                       column(width = 3,
                              sliderInput("rabund.threshold", "Rb. Threshold (%)",min = 0, max = 10, value = 0, step = 0.1)
                       ),
                       column(width = 3,
                              selectInput("diff.abund.select.taxa", "Select Taxon", choices = "all"),
                              div(style="display: inline-block;vertical-align:top; width: 50px;", disabled(actionButton("diff.abund.go", "Go!"))),
                              div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("diff.abund.show.tab"))
                       )
                     )
                   ),
                   box(width = 12, #height = "800px", 
                       title = "Differential abundance", collapsible = T, collapsed = F, solidHeader = T, status = "primary",
                       plotOutput("diff.abund.plot", height = "600px"),
                       uiOutput("param.diff.abund")

                   )
                   )
