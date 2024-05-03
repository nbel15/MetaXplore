coretab <- tabItem(tabName = "core_tab",
                   box(
                     width = 12, title = "Core microbiome parameters", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                     fluidRow(
                       column(width = 4, 
                              selectInput("core.factor", "Factor", choices = c("choose.."))
                       ),
                       column(width = 4, 
                              numericInput("core.rb.detection", "Detection threshold (%)",value =  0.01, min = 0.01, max = 100, step = 0.01)
                       ),
                       column(width = 4, 
                              numericInput("core.prevalence", "Prevalence threshold (%)",value =  75, min = 10, max = 100, step = 5),
                              div(style="display: inline-block;vertical-align:top; width: 50px;", actionButton("core.go", "Go!"))
                       )
                     )
                   ),
                   fluidRow(
                     column(width = 7, offset = 0, style='padding:0px;',
                            box(width = 12, boxToolSize = "sm",
                                title = "Presence/Absence table", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                                fluidRow(
                                  column(12, align="center", 
                                         div(style = 'overflow-x: scroll',DT::dataTableOutput("core.table",width = "100%"))
                                  )
                                ),
                                uiOutput("param.core.table")
                            )),
                     column(width = 5, offset = 0, style='padding:0px;',
                            box(width = 12, boxToolSize = "sm",
                                title = "Venn Diagram", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                                fluidRow(
                                  column(12, 
                                         plotOutput("core.venn.plot"),
                                         uiOutput("core.venn.param")
                                  )
                                )
                            )
                     ))
                   
)