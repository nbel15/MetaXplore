guide <-  tabItem(tabName = "guide",
                  fluidRow(
                    tabBox(
                      width = 12, 
                      tabPanel("MetaXplore Workflow", uiOutput("metaxplore_workflow"))
                    )) 
)