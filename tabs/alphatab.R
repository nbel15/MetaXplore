alphatab <- tabItem(tabName = "alpha_tab",
                    box(
                      width = 12, title = "Alpha diversity parameters", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                      fluidRow(
                        column(width = 4, 
                               selectInput("alpha.div.factor", "Factor", choices = c("choose..")),
                               checkboxInput("alpha.div.change.order", label = "Change X-axis order", value = F),
                               uiOutput("alpha.div.order")
                        ),
                        column(width = 4, 
                               selectInput("alpha.div.indices", "Alpha diversity indice(s)", 
                                           choices = c(#"Goods coverage" = "goods", #"Richness" = "richness",
                                             "Observed OTUs" = "obs.otus",  
                                             "ACE" = "ace", "Chao" = "chao", 
                                             "Simpson" = "simpson", "Shannon" = "shannon", 
                                             "Evenness" = "evenness", "PD"= "pd",
                                             "Simpson effecive" = "simpson_eff", "Shannon effective" = "shannon_eff"), 
                                           selected = "simpson", multiple=TRUE, selectize=TRUE)
                        ),
                        column(width = 4, 
                               numericInput("alpha.pval.cutoff", "P-value cutoff",value =  0.05, min = 0.02, max = 1, step = 0.01),
                               div(style="display: inline-block;vertical-align:top; width: 50px;", actionButton("alpha.div.go", "Go!")),
                               div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("alpha.div.show.tab"))
                        )
                      )
                    ),
                    box(width = 12, boxToolSize = "sm",
                        title = "Boxplot", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                        fluidRow(
                          column(12, align="center", 
                                 plotlyOutput("alpha.div.boxplot", width = "100%", height = "500px")
                          )
                        ),
                        uiOutput("param.alpha.div.boxplot")
                        
                    )
)