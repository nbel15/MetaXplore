betatab <- tabItem(tabName = "beta_div_graph",
                   box(
                     width = 12, title = "Beta diversity parameters", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                     fluidRow(
                       column(width = 3, 
                              selectInput("beta.div.factor", "Factor", choices = "choose")
                       ),
                       column(width = 3,
                              selectInput("dist.method", "Distance metric",
                                          choices = c("GUnifrac" = "d_0.5", "Weighted" = "d_1", "Unweighted" = "d_UW", 
                                                      "VAW Unifract" = "d_VAW", "Bray Curtis"="bray"), selected = "d_0.5")
                       ),
                       column(width = 3,
                              selectInput("ordination.method", "Ordination", 
                                          choices = c("MDS"="mds" , "NMDS"="nmds", 
                                                      "PCoA"="pcoa", "CAP"="cap"), selected = "pcoa"),
                              uiOutput("beta.div.cap.param")
                       ),
                       column(width = 3,
                              selectInput("permut.permanova", "Number of permutation", width = "80%", 
                                          choices = c(99, 999, 9999), selected = 999),
                              div(style="display: inline-block;vertical-align:top; width: 50px;", disabled(actionButton("beta.go", "Go!"))),
                              div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("beta.show.dist.tab"))
                       )
                     )
                   ),
                   fluidRow(
                     column(width = 7, offset = 0, style='padding:0px;',
                            box(width = 12, title = "Beta diversity graph", collapsible = T, collapsed = F, solidHeader = T, status = "primary",
                                plotOutput("beta.div.plot"),
                                uiOutput("beta.div.plot.param")
                            )
                     ),
                     column(width = 5, offset = 0, style='padding:0px;',
                            box(width = 12, 
                                title = "Pairwise Comparison", collapsible = T, collapsed = F, solidHeader = T, status = "primary",
                                plotOutput("beta.div.pair.plot"),
                                #div(style = 'overflow-x: scroll',DT::dataTableOutput("beta.dist.matrix.test",width = "100%")),
                                uiOutput("beta.div.stat1.param")
                            )
                     )
                   )
)