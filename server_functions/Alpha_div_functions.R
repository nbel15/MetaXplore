# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# import requested functions

source("tools_functions.R")
source("analysis_functions/alpha_div.R")

#########################################
# load and change order of factors
observeEvent(input$metamat_input,
             updateSelectInput(session, "alpha.div.factor", choices = colnames(meta_file())
                               )
             )

#Change factor's order
observeEvent(input$metamat_input,
             output$alpha.div.order<- renderUI({
               if(input$alpha.div.change.order == T)
                 selectInput("alpha.div.new.order", 
                             choices = change.factor.order(physeq(), input$alpha.div.factor), 
                             label = NULL, multiple = T
                 )
             })
)

#Show show table button
observeEvent(input$alpha.div.go,
             output$alpha.div.show.tab<- renderUI({
               actionButton("alpha.div.show.tab.key", "Show table")
             })
)

# Parameter section
output$param.alpha.div.boxplot <- renderUI({
  bsCollapse(id = "collapse.param.alpha.boxplot", open = "",
             bsCollapsePanel("Visual parameters",
                             fluidRow(
                               column(width = 2, checkboxInput(inputId = "boxplot.visual.labels", label = "Labels", value = F)),
                               column(width = 2, checkboxInput(inputId = "boxplot.visual.color", label = "color", value = F))
                             ),
                             uiOutput("boxplot.labels.element"),
                             uiOutput("boxplot.color.element")
                             
             ),
             bsCollapsePanel("Download",
                             downloadButton("boxplot.download.png", "PNG"), 
                             downloadButton("boxplot.download.pdf", "PDF")
             )
  )
})

observeEvent(input$boxplot.visual.labels,
             output$boxplot.labels.element<- renderUI({
               if( input$boxplot.visual.labels==T){
                 
                 fluidRow(
                   column(width = 3,
                          sliderInput("boxplot.columns.label.size", "Labels size", min=10, max=24, value=12,step=2, width = "80%")
                          ),
                   column(width = 3,
                          sliderInput("boxplot.columns.label.rot", " Lables rotation", min=0, max=90, value=45,step=10, width = "80%")
                          ),
                   column(width = 3,
                          radioButtons("boxplot.label.hjust", label = "Ajust labels", inline = TRUE, 
                                       choices = list("Center" = 0.5, "Left" = 0, "Right" = 1), selected = 0.5, width = "80%")
                   )
                   )
                 }
               })
             )

observeEvent(input$boxplot.visual.color,
             output$boxplot.color.element<- renderUI({
               if( input$boxplot.visual.color==T){
                 selectInput("boxplot.color.list", label = "Boxplot color", selected = "A", width = "50%",
                             choices = c("magma" = "A", 
                                         "inferno" = "B",
                                         "plasma" = "C",
                                         "viridis" = "D",
                                         "cividis" = "E",
                                         "rocket" = "F",
                                         "mako" = "G", 
                                         "turbo" = "H") )
               }
             })
)
#########################################
# View Alpha Results

observe({
  if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input)) {
    shinyjs::disable("alpha.div.go")
  } else {
    shinyjs::enable("alpha.div.go")
  }
})

alpha.div.plot.fun <- function(){
  p <- alpha.div.plot(isolate(t(physeq())), 
                      method = isolate(input$alpha.div.indices), 
                      grouping_column = isolate(input$alpha.div.factor), 
                      pValueCutoff = isolate(input$alpha.pval.cutoff),
                      neword = isolate(input$alpha.div.new.order),
                      changefac = isolate(input$alpha.div.change.order))
  
  p <- p + theme(text = element_text(family="serif"),
            axis.text = element_text(size =input$boxplot.columns.label.size),
            axis.title = element_text(size =input$boxplot.columns.label.size),
            axis.text.x = element_text(angle = input$boxplot.columns.label.rot, hjust = input$boxplot.lable.hjust),
            legend.text = element_text(size =input$boxplot.columns.label.size),
            legend.title = element_text(size =input$boxplot.columns.label.size),
            strip.text.x = element_text(size = input$boxplot.columns.label.size, face = "bold"))
  
  if(!is.null(input$boxplot.color.list))
    p <- p + scale_colour_viridis(discrete = TRUE, option = input$boxplot.color.list, end = 0.8)
  
  return(p)
}

observeEvent(input$alpha.div.go,
             output$alpha.div.boxplot <-renderPlotly({
               ggplotly(alpha.div.plot.fun())
             })
)

#View alpha diversity indices table Expanded mod
observeEvent(input$alpha.div.show.tab.key, {
  showModal(modalDialog(
    title = "Alpha diversity indices",
    div(style="display: inline-block;vertical-align:top; width: 250px;",checkboxInput(inputId = "alpha.div.add.se", label = "Add standard error", value = F)),
    div(style="display: inline-block;vertical-align:top; width: 100px;",downloadButton("alpha.div.download.tab", "Download")),
    div(style = 'overflow-x: scroll',DT::dataTableOutput("alpha.div.tab.expand",width = "100%")),
    easyClose = TRUE,
    size = "l"
  ))
})


observeEvent(input$alpha.div.go, {
  output$alpha.div.tab.expand <- renderDataTable({
    dataset <- alpha.div.tab.se(physeq = isolate(physeq()), 
                                method = isolate(input$alpha.div.indices),
                                grouping_column = isolate(input$alpha.div.factor), 
                                add.se = input$alpha.div.add.se
                                )
    DT::datatable(dataset, rownames = T, escape = F, options = list(pageLength = 10, columnDefs = list(list(className = 'dt-center', targets = "_all"))))
  })
})


#########################################
#Download Plots

###Download plot PNG
output$boxplot.download.png <- downloadHandler(
  filename = function() { paste(input$alpha.div.factor,'-AlphaDiv-plot', '.png', sep='') },
  content = function(file) {
    if(!is.null(input$alpha.div.indices)) {
      if(length(input$alpha.div.indices) > 2) width = 1100
      else width = 700
    } 
   png(file, height = 500, width = width)
    print(alpha.div.plot.fun())
    dev.off()
  })

###Download plot PDF
output$boxplot.download.pdf <- downloadHandler(
  filename = function() { paste(input$alpha.div.factor,'-AlphaDiv-plot', '.pdf', sep='') },
  content = function(file) {
    if(!is.null(input$alpha.div.indices)) {
      if(length(input$alpha.div.indices) > 3) width = 1300
      else width = 700
    } 
    pdf(file, onefile=T,paper = "a4r", height = 500, width = width)
    print(alpha.div.plot.fun())
    dev.off()
  })


###Download alpha diversity indices table
output$alpha.div.download.tab <- downloadHandler(
  filename = function() { paste(input$alpha.div.factor,'-AlphaDiv-indices', '.txt', sep='') },
  content = function(file) {
    dataset <- alpha.div.tab.se(physeq = isolate(physeq()), 
                                method = isolate(input$alpha.div.indices),
                                grouping_column = isolate(input$alpha.div.factor), 
                                add.se = input$alpha.div.add.se
    )
    write.csv(x = dataset, file = file, quote = F, sep = "\t", row.names = T, col.names = T)
  })
