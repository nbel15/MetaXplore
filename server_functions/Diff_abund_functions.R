# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# import requested functions
source("analysis_functions/diff_adundance.R")

#########################################
# Activate Go button conditions
observe({
  if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input)) {# || is.null(input$tree_file_input)) {
    shinyjs::disable("diff.abund.go")
  } else {
        shinyjs::enable("diff.abund.go")
  }
})

#########################################
#Update UI
# Load list of factors for biomarker discovery
observeEvent(input$metamat_input,
             updateSelectInput(session, "diff.abund.factor", choices = colnames(phyloseq::sample_data(isolate(physeq()))) #colnames(meta_file())
             )
)

observe({
  if (!is.null(input$metamat_input) && !is.null(physeq()))
             updateSelectInput(session, "diff.abund.select.taxa", 
                               choices = c("all" = "all",sort(unique(phyloseq::tax_table(isolate(physeq()))[,input$diff.abund.level.taxa]))) #colnames(meta_file())
             )}
)

#########################################
# Show distance matrix
#Show show graph button
observeEvent(input$diff.abund.go,
             output$diff.abund.show.tab<- renderUI({
               actionButton("diff.abund.tab.key", "Show table")
             })
)
output$param.diff.abund <- renderUI({
  bsCollapse(id = "collapse.Download", open = "",
             bsCollapsePanel("Visual parameters",
                             fluidRow(
                               column(width = 3,
                                      sliderInput("diff.abund.columns.label.size", "Labels size", min=10, max=24, value=14,step=2, width = "80%")
                               ),
                               column(width = 3,
                                      sliderInput("diff.abund.columns.label.rot", " Column lable rotation", min=0, max=90, value=90,step=10, width = "80%")
                               ),
                               column(width = 3,
                                      radioButtons("diff.abund.label.hjust", label = "Ajust labels", inline = TRUE, 
                                                   choices = list("Center" = 0.5, "Left" = 0, "Right" = 1), selected = 0.5, width = "80%")
                               )
                             )
             ),
             bsCollapsePanel("Download",
                             fluidRow(
                               column(width = 3,
                                      sliderInput("diff.abund.width", "width", min=10, max=60, value=30, step=1, width = "100%")
                               ),
                               column(width = 3,
                                      sliderInput("diff.abund.height", "Height", min=10, max=50, value=25, step=1, width = "100%"),
                               )
                             ),
                             downloadButton("diff.abund.download.png", "PNG"), 
                             downloadButton("diff.abund.download.pdf", "PDF")
                             )
  )
})

### Diff. abundance plots
diff.abund.plot.out <- function(){
  showNotification("Initialize the analysis", duration = 3, type = "message")
  t2 <- diff.abund.analysis(physeq = isolate(physeq()),
                            group.fac = isolate(input$diff.abund.factor), 
                            level.taxa = isolate(input$diff.abund.level.taxa),
                            abund.threshold = isolate(input$rabund.threshold))
  
  if (nrow(t2$res_diff) != 0) {
    if(isolate(input$diff.abund.select.taxa) %in% c("all", unique(t2$res_diff$Taxa))){
      diff.abund.plot(t2 = t2, select_taxa = isolate(input$diff.abund.select.taxa),
                       angle = as.numeric(input$diff.abund.columns.label.rot), 
                       ajust = as.numeric(input$diff.abund.label.hjust), 
                       label.size = as.numeric(input$diff.abund.columns.label.size)
                       )
    }else {
      shinyalert("Oops!", "No significant taxa can be used to plot the abundance!", type = "warning")
      return(NULL)
    }
  }else {
    shinyalert("Oops!", "No significant taxa can be used to plot the abundance!", type = "warning")
    return(NULL)
  }
} 

observeEvent(input$diff.abund.go,
             output$diff.abund.plot<- renderPlot({
               diff.abund.plot.out()
             })
)

#########################################
# View Results
# Differential table Expanded mod
observeEvent(input$diff.abund.tab.key, {
  showModal(modalDialog(
    title = paste("Differential abundance table : ", input$diff.abund.factor, sep = ""),
    div(style="display: inline-block;vertical-align:top; width: 100px;",downloadButton("diff.abund.tab", "Download")),
    div(style = 'overflow-x: scroll',DT::dataTableOutput("diff.abund.tab.expand",width = "100%")),
    easyClose = TRUE,
    size = "l"
  ))
})

##Differential abundance
observeEvent(input$diff.abund.go,
             output$diff.abund.tab.expand<- renderDataTable({
               diff.abund.table.export(physeq = isolate(physeq()),
                                       group.fac = isolate(input$diff.abund.factor), 
                                       level.taxa = isolate(input$diff.abund.level.taxa),
                                       select_taxa = isolate(input$diff.abund.select.taxa),
                                       abund.threshold =  isolate(input$rabund.threshold))
             })
)

#### Download
  ##Differential abundance table
  output$diff.abund.tab <- downloadHandler(
    filename = function() { paste(input$diff.abund.factor,input$diff.abund.select.taxa,'Diff_abund_tab.csv', sep='_')},
    content = function(file) {
      diff.tab <- diff.abund.table.export(physeq = isolate(physeq()),
                                          group.fac = isolate(input$diff.abund.factor), 
                                          level.taxa = isolate(input$diff.abund.level.taxa),
                                          select_taxa = isolate(input$diff.abund.select.taxa) )
      
      write.csv(x = diff.tab, file = file, quote = F, sep = "\t", row.names = F)
    })

  ##Differential abundance graph
    ###Download plot PNG
      output$diff.abund.download.png <- downloadHandler(
        filename = function() { paste(input$diff.abund.factor,input$diff.abund.select.taxa,'Diff_abund_plot.png', sep='_') },
        content = function(file) {
          #png(file, height = 500, width = 700, res = 100)
          width = conv_unit(isolate(input$diff.abund.width), "cm", "inch")*1.25
          height = conv_unit(isolate(input$diff.abund.height), "cm", "inch")*1.25
          
          grDevices::png(file, units = "px", height = height*96, width = width*96) #, res = 150
          print(diff.abund.plot.out())
          dev.off()
        })

    ###Download plot PDF
      output$diff.abund.download.pdf <- downloadHandler(
        filename = function() { paste(input$diff.abund.factor,input$diff.abund.select.taxa,'Diff_abund_plot.pdf', sep='_') },
        content = function(file) {
          #pdf(file, onefile=F,paper = "a4r", height = 500, width = 700)
          width = conv_unit(isolate(input$diff.abund.width), "cm", "inch") #15
          height = conv_unit(isolate(input$diff.abund.height), "cm", "inch") #10
          
          pdf(file, onefile = T,  height = height*1.25, width = width*1.25)
          print(diff.abund.plot.out())
          dev.off()
        })