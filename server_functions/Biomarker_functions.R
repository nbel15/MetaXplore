# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# import requested functions

#source("tools_functions.R")
source("analysis_functions/biomarker.R")

#########################################
# Activate Go button conditions

observe({
  if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input)) {# || is.null(input$tree_file_input)) {
    shinyjs::disable("biomarker.go")
  } else {
    shinyjs::enable("biomarker.go")
  }
})

#########################################
#Update UI
# Load list of factors for biomarker discovery
observeEvent(input$metamat_input,
             updateSelectInput(session, "biomarker.factor", choices = colnames(sample_data(isolate(physeq()))) #colnames(meta_file())
             )
)

output$biomarker.plot.download <- renderUI({
  bsCollapse(id = "collapse.Download", open = "",
             bsCollapsePanel("Download",
                             fluidRow(
                               column(width = 3,
                                      sliderInput("biomarker.width", "width", min=10, max=60, value=30, step=1, width = "100%")
                               ),
                               column(width = 3,
                                      sliderInput("biomarker.height", "Height", min=10, max=50, value=25, step=1, width = "100%"),
                               )
                             ),
                             downloadButton("biomarker.download_png", "PNG"), 
                             downloadButton("biomarker.download_pdf", "PDF")
                             #downloadButton("heatmap.download_svg", "SVG")
             )
  )
})
###Biomarker plots
biomarker.plot.out <- function(){
  showNotification("Initialize the biomarker analysis", duration = 3, type = "message")
  biomarker.plot.fun(physeq = isolate(physeq()), 
                     group.fac = isolate(input$biomarker.factor), 
                     rb.factor = isolate(input$biomarker.rb.factor), 
                     lda.threshold = isolate(input$biomarker.lda.threshold), 
                     level.taxa = isolate(input$biomarker.level.taxa),
                     strict =  isolate(input$biomarker.strict))
  
} 

observeEvent(input$biomarker.go,
             output$biomarker.plot<- renderPlot({
               biomarker.plot.out()
             })
)

#### Download
output$biomarker.download_pdf <- downloadHandler(
  filename = function() { paste('biomarker_', isolate(input$biomarker.factor), '.pdf', sep='') },
  content = function(file) {
    biomarker.width = conv_unit(isolate(input$biomarker.width), "cm", "inch") #15
    biomarker.height = conv_unit(isolate(input$biomarker.height), "cm", "inch") #10
    
    pdf(file, onefile = T,  height = biomarker.height*1.25, width = biomarker.width*1.25)
    print(biomarker.plot.out())
    dev.off()
    
  }
)

output$biomarker.download_png <- downloadHandler(
  filename = function() { paste('biomarker_', isolate(input$biomarker.factor), '.png', sep='') },
  content = function(file) {
    #biomarker.width = 15
    #biomarker.height = 10
    
    biomarker.width = conv_unit(isolate(input$biomarker.width), "cm", "inch")*1.25
    biomarker.height = conv_unit(isolate(input$biomarker.height), "cm", "inch")*1.25
    
    grDevices::png(file, units = "px", height = biomarker.height*96, width = biomarker.width*96) #, res = 150

    #grDevices::png(file, units = "px", height = 1000, width = 800)
    print(biomarker.plot.out())
    dev.off()
    
  }
)