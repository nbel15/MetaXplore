# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# import requested functions

source("tools_functions.R")
source("analysis_functions/core_microbiome.R")

#########################################
# load and change order of factors
observeEvent(input$metamat_input,
             updateSelectInput(session, "core.factor", choices = colnames(meta_file())
             )
)

# Parameter section
output$param.core.table <- renderUI({
  bsCollapse(id = "collapse.param.alpha.boxplot", open = "",
             bsCollapsePanel("Download",
                             downloadButton("core.download.tab", "Download")
             )
  )
})

# Parameter core venn Diagram
output$core.venn.param <- renderUI({
  bsCollapse(id = "collapse.param.core.venn", open = "",
             bsCollapsePanel("Visual parameters",
                             fluidRow(
                               column(width = 5,
                                      sliderInput("core.label.size", "Labels size", min=4, max=10, value=6,step=1, width = "80%")),
                               column(width = 5,
                                      sliderInput("core.stroke.size", "Stroke size", min=0, max=1, value=1,step=0.1, width = "100%"))
                             )
             ),
             bsCollapsePanel("Download",
                             downloadButton("core.download.png", "PNG"), 
                             downloadButton("core.download.pdf", "PDF")
             )
  )
})

#Enable core go button when input files are uploaded
observe({
  if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input)) {
    shinyjs::disable("core.go")
  } else {
    shinyjs::enable("core.go")
  }
})

#Render core microbiome table
observeEvent(input$core.go, {
  output$core.table <- renderDataTable({
    dataset <- core.microbiome(physeq = isolate(physeq()),
                               core.factor = isolate(input$core.factor), 
                               prevalence = isolate(input$core.prevalence),
                               rb.threshold = isolate(input$core.rb.detection)
    )
    
    core.table <- data.frame(dataset$core.combined.taxa, check.names = F)
    DT::datatable(core.table, rownames = T, escape = F, options = list(pageLength = 10, columnDefs = list(list(className = 'dt-center', targets = "_all"))))
  })
})

core.venn.plot.fun <- function(){
  dataset <- core.microbiome(physeq = isolate(physeq()),
                             core.factor = isolate(input$core.factor), 
                             prevalence = isolate(input$core.prevalence),
                             rb.threshold = isolate(input$core.rb.detection)
  )
  if(is.null(input$core.stroke.size)) core.stroke.size = 1
  else core.stroke.size = input$core.stroke.size
  
  if(is.null(input$core.label.size)) core.label.size = 6
  else core.label.size = input$core.label.size
  
  ggvenn(data = dataset$list.core.taxa, 
         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#FFFFFF"),
         stroke_size = core.stroke.size, set_name_size = core.label.size, stroke_linetype = "solid"
  )
}
#Render core microbiome venn Diagram
observeEvent(input$core.go, {
  output$core.venn.plot <- renderPlot({
    showNotification("Initialize plot.", duration = 2, type = "message")
    core.venn.plot.fun()
  })
})


###Download plot Venn Diagram PNG
output$core.download.png <- downloadHandler(
  filename = function() { paste( input$core.factor,'-Core_venn','-plot', '.png', sep='') },
  content = function(file) {
    png(file, height = 700, width = 900, res = 100)
    print(core.venn.plot.fun())
    dev.off()
  })

###Download plot Venn Diagram PDF
output$core.download.pdf <- downloadHandler(
  filename = function() { paste(input$core.factor,'-Core_venn','-plot', '.pdf', sep='') },
  content = function(file) {
    pdf(file, onefile=F,paper = "a4r", height = 500, width = 700)
    print(core.venn.plot.fun())
    dev.off()
  })

##Download Core table
output$core.download.tab <- downloadHandler(
  filename = function() { paste(input$core.factor,'core_table', '.csv', sep='')},
  content = function(file) {
    dataset <- core.microbiome(physeq = isolate(physeq()),
                               core.factor = isolate(input$core.factor), 
                               prevalence = isolate(input$core.prevalence),
                               rb = isolate(input$core.rb.detection)
    )
    
    core.table <- data.frame(dataset$core.combined.taxa, check.names = F)
    write.csv(x = core.table, file = file, quote = F, sep = "\t", row.names = T)
  })