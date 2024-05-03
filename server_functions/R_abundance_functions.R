# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# import requested functions

source("tools_functions.R")
source("analysis_functions/r_abundance.R")

#########################################
# Load list of factors

observeEvent(input$metamat_input,
             updateSelectInput(session, "abun.fac", choices = colnames(meta_file())
             )
)

#########################################
# Activate Go button conditions

observe({
  if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input) || is.null(input$tree_file_input)) {
    shinyjs::disable("goAbun")
  } else {
    shinyjs::enable("goAbun")
  }
})

######################################################
#################  Stacked barchart  #################
######################################################
observe({
  if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input) || is.null(input$tree_file_input)) {
    shinyjs::disable("st.barchart.avg")
  } else {
    shinyjs::enable("st.barchart.avg")
  }
})

output$param.st.barchart <- renderUI({
  bsCollapse(id = "collapse.param.st.barchart", open = "",
             bsCollapsePanel("Visual parameters",
                             #checkboxGroupInput(inputId = "st.barchart.visual", label = NULL, inline = T,
                              #                  choices = c("Labels" = "Labels", "Color" = "Color")
                               #                 ),
                             fluidRow(
                               column(width = 2, checkboxInput(inputId = "st.barchart.visual.labels", label = "Labels", value = F)),
                               column(width = 2, checkboxInput(inputId = "st.barchart.visual.color", label = "Color", value = F))
                             ),
                             
                             uiOutput("st.barchart.labels.element"),
                             uiOutput("st.barchart.color.element")
                        
             ),
             bsCollapsePanel("Download",
                             fluidRow(
                               column(width = 4,
                                      selectInput("st.barchart.plot.format", "Format", choices = c("PNG" = "png", "PDF" = "pdf"), width = "80%")
                               ),
                               column(width = 4,
                                      sliderInput("st.barchart.download.height", "Height", min=500, max=1200, value=500,step=20, width = "80%")
                               ),
                               column(width = 4,
                                      sliderInput("st.barchart.download.width", "Width", min=500, max=1200, value=500,step=20, width = "80%")
                               )
                             ),
                             downloadButton("st.barplot.download.plot", "")
                             #downloadButton("st.barplot.download_png", "PNG"), 
                             #downloadButton("st.barplot.download_pdf", "PDF")
             )
  )
})

observeEvent(input$st.barchart.visual.labels,
             output$st.barchart.labels.element<- renderUI({
               if( input$st.barchart.visual.labels==T){
                 
                 fluidRow(
                   column(width = 4,
                          sliderInput("st.barchart.columns.label.size", "Labels size", min=10, max=24, value=12,step=2, width = "80%")
                        ),
                   column(width = 4,
                          sliderInput("st.barchart.columns.label.rot", " Column lable rotation", min=0, max=90, value=45,step=10, width = "80%")
                        ),
                   column(width = 4,
                          textInput("st.barplot.title", "Plot title", "Enter title...")
                        )
                 )
                 }
               })
             )

observeEvent(input$st.barchart.visual.color,
             output$st.barchart.color.element<- renderUI({
               if(input$st.barchart.visual.color == T){
                
                 selectInput("st.barchart.color.list", label = "Barchart color", selected = "A", width = "50%",
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

observeEvent(input$metamat_input,
             updateSelectInput(session, "st.barchart.factor", choices = colnames(meta_file())
             )
)

#Change factor's order
observeEvent(input$metamat_input,
             output$st.barchart.order<- renderUI({
               if(input$st.barchart.change.order == T)
                 selectInput("st.barchart.new.order", 
                             choices = change.factor.order(physeq(), input$st.barchart.factor), 
                             label = NULL, multiple = T
                 )
             })
)

#Show show table button
observeEvent(input$st.barchart.go,
             output$st.barchart.show.tab<- renderUI({
               actionButton("st.barchart.show.tab.key", "Show table")
             })
)

st.barchart.plot <- function(){
  showNotification("Initialize the stucked Barchart.", duration = 2, type = "message")
  
  dataset <- data.plot(isolate(physeq()),
                       group_name = isolate(input$st.barchart.factor),
                       change_fac = isolate(input$st.barchart.change.order),
                       ordered_fac = isolate(input$st.barchart.new.order),
                       threshold = isolate(input$st.barchart.threshold),
                       level_taxa = isolate(input$st.barchart.level.taxa))
  
  
  #plot_abd(average5[,1:groupsnum], upper, lower, isolate(input$OTUs))
  p1 <- st.barchart.plot.main(dataset$AVG_tab)  
  if(is.null(input$st.barchart.color.list)) viridiscolor = "H"
  else viridiscolor = input$st.barchart.color.list
  p1 + 
    ggtitle(input$st.barplot.title) + 
    scale_fill_viridis(discrete = TRUE, option = viridiscolor, direction = -1)+
    theme(axis.text = element_text(size = input$st.barchart.columns.label.size),
          axis.title = element_text(size = input$st.barchart.columns.label.size),
          axis.text.x = element_text(angle = input$st.barchart.columns.label.rot),
          legend.text = element_text(size = input$st.barchart.columns.label.size),
         # legend.title = element_text(size = input$st.barchart.columns.label.size),
          text = element_text(family="serif"), 
          plot.title = element_text(size = input$st.barchart.columns.label.size, hjust = 0.5),
          legend.position = "right", legend.title = element_blank()
          )
}
#View Relative abundance Plot
observeEvent(input$st.barchart.go, 
             {
               output$st.barchart <- renderPlotly({
                 ggplotly(st.barchart.plot())
                 })
             })

#View relative abundance table --Heatmap-- Expanded mod
observeEvent(input$st.barchart.show.tab.key, {
  showModal(modalDialog(
    title = "Relative abundance table",
    div(style="display: inline-block;vertical-align:top; width: 250px;",checkboxInput(inputId = "st.barchart.add.se", label = "Add standard error", value = F)),
    div(style="display: inline-block;vertical-align:top; width: 100px;",downloadButton("st.barchart.download.rb.tab", "Download")),
    div(style = 'overflow-x: scroll',DT::dataTableOutput("st.barchart.avg.tab.expand",width = "100%")),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$st.barchart.go, {
  output$st.barchart.avg.tab.expand <- renderDataTable({
    dataset <- data.plot(isolate(physeq()),
                         group_name = isolate(input$st.barchart.factor),
                         change_fac = isolate(input$st.barchart.change.order),
                         ordered_fac = isolate(input$st.barchart.new.order),
                         threshold = isolate(input$st.barchart.threshold),
                         level_taxa = isolate(input$st.barchart.level.taxa))
    AVG_tab <- dataset$AVG_tab
    AVG_tab <- round(AVG_tab[2:(dim(AVG_tab)[2]-1)], digits = 3)
    SE_tab <- dataset$SE_tab
    SE_tab <- round(SE_tab[2:dim(SE_tab)[2]], digits = 2)
    merged.tab <- merge.tab.se.fun(tab = AVG_tab, se.tab = SE_tab, add.se = input$st.barchart.add.se)
    DT::datatable(merged.tab, rownames = T, escape = F, options = list(pageLength = 10))
  })
})

#### Download
output$st.barplot.download.plot <- downloadHandler(
  filename = function() { paste('st.barchart', isolate(input$st.barchart.factor), '.', isolate(input$st.barchart.plot.format) , sep='') },
  content = function(file) {
    switch(input$st.barchart.plot.format,
           "png" = png(file, height = input$st.barchart.download.height, width = input$st.barchart.download.width),
           "pdf" = pdf(file, onefile=T, height = input$st.barchart.download.height, width = input$st.barchart.download.width),
           "svg" = svg(file, height = input$st.barchart.download.height, width = input$st.barchart.download.width)
           )
      print(st.barchart.plot())
    dev.off()
  }
)

output$st.barchart.download.rb.tab <- downloadHandler(
  filename = function() { paste('st.barchart_Average_table', isolate(input$st.barchart.factor), '.csv', sep='') },
  content = function(file) {
    dataset <- data.plot(isolate(physeq()),
                         group_name = isolate(input$st.barchart.factor),
                         change_fac = isolate(input$st.barchart.change.order),
                         ordered_fac = isolate(input$st.barchart.new.order),
                         threshold = isolate(input$st.barchart.threshold),
                         level_taxa = isolate(input$st.barchart.level.taxa))
    AVG_tab <- dataset$AVG_tab
    AVG_tab <- round(AVG_tab[2:(dim(AVG_tab)[2]-1)], digits = 3)
    SE_tab <- dataset$SE_tab
    SE_tab <- round(SE_tab[2:dim(SE_tab)[2]], digits = 2)
    save.rb.table <- merge.tab.se.fun(tab = AVG_tab, se.tab = SE_tab, add.se = input$st.barchart.add.se)
    
    write.csv(x = save.rb.table, file = file, quote = F, sep = "\t", row.names = T, col.names = T)
  }
)

######################################################
#################  Heatmap  #################
######################################################
output$param.heatmap <- renderUI({
  bsCollapse(id = "collapse.param.heatmap", open = "",
             bsCollapsePanel("Visual parameters",
                             #checkboxGroupInput(inputId = "heatmap.visual", label = NULL, inline = T, 
                              #                  choices = c("Labels" = "Labels", "Color" = "Color", "Cluster" = "Cluster", "Body" = "Body")
                             #),
                             fluidRow(
                             column(width = 2, checkboxInput(inputId = "heatmap.visual.labels", label = "Labels", value = F)),
                             column(width = 2, checkboxInput(inputId = "heatmap.visual.color", label = "Color", value = F)),
                             column(width = 2, checkboxInput(inputId = "heatmap.visual.cluster", label = "Cluster", value = F)),
                             column(width = 2, checkboxInput(inputId = "heatmap.visual.resize", label = "Re-size", value = F)),
                             column(width = 2, checkboxInput(inputId = "heatmap.visual.body", label = "Body", value = F))
                             ),
                             
                             uiOutput("heatmap.labels.element"),
                             uiOutput("heatmap.color.element"),
                             uiOutput("heatmap.cluster.element"),
                             uiOutput("heatmap.resize.element"),
                             uiOutput("heatmap.body.element")
                             
             ),
             bsCollapsePanel("Download",
                             downloadButton("heatmap.download_png", "PNG"), 
                             downloadButton("heatmap.download_pdf", "PDF")
                             #downloadButton("heatmap.download_svg", "SVG")
             )
  )
})

observeEvent(input$heatmap.visual.labels,
             output$heatmap.labels.element<- renderUI({
               if( input$heatmap.visual.labels==T){
                 fluidRow(
                   column(width = 4,
                          sliderInput("heatmap.rows.label.size", "Rows labels size", min=10, max=24, value=12, step=2, width = "80%"),
                          sliderInput("heatmap.columns.label.size", "columns labels size", min=10, max=24, value=12, step=2, width = "80%")
                   ),
                   column(width = 4,
                          sliderInput("heatmap.rows.label.rot", "Rows labels rotation", min=0, max=90, value=0, step=5, width = "80%"),
                          sliderInput("heatmap.columns.label.rot", "columns labels rotation", min=0, max=90, value=0, step=5, width = "80%")
                          
                   )
                 )
               }
             })
)

observeEvent(input$heatmap.visual.color,
             output$heatmap.color.element<- renderUI({
               if(input$heatmap.visual.color==T){
                 
                 selectInput("heatmap.color.list", label = "Heatmap color", selected = "A", width = "50%",
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
observeEvent(input$heatmap.visual.cluster,
             output$heatmap.cluster.element<- renderUI({
               if(input$heatmap.visual.cluster == T){
                 fluidRow(
                   checkboxInput("heatmap.cluster.rows", label = "Cluster rows", value = T),
                   checkboxInput("heatmap.cluster.columns", label = "Cluster columns", value = T)
                 )
               }
             })
)

observeEvent(input$heatmap.visual.resize,
             output$heatmap.resize.element<- renderUI({
               if( input$heatmap.visual.resize==T){
                 
                 fluidRow(
                   column(width = 3,
                          sliderInput("heatmap.resize.width", "Width", min=4, max=50, value=15, step=1, width = "80%")
                   ),
                   column(width = 3,
                          sliderInput("heatmap.resize.height", "Height", min=4, max=20, value=15, step=1, width = "80%")
                   )
                 )
               }
             })
)
observeEvent(input$heatmap.visual.body,
             output$heatmap.body.element<- renderUI({
               if(input$heatmap.visual.body == T){
                 checkboxInput("heatmap.cells.value", label = "Cells values", value = T)
                 }
             })
)

observeEvent(input$metamat_input,
             updateSelectInput(session, "heatmap.factor", choices = colnames(meta_file())
             )
)

#Change factor's order
observeEvent(input$metamat_input,
             output$heatmap.order<- renderUI({
               if(input$heatmap.change.order == T)
                 selectInput("heatmap.new.order", 
                             choices = change.factor.order(physeq(), input$heatmap.factor), 
                             label = NULL, multiple = T
                 )
             })
)

#Show show table button
observeEvent(input$heatmap.go,
             output$heatmap.show.tab<- renderUI({
               actionButton("heatmap.show.tab.key", "Show table")
             })
)

#Heatmap
heatmap.plot <- function(){
  showNotification("Initialize the Heatmap.", duration = 3, type = "message")
  dataset <- data.plot(isolate(physeq()),
                       group_name = isolate(input$heatmap.factor),
                       change_fac = isolate(input$heatmap.change.order),
                       ordered_fac = isolate(input$heatmap.new.order),
                       threshold = isolate(input$heatmap.threshold),
                       level_taxa = isolate(input$heatmap.level.taxa))
  
  # Heatmap with abundance
  heatmap.data <- as.matrix(dataset$AVG_tab[,2:(dim(dataset$AVG_tab)[2]-1)])
  if(input$heatmap.rot == T)
    heatmap.data <- t(heatmap.data)
  
  if(is.null(input$heatmap.color.list)) heatmap.viridiscolor = "F"
  else heatmap.viridiscolor = input$heatmap.color.list
  
  coll_fun = colorRamp2(c(0,100), viridis(2, direction = -1,option = heatmap.viridiscolor, begin = 0.4))
  
  if(!is.null(input$heatmap.cells.value))
    if(input$heatmap.cells.value == T)
      heatmap.cell.fun = function(j, i, x, y, width, height, fill) {
        if(heatmap.data[i, j] >= 0.1)
          grid.text(sprintf("%.1f", heatmap.data[i, j]), x, y, gp = gpar(fontsize = 10))
        } 
  else heatmap.cell.fun = NULL
  
  if(is.null(input$heatmap.cluster.rows)) cluster.rows = F
  else cluster.rows = input$heatmap.cluster.rows
  
  if(is.null(input$heatmap.cluster.columns)) cluster.columns = F
  else cluster.columns = input$heatmap.cluster.columns
  
  if(is.null(input$heatmap.rows.label.size)) heatmap.rows.size = 10
  else heatmap.rows.size = input$heatmap.rows.label.size
  
  if(is.null(input$heatmap.columns.label.size)) heatmap.columns.size = 10
  else heatmap.columns.size = input$heatmap.columns.label.size
  
  if(is.null(input$heatmap.rows.label.rot)) heatmap.rows.rot = 0
  else heatmap.rows.rot = input$heatmap.rows.label.rot
  
  if(is.null(input$heatmap.columns.label.rot)) heatmap.columns.rot = 0
  else heatmap.columns.rot = input$heatmap.columns.label.rot
  
  if(is.null(input$heatmap.resize.width)) heatmap.width = 15
  else heatmap.width = input$heatmap.resize.width
  
  if(is.null(input$heatmap.resize.height)) heatmap.height = 10
  else heatmap.height = input$heatmap.resize.height
  
  showNotification("Initialize the heatmap.", duration = 2, type = "message")
  
  ht <- Heatmap(heatmap.data, name = "RB (%)", clustering_distance_rows = "pearson", clustering_distance_columns = "pearson", border = TRUE,
                cluster_rows =  cluster.rows, cluster_columns = cluster.columns,
                row_names_gp = gpar(fontsize = heatmap.rows.size), column_names_gp = gpar(fontsize = heatmap.columns.size),
                row_names_rot = heatmap.rows.rot, column_names_rot = heatmap.columns.rot,
                column_names_centered = T,
                col = coll_fun,
                heatmap_width = unit(heatmap.width, "cm"), 
                #heatmap_height = unit(heatmap.height, "cm"),
                cell_fun = switch(exists("heatmap.cell.fun"), heatmap.cell.fun, NULL)
  )  
  
 
  ht
}

#View Relative abundance Plot -Heatmap-
shiny_env = new.env()
observeEvent(input$heatmap.go, { #makeInteractiveComplexHeatmap(input, output, session, ht_list = draw(heatmap.plot()), heatmap_id = "heatmap")
               output$heatmap <- renderPlot({heatmap.plot()})
               }
)
#View relative abundance table --Heatmap-- Expanded mod
observeEvent(input$heatmap.show.tab.key, {
  showModal(modalDialog(
    title = "Relative abundance table",
    div(style="display: inline-block;vertical-align:top; width: 250px;",checkboxInput(inputId = "heatmap.add.se", label = "Add standard error", value = F)),
    div(style="display: inline-block;vertical-align:top; width: 100px;",downloadButton("heatmap.download.rb.tab", "Download")),
    div(style = 'overflow-x: scroll',DT::dataTableOutput("heatmap.avg.tab.expand",width = "100%")),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$heatmap.go, {
  output$heatmap.avg.tab.expand <- renderDataTable({
    dataset <- data.plot(isolate(physeq()),
                         group_name = isolate(input$heatmap.factor),
                         change_fac = isolate(input$heatmap.change.order),
                         ordered_fac = isolate(input$heatmap.new.order),
                         threshold = isolate(input$heatmap.threshold),
                         level_taxa = isolate(input$heatmap.level.taxa))
    AVG_tab <- dataset$AVG_tab
    AVG_tab <- round(AVG_tab[2:(dim(AVG_tab)[2]-1)], digits = 3)
    SE_tab <- dataset$SE_tab
    SE_tab <- round(SE_tab[2:dim(SE_tab)[2]], digits = 2)
    merged.tab <- merge.tab.se.fun(tab = AVG_tab, se.tab = SE_tab, add.se = input$heatmap.add.se)
    DT::datatable(merged.tab, rownames = T, escape = F, options = list(pageLength = 10))
  })
})

#####Download Heatmap
output$heatmap.download_pdf <- downloadHandler(
  filename = function() { paste('heatmap_', isolate(input$heatmap.factor), '.pdf', sep='') },
  content = function(file) {
    heatmap.width = 15
    heatmap.height = 10
    if(!is.null(input$heatmap.resize.width) & !is.null(input$heatmap.resize.height)) {
      heatmap.width = input$heatmap.resize.width
      heatmap.height = input$heatmap.resize.height
    }
    #else {} 
    library(measurements)
    heatmap.width = conv_unit(heatmap.width, "cm", "inch")
    heatmap.height = conv_unit(heatmap.height, "cm", "inch")
    
    pdf(file, onefile = T,  height = heatmap.height*1.25, width = heatmap.width*1.25)
    print(heatmap.plot())
    dev.off()
    
  }
)

output$heatmap.download_png <- downloadHandler(
  filename = function() { paste('heatmap_', isolate(input$heatmap.factor), '.png', sep='') },
  content = function(file) {
    heatmap.width = 15
    heatmap.height = 10
    if(!is.null(input$heatmap.resize.width) & !is.null(input$heatmap.resize.height)) {
      heatmap.width = input$heatmap.resize.width
      heatmap.height = input$heatmap.resize.height
    }
    library(measurements)
    heatmap.width.png = conv_unit(heatmap.width, "cm", "inch")*1.25
    heatmap.height.png = conv_unit(heatmap.height, "cm", "inch")*1.25
    
    grDevices::png(file, units = "px", height = heatmap.height.png*96, width = heatmap.width.png*96)
    print(heatmap.plot())
    dev.off()
    
  }
)

output$heatmap.download.rb.tab <- downloadHandler(
  filename = function() { paste('heatmap_Average_table', isolate(input$heatmap.factor), '.csv', sep='') },
  content = function(file) {
    dataset <- data.plot(isolate(physeq()),
                         group_name = isolate(input$heatmap.factor),
                         change_fac = isolate(input$heatmap.change.order),
                         ordered_fac = isolate(input$heatmap.new.order),
                         threshold = isolate(input$heatmap.threshold),
                         level_taxa = isolate(input$heatmap.level.taxa))
    AVG_tab <- dataset$AVG_tab
    AVG_tab <- round(AVG_tab[2:(dim(AVG_tab)[2]-1)], digits = 3)
    SE_tab <- dataset$SE_tab
    SE_tab <- round(SE_tab[2:dim(SE_tab)[2]], digits = 2)
    save.rb.table <- merge.tab.se.fun(tab = AVG_tab, se.tab = SE_tab, add.se = input$heatmap.add.se)
    
    write.csv(x = save.rb.table, file = file, quote = F, sep = "\t", row.names = T, col.names = T)
  }
)


######################################################
#################  Barchart  #################
######################################################
output$param.barchart <- renderUI({
  bsCollapse(id = "collapse.param.barchart", open = "",
             bsCollapsePanel("Visual parameters",
                             fluidRow(
                               column(width = 2, checkboxInput(inputId = "barchart.visual.labels", label = "Labels", value = F)),
                               column(width = 2, checkboxInput(inputId = "barchart.visual.color", label = "Color", value = F))
                             ),
                             
                             uiOutput("barchart.labels.element"),
                             uiOutput("barchart.color.element")
                             
             ),
             bsCollapsePanel("Download",
                             fluidRow(
                               column(width = 4,
                                      selectInput("barchart.plot.format", "Format", choices = c("PNG" = "png", "PDF" = "pdf"), width = "80%")
                               ),
                               column(width = 4,
                                      sliderInput("barchart.download.height", "Height (pixel)", min=500, max=1200, value=500,step=20, width = "80%")
                               ),
                               column(width = 4,
                                      sliderInput("barchart.download.width", "Width (pixel)", min=500, max=1200, value=500,step=20, width = "80%")
                               )
                             ),
                             downloadButton("barplot.download.plot", "")
                             #downloadButton("barplot.download_png", "PNG"), 
                             #downloadButton("barplot.download_pdf", "PDF")
             )
  )
})

observeEvent(input$barchart.visual.labels,
             output$barchart.labels.element<- renderUI({
               if( input$barchart.visual.labels==T){
                 
                 fluidRow(
                   column(width = 4,
                          sliderInput("barchart.columns.label.size", "Labels size", min=10, max=24, value=12,step=2, width = "80%")
                   ),
                   column(width = 4,
                          sliderInput("barchart.columns.label.rot", " Column lable rotation", min=0, max=90, value=45,step=10, width = "80%")
                   ),
                   column(width = 4,
                          textInput("barplot.title", "Plot title", "Enter title...")
                   )
                 )
               }
             })
)

observeEvent(input$barchart.visual.color,
             output$barchart.color.element<- renderUI({
               if(input$barchart.visual.color == T){
                 
                 selectInput("barchart.color.list", label = "Barchart color", selected = "A", width = "50%",
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

observeEvent(input$metamat_input,
             updateSelectInput(session, "barchart.factor", choices = colnames(meta_file())
             )
)

#Change factor's order
observeEvent(input$metamat_input,
             output$barchart.order<- renderUI({
               if(input$barchart.change.order == T)
                 selectInput("barchart.new.order", 
                             choices = change.factor.order(physeq(), input$barchart.factor), 
                             label = NULL, multiple = T
                 )
             })
)

#Show show table button
observeEvent(input$barchart.go,
             output$barchart.show.tab<- renderUI({
               actionButton("barchart.show.tab.key", "Show table")
             })
)

# Prepare data table for barchart graph
barchart.plot <- function(){
  showNotification("Initialize the Barchart.", duration = 2, type = "message")
  dataset <- barchart.data(isolate(physeq()),
                        group_name = isolate(input$barchart.factor),
                        change_fac = isolate(input$barchart.change.order),
                        ordered_fac = isolate(input$barchart.new.order),
                        n = isolate(input$barchart.num.OTUs),
                        level_taxa = isolate(input$barchart.level.taxa))
  average_table <- dataset$average_table
  se_table <- dataset$se_table
  
  groupsnum <- dim(average_table)[2]-1
  upper = as.matrix(average_table[,1:groupsnum, drop=F] + se_table)
  lower = as.matrix(average_table[,1:groupsnum, drop=F])
  p1 <- barchart.plot.main(average_table[,1:groupsnum, drop=F], upper, lower)  
  
  if(is.null(input$barchart.color.list)) barchart.viridiscolor = "H"
  else barchart.viridiscolor = input$barchart.color.list
  p1 + 
    ggtitle(input$barplot.title) + 
    scale_fill_viridis(discrete = TRUE, option = barchart.viridiscolor, direction = -1)+
    theme(axis.text = element_text(size = input$barchart.columns.label.size),
          axis.title = element_text(size = input$barchart.columns.label.size),
          axis.text.x = element_text(angle = input$barchart.columns.label.rot),
          legend.text = element_text(size = input$barchart.columns.label.size),
          legend.title = element_text(size = input$barchart.columns.label.size),
          text = element_text(family="serif"),
          plot.title = element_text(size = input$barchart.columns.label.size, hjust = 0.5) 
    ) 
}

#View relative abundance Plot
observeEvent(input$barchart.go, 
             {
               output$barchart <- renderPlotly({
                 ggplotly(barchart.plot())
                 })
               }
             )

#View relative abundance table --barchart-- Expanded mod
observeEvent(input$barchart.show.tab.key, {
  showModal(modalDialog(
    title = "Relative abundance table",
    #div(style="display: inline-block;vertical-align:top; width: 250px;",p("Relative abundance table", style = "font-size: 18px")), 
    div(style="display: inline-block;vertical-align:top; width: 250px;",checkboxInput(inputId = "barchart.add.se", label = "Add standard error", value = F)),
    div(style="display: inline-block;vertical-align:top; width: 100px;",downloadButton("barplot.download.rb.tab", "Download")),
    div(style = 'overflow-x: scroll',DT::dataTableOutput("barchart.avg.tab.expand",width = "100%")),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$barchart.go, {
  output$barchart.avg.tab.expand <- renderDataTable({
    dataset <- barchart.data(isolate(physeq()),
                          group_name = isolate(input$barchart.factor),
                          change_fac = isolate(input$barchart.change.order),
                          ordered_fac = isolate(input$barchart.new.order),
                          n = isolate(input$barchart.num.OTUs),
                          level_taxa = isolate(input$barchart.level.taxa))
    average_table <- dataset$average_table
    average_table <- round(average_table[1:dim(average_table)[2]-1], digits = 3)
    se_table <- round(dataset$se_table, digits = 3)
    merged.tab <- merge.tab.se.fun(tab = average_table, se.tab = se_table, add.se = input$barchart.add.se)
    DT::datatable(merged.tab, rownames = T, escape = F, options = list(pageLength = 10))
  })
})


#Downloads
#Download plot
output$barplot.download.plot<- downloadHandler(
  filename = function() { paste('barchart', isolate(input$barchart.factor), '.', isolate(input$barchart.plot.format) , sep='') },
  content = function(file) {
    switch(input$barchart.plot.format,
           "png" = png(file, height = input$barchart.download.height, width = input$barchart.download.width),
           "pdf" = pdf(file, onefile=T, height = input$barchart.download.height, width = input$barchart.download.width),
           "svg" = svg(file, height = input$barchart.download.height, width = input$barchart.download.width)
    )
    print(barchart.plot())
    dev.off()
  }
)

#Download relative abundance table
output$barplot.download.rb.tab <- downloadHandler(
  filename = function() { paste('barchart_rb_table', isolate(input$barchart.factor), '.csv', sep='') },
  content = function(file) {
    dataset <- barchart.data(isolate(physeq()),
                          group_name = isolate(input$barchart.factor),
                          change_fac = isolate(input$barchart.change.order),
                          ordered_fac = isolate(input$barchart.new.order),
                          n = isolate(input$barchart.num.OTUs),
                          level_taxa = isolate(input$barchart.level.taxa))
    average_table <- dataset$average_table
    average_table <- round(average_table[1:dim(average_table)[2]-1], digits = 3)
    se_table <- round(dataset$se_table, digits = 3)
    save.rb.table <- merge.tab.se.fun(tab = average_table, se.tab = se_table, add.se = input$barchart.add.se)
    
    write.csv(x = save.rb.table, file = file, quote = F, sep = "\t", row.names = T, col.names = T)
  }
)