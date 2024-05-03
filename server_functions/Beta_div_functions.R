# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# import requested functions

source("tools_functions.R")
source("analysis_functions/beta_div.R")

#########################################
# Activate Go button conditions

observe({
  if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input)) {# || is.null(input$tree_file_input)) {
    shinyjs::disable("beta.go")
  } else {
    shinyjs::enable("beta.go")
    }
  })

#observeEvent(input$metamat_input, {
  
 # metadata <- sample_data(isolate(physeq()))
  
  #if(unique(metadata[,beta.div.factor]) > 2){
   # enable("beta.go")}
#})

#########################################
#Update UI
# Load list of factors for permanova
observeEvent(input$metamat_input,
             updateSelectInput(session, "beta.div.factor", choices = colnames(sample_data(isolate(physeq()))) #colnames(meta_file())
             )
)

observeEvent(input$collapse.param.beta.ordi,
             updateSelectInput(session, "Ordi.select.samples", choices = rownames(sample_data(isolate(physeq())))
             )
)

# Select m for CAP
observeEvent(input$metamat_input,
             output$beta.div.cap.param<- renderUI({
               if(input$ordination.method == "cap"){
                 fluidRow(
                 numericInput("beta.div.cap.m", "m (lda)", value = 11, min = 0, max = 30, step = 1, width = "30%")#,
                 #actionButton("beta.show.best.m", "Best m value")
                 )
               }
             })
)

# Select shape factor 
observeEvent(input$beta.add.shape,
             output$beta.shape.factor<- renderUI({
               if(input$beta.add.shape == T)
                 selectInput("beta.factor.shape", "Shape mapping factor", choices = "choose")
             })
)
# Load list for shape factor
shape.factor <- function(){
  list.factot2 <- colnames(meta_file())
  list <- list.factot2[!list.factot2 %in% input$beta.factor.color]
  return(list)
}
observeEvent(input$beta.factor.color,
             updateSelectInput(session, "beta.factor.shape", choices = colnames(meta_file())
             )
)
#########################################
# Show distance matrix
#Show show graph button
observeEvent(input$beta.go,
             output$beta.show.dist.tab<- renderUI({
               actionButton("beta.show.dist.key", "Dissimilarity Matrix")
             })
)

#########################################
# Show best m in case of CAP analysis
#Show show table button
observeEvent(input$beta.show.best.m, {
  showModal(modalDialog(
    title = paste("Number of PCO axes (m) : ", sep = ""),
    #div(style="display: inline-block;vertical-align:top; width: 100px;",downloadButton("beta.dist.matrix", "Download")),
    #div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput(inputId = "beta.maxm", label = "max (m)", value = 5, min = 2, max = 50,step = 1)),
    plotOutput("beta.m.graph"),
    easyClose = TRUE,
    size = "l"
  ))
})

##Number of PCO axes (m) graph
observeEvent(input$beta.show.best.m,
             output$beta.m.graph<- renderPlot({
               best.m.graph(physeq = isolate(physeq()), 
                            div.factor = isolate(input$beta.div.factor),
                            permut = isolate(input$permut.permanova),
                            dist.method = isolate(input$dist.method)
                              )
               
             })
)



#########################################
# View Beta Results
# Dissimilarity matrix Expanded mod
observeEvent(input$beta.show.dist.key, {
  showModal(modalDialog(
    title = "Dissimilarity matrix: ",
    div(style="display: inline-block;vertical-align:top; width: 100px;",downloadButton("beta.dist.matrix", "Download")),
    div(style = 'overflow-x: scroll',DT::dataTableOutput("beta.dist.matrix.expand",width = "100%")),
    easyClose = TRUE,
    size = "l"
  ))
})

##Distance matrix
observeEvent(input$beta.go,
             output$beta.dist.matrix.expand<- renderDataTable({
               beta.dist(isolate(physeq()),
                         dist.method = isolate(input$dist.method)
                         )
               })
             )

##Download Distance matrix
##Distance matrix
output$beta.dist.matrix <- downloadHandler(
  filename = function() { paste(input$weight,'Dist_table', '.csv', sep='')},
  content = function(file) {
    dist_mat <- beta.dist(isolate(physeq()),
                          dist.method = isolate(input$dist.method)
                          )
    write.csv(x = dist_mat, file = file, quote = F, sep = "\t", row.names = T, col.names = NA)
  })


# Parameter section beta diversity plot
output$beta.div.plot.param <- renderUI({
  bsCollapse(id = "collapse.param.beta.ordi", open = "",
             bsCollapsePanel("Visual parameters",
                             fluidRow(
                               selectInput("Ordi.select.samples", label = "Select samples", choices = c("choose.."), multiple = T)),
                             fluidRow(
                               column(width = 3, checkboxInput(inputId = "ordi.ellipse", label = "Add ellipse", value = F)),
                               column(width = 3, checkboxInput(inputId = "ordi.star", label = "Add star", value = F))
                             )
             ),
             bsCollapsePanel("Download",
                             downloadButton("ordi.download.png", "PNG"), 
                             downloadButton("ordi.download.pdf", "PDF")
             )
  )
})

# Parameter section Pairwise comparison plot
output$beta.div.stat1.param <- renderUI({
  bsCollapse(id = "collapse.param.beta.stat1", open = "",
             bsCollapsePanel("Visual parameters",
                             fluidRow(
                               column(width = 4,
                                      sliderInput("p.val.plot.label.size", "Labels size", min=10, max=24, value=14,step=2, width = "80%")),
                               column(width = 4,
                                      sliderInput("p.val.body.label.size", "Labels size", min=2, max=10, value=3,step=1, width = "80%")),
                               column(width = 4,
                                      sliderInput("p.val.plot.xlabel.rotation", "Lable rotation (x.axis)", min=0, max=90, value=0,step=10, width = "100%"))
                             )
             ),
             bsCollapsePanel("Download",
                             downloadButton("pair.stat1.download.png", "PNG"), 
                             downloadButton("pair.stat1.download.pdf", "PDF")
             )
  )
})




beta.div.plot.fun <- function(){
  showNotification("Initialize plot.", duration = 2, type = "message")
  beta.ordi.list <- beta.ordi(physeq = isolate(physeq()), 
                              dist.method = isolate(input$dist.method), 
                              div.factor = isolate(input$beta.div.factor), 
                              ordination.method = isolate(input$ordination.method), 
                              m = isolate(input$beta.div.cap.m),
                              permut = isolate(input$permut.permanova)
                              )
  if(!is.null(beta.ordi.list)){
  
  P <- ggscatter(beta.ordi.list[["ordi.tab"]], x = "axis1", y = "axis2",# facet.by = "Sex",
            color =  paste(isolate(input$beta.div.factor)), palette = get_palette(palette = "jco", 14), #"jco", 
            size = 5, ellipse.alpha = 0.2, ggtheme = theme_bw(),
            ellipse = ifelse(!is.null(input$ordi.ellipse), yes = input$ordi.ellipse, no = F), 
            ellipse.type = "confidence", #ellipse.level = 0.5,
            mean.point = F, 
            star.plot = ifelse(!is.null(input$ordi.star), yes = input$ordi.star, no = F),
            xlab = beta.ordi.list[["xlabel"]], 
            ylab = beta.ordi.list[["ylabel"]], 
            label = beta.ordi.list[["ordi.tab"]]$labels, 
            label.select = ifelse(!is.na(input$Ordi.select.samples), yes = input$Ordi.select.samples, no = c("")),
            font.label = c(14, "bold"), repel = T,
            font.tickslab="bold", #facet.by = "Sex",
            #title = "title", 
            caption = paste("Plot of Microbial Profiles","\n","(p-value: ", beta.ordi.list$adonis,")",sep="")
  )
  return(P)
  } else {
    shinyalert("Oops!", "Select factor with more than 2 elements", type = "warning")
    return(NULL)
  }
}


###Beta Diversity plots
observeEvent(input$beta.go,
             output$beta.div.plot<- renderPlot({
               beta.div.plot.fun()
               })
)

#########################################
#Download Plots

###Download plot PNG
output$ordi.download.png <- downloadHandler(
  filename = function() { paste(input$ordination.method, '-', input$beta.div.factor,'-plot', '.png', sep='') },
  content = function(file) {
    png(file, height = 500, width = 700, res = 100)
    print(beta.div.plot.fun())
    dev.off()
  })

###Download plot PDF
output$ordi.download.pdf <- downloadHandler(
  filename = function() { paste(input$ordination.method, '-', input$beta.div.factor,'-plot', '.pdf', sep='') },
  content = function(file) {
    pdf(file, onefile=F,paper = "a4r", height = 500, width = 700)
    print(beta.div.plot.fun())
    dev.off()
  })



beta.div.pair.plot.fun <- function(){
  showNotification("Initialize plot.", duration = 2, type = "message")
  beta.pair.plot <- beta.pair.basic(physeq = isolate(physeq()),
                                    dist.method = isolate(input$dist.method),
                                    div.factor = isolate(input$beta.div.factor),
                                    permut = isolate(input$permut.permanova)
                                    )
  if(!is.null(beta.pair.plot)){
    
    P <- ggplot(data = beta.pair.plot[["df.select"]], aes(Groups1, Groups2, fill = p_range))+
      geom_tile(linetype = 1, color = "white", lwd = 1)+
      scale_fill_manual(name = "p-value \n range", values = beta.pair.plot[["col_breaks"]])+
                        #labels = c("p\u22640.001", "0.001<p\u22640.01", "0.01<p\u22640.05", "p>0.05"))+
      theme_minimal()+ 
      geom_text(aes(label = Pval), color = "black", size = input$p.val.body.label.size)+
      coord_fixed()+
      labs(title = "PERMANOVA: Multiple comparison")+
      theme(axis.text.x = element_text(vjust = 1, hjust = 0.5, face="bold",
                                       size = input$p.val.plot.label.size, 
                                       angle = input$p.val.plot.xlabel.rotation),
            axis.text.y = element_text(vjust = 1, hjust = 1, face="bold", 
                                       size = input$p.val.plot.label.size),
            axis.title = element_blank(),
            panel.grid.major = element_blank(), 
            title = element_text(size = 14, face="bold"),
            legend.text = element_text(size = 14), 
            legend.position = "bottom", 
            legend.title = element_blank())
    return(P)
  } else {
    #shinyalert("Oops!", "Select factor with more than 2 elements", type = "warning")
    return(NULL)
  }
}


###Beta Diversity pairwise comparison plots
observeEvent(input$beta.go,
             output$beta.div.pair.plot<- renderPlot({
               beta.div.pair.plot.fun()
             })
)

###Download plot pairwise comparison PNG
output$pair.stat1.download.png <- downloadHandler(
  filename = function() { paste( input$beta.div.factor,'-Pairwise-stat','-plot', '.png', sep='') },
  content = function(file) {
    png(file, height = 700, width = 900, res = 100)
    print(beta.div.pair.plot.fun())
    dev.off()
  })

###Download plot pairwise comparison PDF
output$pair.stat1.download.pdf <- downloadHandler(
  filename = function() { paste(input$beta.div.factor,'-Pairwise-stat','-plot', '.pdf', sep='') },
  content = function(file) {
    pdf(file, onefile=F,paper = "a4r", height = 500, width = 700)
    print(beta.div.pair.plot.fun())
    dev.off()
  })































###phylogram
observeEvent(input$goBeta,
             output$beta_phylo_view<- renderPlot({
               beta.phylo.tree(isolate(physeq()),
                                        group_name = isolate(input$beta_fac),
                                        weight = isolate(input$weight)
                               )
               })
             )



###Pairwise Beta diversity
observeEvent(input$goBeta,
             output$select_comb <- renderUI({
               selectInput(inputId = "comb", label = "",
                           choices = comb.select(physeq(), input$beta_fac), 
                           width = "40%"
                           )
               })
             )

observeEvent(input$goBeta,
             output$beta_pairwise_view<- renderPlot({
               df <- colnames(meta_file()[,sapply(meta_file(), function(x){length(unique(x))})>2, drop = F])
               testpair <- input$beta_fac %in% df
               
               if(!is.null(input$comb) && testpair){
                 showNotification(paste("Please wait ...."), id = "beta2", 
                                  duration = NULL, closeButton = F)
                 beta.pair.plot(isolate(physeq()),
                                group_name = isolate(input$beta_fac),
                                weight = isolate(input$weight),
                                comb = isolate(input$comb),
                                view_meth = isolate(input$view_meth)
                                )
                   removeNotification("beta2")
                  }
               })
             )

###Pairwise comparison p values
observeEvent(input$goBeta,
             output$pair_comp_table<- renderDataTable({
               out <- beta.pval.parwise(isolate(physeq()), 
                                        isolate(input$beta_fac), 
                                        isolate(input$weight))
               data.frame("Group1" = out$pair_1_list, 
                          "Group2" = out$pair_2_list,
                          "Pval" = out$pVal,
                          "PvalBH" = out$pVal_BH)
               
             })
)


#########################################
#Download Plots NMDS - MDS - PCoA - CAP

###Download area
observeEvent(input$goBeta,
             output$download_beta_plot <- renderUI({
               fluidRow(        
                 "Download current Plots:", 
                 br(),
                 downloadButton("down_beta_plots", "PDF"),
                 br(),
                 br(),
                 "Download All factors:", 
                 br(),
                 downloadButton("down_beta_all", "All")
                 )
               })
             )
###Download plot PDF
output$down_beta_plots <- downloadHandler(
  filename = function() { paste(input$beta_fac,'-beta-diversity','.pdf', sep='')},
  content = function(file) {
    pdf(file)
    print(beta.plot(isolate(physeq()),
                    group_name = isolate(input$beta_fac),
                    weight = isolate(input$weight),
                    view_meth = isolate(input$view_meth)
                    )
          )
    dev.off()
  })

###Download plot ALL PDF
output$down_beta_all <- downloadHandler(
  filename = function() { paste('Beta_div_plots','.pdf', sep='') },
  content = function(file) {
    pdf(file)
    withProgress(message = 'Download plots', value = 0, {
    for (i in 1:length(colnames(meta_file())))
      {
      group_name <- colnames(meta_file())[i]
      incProgress(1/length(colnames(meta_file())), detail = paste(": ", group_name))
      beta.plot(isolate(physeq()),
                group_name = group_name,
                weight = isolate(input$weight),
                view_meth = isolate(input$view_meth)
                )
      Sys.sleep(0.1)
    }
    })
    dev.off()
    })


#########################################
#Download Distance Matrix

###Download area
observeEvent(input$goBeta,
             output$down_mat_button<- renderUI({
               tagList(
                 hr(),
                 downloadButton("down_beta_mat", "Download table")
                 )
               })
             )
###Download Beta-div matrix
output$down_beta_mat <- downloadHandler(
  filename = function() { paste(input$weight,'Dist_table', '.csv', sep='')},
  content = function(file) {
    dist_mat <- beta.dist(isolate(physeq()),
                          weight = isolate(input$weight)
                          )
    write.csv(x = dist_mat, file = file, quote = F, sep = "\t", row.names = T, col.names = NA)
    })

#########################################
#Download pairwise Plots NMDS/MDS

###Download area
observeEvent(input$goBeta,
             output$download_beta_pairwise <- renderUI({
               fluidRow(
                 "Download Pairwise comparison plots:",
                 br(),
                 downloadButton("down_pair_plots", "PDF"),
                 br()
               )
               })
             )

###Download plot PDF
output$down_pair_plots <- downloadHandler(
  filename = function() { paste(input$beta_fac,'-Pairwise-beta-div','.pdf', sep='')},
  content = function(file) {
    pdf(file)
    df <- colnames(meta_file()[,sapply(meta_file(), function(x){length(unique(x))})>2])
    
    if(input$beta_fac %in% df){
      list_comb <- comb.select(physeq(), input$beta_fac)
      withProgress(message = 'Download plot', value = 0, {
            for(i in 1:length(list_comb)){
              incProgress(1/length(list_comb), detail = paste(": ", list_comb[i]))
              current_comb <- list_comb[i]
              beta.pair.plot(isolate(physeq()),
                             group_name = isolate(input$beta_fac),
                             weight = isolate(input$weight),
                             comb = current_comb,
                             view_meth = isolate(input$view_meth)
              )
                Sys.sleep(0.1)
              }
        })
    }
    dev.off()
  })

