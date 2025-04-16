# MetaXplore Version 1.0
# Last modified on 14/042025
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# import requested functions
source("tools_functions.R")

#########################################
# Import Otu table
#########################################
otu_file <- reactive({
  if (is.null(input$otumat_input)) {
    return(NULL)
  } else {
    ## check whether a .txt file is uploaded
    validate(
      need(tools::file_ext(input$otumat_input$name) %in% c('txt','tsv'), 
           "Wrong File Format Uploaded")
    )	
    ## Import file
    otumat <- read.table(input$otumat_input$datapath, header=TRUE, comment.char = "", stringsAsFactors = F, sep = "\t", row.names = 1, check.names = FALSE)
    # Clean table from empty lines
    otumat <- otumat[!apply(is.na(otumat) | otumat =="",1,all),]
    # Clean table from empty columns
    otumat <- otumat %>% select_if(colSums(.) != 0)
    ##order row && column names
    otumat <- otumat[order(rownames(otumat)), order(colnames(otumat))]
    return(otumat)
  }
})

#########################################
# Import Taxonomy table
#########################################
taxa_file <- reactive({
  if (is.null(input$taxamat_input)) {
    return(NULL)
  } else {
    ## check whether a .CSV or .txt file is uploaded
    validate(
      need(tools::file_ext(input$taxamat_input$name) %in% c('tsv', 'txt'), 
           "Wrong File Format Uploaded")
    )	
    ## Import file
    taxamat <- read.table(input$taxamat_input$datapath, header=TRUE, comment.char = "", stringsAsFactors = F, sep = "\t", check.names = FALSE)
    # Clean table from empty lines
    taxamat <- taxamat[!apply(is.na(taxamat) | taxamat =="",1,all),]
    taxamat$Taxon <- str_replace_all(taxamat$Taxon,'d__','d:')
    taxamat$Taxon <- str_replace_all(taxamat$Taxon,' p__','p:')
    taxamat$Taxon <- str_replace_all(taxamat$Taxon,' c__','c:')
    taxamat$Taxon <- str_replace_all(taxamat$Taxon,' o__','o:')
    taxamat$Taxon <- str_replace_all(taxamat$Taxon,' f__','f:')
    taxamat$Taxon <- str_replace_all(taxamat$Taxon,' g__','g:')
    taxamat$Taxon <- str_replace_all(taxamat$Taxon,'_','-')
    ##order row names
    taxamat <- taxamat[order(rownames(taxamat)), ]
    taxatab <- split.taxon(taxamat)
    taxatab <- taxatab[order(rownames(taxatab)),]
    taxatab
  }
})

#########################################
# Import Meta table
#########################################
meta_file <- reactive({
  if (is.null(input$metamat_input)) {
    return(NULL)
  } else {
    ## check whether a .CSV or .txt file is uploaded
    validate(
      need(tools::file_ext(input$metamat_input$name) %in% c('txt'), 
           "Wrong File Format Uploaded")
    )	
    ## Import file
    metamat <- read.table(input$metamat_input$datapath, header=TRUE, comment.char = "", stringsAsFactors = F, sep = "\t", row.names = 1, check.names = FALSE)
    # Clean table from empty lines
    metamat <- metamat[!apply(is.na(metamat) | metamat =="",1,all),]
    #order row names
    metamat <- metamat[order(rownames(metamat)),]
    #metamat_filtered <- metamat #%>%  
    #select(where(~ n_distinct(.) > 1)) %>%
    # select(where(~ n_distinct(.) < length(rownames(metamat))))
    # replace special characters with underscores to entire data.frame
    metamat_currated <- as.data.frame(lapply(metamat, replace_special_chars), stringsAsFactors = FALSE)
    rownames(metamat_currated) <- rownames(metamat)
    metamat <- metamat_currated
    metamat
  }
})

#########################################
# Import Tree table
#########################################
tree_file <- reactive({
  if (is.null(input$tree_file_input)) {
    return(NULL)
  } else {
    ## check whether a .CSV or .txt file is uploaded
    validate(
      need(tools::file_ext(input$tree_file_input$name) %in% c('nwk','NWK', 'tree', 'tre'), 
           "Wrong File Format Uploaded")
    )	
    ## Import file
    tree <- ape::read.tree(input$tree_file_input$datapath)
    #tree <- phyloseq::read_tree_greengenes(fixUploadedFilesNames(input$tree_file_input)$datapath)
    tree
  }
})

#########################################
# Build phyloseq data
#########################################
physeq <- reactive({
  #if (is.null(input$otumat_input) && is.null(input$taxamat_input) && is.null(input$metamat_input)) {
  if (!is.null(input$read_input)) {
    OTU = phyloseq::otu_table(as.matrix(otu_file()), taxa_are_rows = T)
    TAX = phyloseq::tax_table(as.matrix(taxa_file()))
    SAM = phyloseq::sample_data(as.data.frame(meta_file()))
    
    if(!is.null(tree_file())){
      TRE = phyloseq::phy_tree(tree_file())
      physeq = phyloseq(OTU, TAX, SAM, TRE)
    } else physeq = phyloseq(OTU, TAX, SAM)
    
    physeq <- prune_samples(sample_sums(physeq)>0, physeq)
    physeq <- otus.freq.filter(physeq, min.freq = isolate(input$min.freq))
    if(!is.null(physeq)){
      physeq <- otus.norm(physeq, method = isolate(input$normalize))
      physeq
    } else{
      shinyalert("Oops!", paste0("No data available! Consider reducing the frequency cutoff (< ", isolate(input$min.freq), "%)", sep = ""), 
                 type = "warning")
      return(NULL)
    } 
  } else {
    return(NULL)
  }
})


#########################################
# View Input Table
#########################################
observe({
  if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input)) {# || is.null(input$tree_file_input)) {
    shinyjs::disable("read_input")
  } else {
    shinyjs::enable("read_input")
  }
})

observeEvent(input$read_input,{
  if(!is.null(physeq())){
    output$otu_tab_view <- renderDataTable({phyloseq::otu_table(physeq())})
    output$taxa_tab_view <- renderDataTable({phyloseq::tax_table(physeq())})
    output$meta_tab_view <- renderDataTable({phyloseq::sample_data(physeq())})
  }else{
    output$otu_tab_view <- renderDataTable({datatable(data.frame("Empty data" = "No data available! Consider reducing the frequency cutoff."))})
    output$taxa_tab_view <- renderDataTable({datatable(data.frame("Empty data" = "No data available! Consider reducing the frequency cutoff."))})
    output$meta_tab_view <- renderDataTable({datatable(data.frame("Empty data" = "No data available! Consider reducing the frequency cutoff."))})
  }
  
  
  output$nsamples <- renderInfoBox({
    raw_samples <- length(rownames(meta_file()))
    if(!is.null(physeq())){
      n_samples <- length(rownames(phyloseq::sample_data(physeq())))
    }else n_samples = 0
    
    infoBox(
      title = "N. Samples",
      value = paste0(n_samples, " / ", raw_samples),
      icon = icon(name = "list"),
      color = "blue", fill = TRUE
    )
  })
  
  output$ntaxa <- renderInfoBox({
    raw_taxa <- length(rownames(taxa_file())) 
    if(!is.null(physeq())){
      taxa_file <- phyloseq::tax_table(physeq())
      taxa_file <- data.frame(taxa_file@.Data)
      n_taxa <- length(rownames(taxa_file))
    }else n_taxa = 0
    
    infoBox(
      title = "N. Taxa (Clusters/OTUs)", 
      value = paste0(n_taxa, " / ", raw_taxa), 
      icon = icon("bacterium"),
      color = "purple", fill = TRUE
    )
  })
  
  output$info_input <- renderInfoBox({
    factors <- length(colnames(meta_file()))
    infoBox(
      title = "N. Factors",
      value = paste0(factors),  
      icon = icon("info"),
      color = "orange", fill = TRUE
    )
  })
})


