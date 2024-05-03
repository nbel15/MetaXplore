#
source("auto_install_load.R")

library(jtools)
tags$head(tags$link(rel="shortcut icon", href="favicon.ico"))
dashboardPage(
  dashboardHeader(
    title = "MetaXplore"
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Guide", tabName = "guide", icon = icon("book-reader")),
      menuItem("Import Data", tabName = "import_data", icon = icon("file-import")),
      menuItem("Alpha diversity", tabName = "alpha_tab", icon = icon("ruler-combined")),
      menuItem("Beta diversity", tabName = "beta_div_graph", icon = icon("project-diagram")),
      menuItem("Relative Abundance", tabName = "rb_tab", icon = icon("poll")),
      menuItem("Differential abundance", tabName = "diff_abund_tab", icon = icon("object-ungroup")),
      menuItem("Core microbiome", tabName = "core_tab", icon = icon("linode")),
      menuItem("Biomarker discovery", tabName = "biomarker_tab", icon = icon("fire")),
      img(src = "mgemlab3-outline.png", height = 150, width = 220,align = "center", style="position:absolute;bottom:0;")
      
    )
  ), 
  dashboardBody(
    useShinyjs(),
    useShinyalert(),
    
    tabItems(
      guide,
      import,
      alphatab,
      betatab,
      rbtab,
      diffabundtab,
      coretab,
      biomarkertab
      )
    )
  ) 

