# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########
# Guide #
#########
output$metaxplore_workflow <- renderUI({
  bsCollapse(id = "collapseqplot", open = "",
             bsCollapsePanel("Overview", htmlOutput("overview_pipeline"), style = "primary"),
             bsCollapsePanel("Import Data", htmlOutput("import_data"), style = "primary"),
             bsCollapsePanel("Alpha Diversity", htmlOutput("alphaguide"),style = "warning"),
             bsCollapsePanel("Beta Diversity", htmlOutput("betaguide"),style = "success"),
             bsCollapsePanel("Relative Abundance", htmlOutput("rb_guide"), style = "danger"),
             bsCollapsePanel("Differential abundance", htmlOutput("diffabun_guide"), style = "danger"),
             bsCollapsePanel("Core microbiome", htmlOutput("coremicro_guide"), style = "info"),
             bsCollapsePanel("Biomarker discovery", htmlOutput("biomarker_guide"), style = "success")
             )
})

output$overview_pipeline <- renderUI(includeHTML("www/overview.html"))
output$import_data <- renderUI(includeHTML("www/import_data.html"))
output$alphaguide <- renderUI(includeHTML("www/alpha_div.html"))
output$betaguide <- renderUI(includeHTML("www/beta_div.html"))
output$rb_guide <- renderUI(includeHTML("www/r_abundance.html"))
output$diffabun_guide <- renderUI(includeHTML("www/diff_abundance.html"))
output$coremicro_guide <- renderUI(includeHTML("www/core_micro.html"))
output$biomarker_guide <- renderUI(includeHTML("www/biomarker.html"))

