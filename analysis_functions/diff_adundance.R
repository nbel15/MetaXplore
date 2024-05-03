# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
library(microeco)

# Calculate MPSE abundance and Perform differential analysis 
  diff.abund.analysis <- function(physeq, group.fac, level.taxa, abund.threshold){
    # Filter taxa with low abundance 
    physeq.filt <- otus.abund.filt(physeq = physeq, abund.threshold = abund.threshold)
    
    # convert phyloseq object to mpse
    data.micro <- file2meco::phyloseq2meco(physeq = physeq.filt)
    
    #calculate abundance
    t2 <- trans_diff$new(dataset = data.micro, method = "wilcox", group = group.fac, taxa_level = level.taxa, filter_thres = 0.001)
    
    #simplify name
    t2$res_diff$Taxa %<>% gsub(".*\\|", "", .)
    t2$res_abund$Taxa %<>% gsub(".*\\|", "", .)
    
    #remove prefix
    t2$res_diff$Taxa %<>% gsub(".__", "", .)
    t2$res_abund$Taxa %<>% gsub(".__", "", .)
    
    # filter something not needed to show
    t2$res_diff %<>% subset(Significance %in% c("***", "**", "*"))
    
    return(t2)
  }

# Plot results 
  diff.abund.plot<- function(t2, select_taxa, angle = 90, ajust = 0.5, label.size =14){
    
    #select_taxa = "g:Fusibacter"
   # diff.abund.select.taxa= "g:Fusibacter"
    #diff.abund.factor="Age2"
    #diff.abund.level.taxa = "Genus"
    #t2 <- diff.abund.analysis(physeq = physeq, group.fac = diff.abund.factor, level.taxa = diff.abund.level.taxa) 
    
    #extract abundance and differential results
    abund_data <- t2$res_abund
    diff_data <- t2$res_diff
  
    # Reorder taxa according to P.adj in res_diff from low to high
    diff_data %<>% .[order(.$P.adj, decreasing = FALSE), ]
   
    #select_taxa_list <- unique(diff_data$Taxa)
    if (select_taxa == "all") {
      all_taxa <- 1:length(unique(as.character(diff_data$Taxa)))
      diff_data %<>% .[.$Taxa %in% unique(as.character(diff_data$Taxa))[all_taxa], ]
      diff_data$Taxa %<>% factor(., levels = rev(unique(as.character(.))))
    } else {
      message("Use provided select_taxa to filter and reorder taxa ...")
      diff_data %<>% .[.$Taxa %in% select_taxa, ]
      diff_data$Taxa %<>% factor(., levels = rev(select_taxa))
    }  
    
    abund_data %<>% .[.$Taxa %in% levels(diff_data$Taxa), ]
    abund_data$Taxa %<>% factor(., levels = sort(levels(diff_data$Taxa)))
    
    diff_data[, "Significance"] %<>% as.character
    
    all_taxa <- levels(abund_data$Taxa)
    
    x_axis_order <- unique(abund_data$Group)
    annotations <- c()
    x_min <- c()
    x_max <- c()
    y_position <- c()
    start_bar_mid <- 1 - (0.8/2 - 0.8/(length(x_axis_order) * 2))
    increase_bar_mid <- 0.8/length(x_axis_order)
    
    for (j in all_taxa) {
      select_use_diff_data <- diff_data %>% dropallfactors %>% 
        .[.$Taxa == j, ]
      select_use_abund_data <- abund_data %>% dropallfactors %>% 
        .[.$Taxa == j, ]
      y_start_use <- max((select_use_abund_data$Mean + select_use_abund_data$SE) )# * 1.5
      for (i in seq_len(nrow(select_use_diff_data))) {
        mid_num <- match(j, all_taxa) - 1
        annotations %<>% c(., select_use_diff_data[i, 
                                                   "Significance"])
        x_min %<>% c(., mid_num + (start_bar_mid + 
                                     (match(gsub("(.*)\\s-\\s(.*)", "\\1", 
                                                 select_use_diff_data[i, "Comparison"]), 
                                            x_axis_order) - 1) * increase_bar_mid))
        x_max %<>% c(., mid_num + (start_bar_mid + 
                                     (match(gsub("(.*)\\s-\\s(.*)", "\\2", 
                                                 select_use_diff_data[i, "Comparison"]), 
                                            x_axis_order) - 1) * increase_bar_mid))
        y_position %<>% c(., y_start_use + 0.01 * i)
        
      }
    }
    
    color_values = viridis::viridis(n = length(unique(abund_data$Group)) , option = "D", direction = -1 )
    
    p <- ggplot(abund_data, aes(x = Taxa, y = Mean, color = Group, fill = Group)) + 
      theme_bw() + 
      geom_bar(stat = "identity", position = position_dodge(), width = 0.8, colour="black") +  
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), 
                    width = 0.45, 
                    position = position_dodge(0.8), color = "black") +
      scale_color_manual(values = color_values) + scale_fill_manual(values = color_values) + 
      ylab("Relative abundance") + theme(legend.position = "right") +
      theme(axis.text.x = element_text(angle=angle, vjust=ajust, hjust=ajust),
            axis.text = element_text(size = label.size), 
            axis.title = element_text(size = label.size), 
            legend.text = element_text(size = label.size))+
      scale_y_continuous(expand =  expansion(mult = c(0, .03))) 
    
    p <- p + ggsignif::geom_signif(annotations = annotations, 
                                   y_position = y_position, xmin = x_min, xmax = x_max, 
                                   color = "black", tip_length = 0.01, vjust = 0.9)
    
    return(p)
  }

# Export differential abundance table
diff.abund.table.export <- function(physeq, group.fac, level.taxa, select_taxa, abund.threshold){
    t2 <- diff.abund.analysis(physeq = physeq,
                              group.fac = group.fac, 
                              level.taxa = level.taxa,
                              abund.threshold = abund.threshold)  
    t2$res_diff$P.adj <- signif(t2$res_diff$P.adj, digits=2)
    out.tab <- t2$res_diff %>%
      arrange(P.adj) %>% 
      select(c("Comparison", "Taxa", "P.adj", "Significance"))
      
      if(select_taxa != "all")
        out.tab <- out.tab %>%
                      filter(Taxa == select_taxa)
    
    return(out.tab)
}
