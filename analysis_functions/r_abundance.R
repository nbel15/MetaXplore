# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#############################
# Import requested functions
#############################
source(file = "tools_functions.R")

# Filter OTUs below threshold
  rb.otus.all <- function(physeq){
  if(!taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  ##Store input data
  otu_file <- phyloseq::otu_table(physeq)
  otu_file <- data.frame(otu_file@.Data)
  otu_file <- data.frame(otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),])
  otu_file <- otu_file[sort(rownames(otu_file)), ]
  
  rel.otu.tab.all <- decostand(otu_table, method = "total", MARGIN = 2)*100
  rel.otu.tab.all <- rel.otu.tab.all[order(rownames(rel.otu.tab.all)),]
  
  }
  
# Get Abundance data from phyloseq object
  abundance_table <- function(physeq, level_taxa="Phylum"){
    if(!taxa_are_rows(physeq)){
      physeq <- t(physeq)
    }
    ##Store input data
    otu_file <- phyloseq::otu_table(physeq)
    otu_file <- data.frame(otu_file@.Data)
    otu_file <- data.frame(otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),])
    otu_file <- otu_file[sort(rownames(otu_file)), ]
    
    meta_file <- phyloseq::sample_data(physeq)
    mappingVar <- names(meta_file)
    meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),],row.names=row.names(meta_file))
    meta_file <- meta_file[sort(rownames(meta_file)), ]
    
    taxa_file <- phyloseq::tax_table(physeq)
    taxa_file <- data.frame(taxa_file@.Data)
    taxa_file <- taxa_file[sort(rownames(taxa_file)), ]
    
    if(length(level_taxa)==1)
      temp_otu_file <- data.frame("OTUs"=taxa_file[, level_taxa], otu_file)
    else 
      temp_otu_file <- data.frame("OTUs"=col_concat(taxa_file[, level_taxa], sep = "|"), otu_file)
    
    temp_otu_file <- aggregate(temp_otu_file[,2:length(temp_otu_file)],temp_otu_file["OTUs"],sum)
    otu_file <- temp_otu_file[,2:length(temp_otu_file)]
    rownames(otu_file) <- temp_otu_file[,"OTUs"]
    ##Calculate Relative abundance
    otu_table<- apply(X = otu_file,2, FUN = as.integer)
    rownames(otu_table) <- rownames(otu_file)
    
    rel_otu_table <- decostand(otu_table, method = "total", MARGIN = 2)*100
    rel_otu_table <- rel_otu_table[order(rownames(rel_otu_table)),]
    
    out <- list("otu_file"= rel_otu_table, "meta_file"= meta_file, "taxa_file"= taxa_file)
    
    return(out)
}

# Calculate average and select the higher values
  average.tab <- function(Data4){
    #Calculate the average 
    df.average <- aggregate(Data4[,2:length(Data4)],Data4["samples"],mean) 
  
    #Select most abundance OTUs
    df.final <- as.data.frame(t(df.average[,2:length(df.average)]))
    colnames(df.final) <- df.average[,1]
    df.final <- data.matrix(df.final)
    return(df.final)
  }

# Calculate Standard error
  serror.tab <- function(Data4, average_table, list_fac){
    #Standard Error
    sd <- aggregate(Data4[,rownames(average_table)],Data4["samples"],sd)
    se <- sd[,-1]/sqrt(max(sort(table(list_fac))))
    se <- as.data.frame(t(se))
    colnames(se) <- sd$samples
  return(se)
  }

# Combine average and SE data
  barchart.data <- function(physeq, group_name, change_fac, ordered_fac, n, level_taxa="Genus" ){
    #Load Data
    dataset <- abundance_table(physeq, level_taxa=level_taxa)
    otu_file<- dataset$otu_file
    meta_file <- dataset$meta_file
    taxa_file <- dataset$taxa_file
    
    #Create table with samples and abundance table
    Data <- data.frame("samples" = meta_file[, group_name], t(otu_file))
    
    #Calculate average based on the selected factor
    average_table <- average.tab(Data)
    sums<-  rowSums(average_table)
    average_table <- data.frame(cbind(average_table, "sums"=sums), check.names = F)
    average_table_order<-average_table[order(average_table$sums),]
    average_table <- tail(average_table_order, n)
    
    #Calculate standard error
    se_table <- serror.tab(Data, average_table, list_fac = meta_file[, group_name])
    
    if(change_fac == T){
      #change factor Order
      level_position <- match(ordered_fac, levels(as.factor(meta_file[,group_name])))
      new_levels <- levels(sort.levels(as.factor(meta_file[,group_name]), level_position))
      #Order columns in average table
      average_table <- average_table[,c(new_levels, "sums")]
      se_table <- se_table[,c(new_levels)]
    }
    
    out <- list("average_table"= average_table, "se_table"= se_table)
    return(out)
}

# Average table for stacked bar plot and heatmap
  data.plot <- function(physeq, level_taxa, group_name, change_fac, ordered_fac, threshold ){
  
    #Calculate relative abundance table and average per group
    rb_tab <- abundance_table(physeq, level_taxa = level_taxa)
    otu_file<- rb_tab$otu_file
    meta_file <- rb_tab$meta_file
    taxa_file <- rb_tab$taxa_file
    
    #combine samples metadata and count matrix (otu_file)
    Data <- data.frame("samples" = meta_file[, group_name], t(otu_file))
    
    colnames(Data) <- c("samples", rownames(otu_file))
    
    #Calculate average per samples group
    average <- data.frame(average.tab(Data), check.names = F)
    
    #calculate average of otus across all the samples
    average$AVG <- apply(average[1:dim(average)[2]], 1, mean)
    average <- average[order(average$AVG),]
    
    #####Select most abundant OTUs and gather the others
    if(threshold != 0){
      #1- Extract data based on the average relative abundance of OTUs 
      top_avg <- average %>%
        filter(AVG >= threshold)
      bot_avg <- subset(average , !(rownames(average)  %in% rownames(top_avg)))
      bot_avg <- colSums(bot_avg)
      
      #combine the top abundant and the gathered Otus
      AVG_tab <- rbind(top_avg, Others = bot_avg)
      rownames(AVG_tab)[length(rownames(AVG_tab))] <- paste("Others (av < ",threshold,"%)", sep = "")
      
      #Calculate standard error
      se_table <- serror.tab(Data, AVG_tab[1:(dim(AVG_tab)[1]-1), ], list_fac = meta_file[, group_name])
    } else {
      AVG_tab <- average
      se_table <- serror.tab(Data, AVG_tab, list_fac = meta_file[, group_name])
    }
    
    if(change_fac == T){
      #change factor Order
      level_position <- match(ordered_fac, levels(as.factor(meta_file[,group_name])))
      new_levels <- levels(sort.levels(as.factor(meta_file[,group_name]), level_position))
      #Order columns in average and standard error tables
      AVG_tab <- AVG_tab[,c(new_levels, "AVG")]
      se_table <- se_table[, c(new_levels)]
    }
    
    AVG_tab <- cbind("OTUs" = rownames(AVG_tab), AVG_tab)
    se_table <- cbind("OTUs" = rownames(se_table), se_table)
    out <- list("AVG_tab"= AVG_tab, "SE_tab"= se_table)
    return(out)
  }

# Relative abundance Plot
  st.barchart.plot.main <- function(AVG_tab){
  
  data5 <- melt(AVG_tab[,-dim(AVG_tab)[2]], id.vars = "OTUs")
  
  ggplot(data5, aes(x = variable, y = value, fill= OTUs)) +
    geom_bar(stat = "identity", position="stack", width = 0.8, colour="black", size = 0.2)+
    theme_classic()+
    scale_y_continuous(name = "Relative abundance (%)", expand = c(0,0))+
    scale_x_discrete(name = "")+
    #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
          axis.text = element_text(size=12, color = "black"),
          legend.text = element_text(size=11))+
    theme( panel.grid.major.x = element_blank() ,
           panel.grid.major.y = element_line( size=.5, color="grey" ),
           legend.position = "bottom", legend.title = element_blank()
           )
    #scale_fill_viridis(discrete = TRUE, option = "H", direction = -1)
}

  barchart.plot.main <- function(AVG_tab, upper, lower){
  data <- data.frame("OTUs" = rownames(AVG_tab), AVG_tab, row.names = NULL, check.names = F)
  new_data1 <- melt(data = data, id.vars = "OTUs")
  upper2  <- data.frame("OTUs" = rownames(upper), upper, row.names = NULL)
  upper2 <- melt(data = upper2, id.vars = "OTUs")
  lower2  <- data.frame("OTUs" = rownames(lower), lower, row.names = NULL)
  lower2 <- melt(data = lower2, id.vars = "OTUs")
  new_data <- data.frame(new_data1, "upper" = upper2$value, "lower" = lower2$value)
  colnames(new_data)<- c("OTUs", "variable", "Abundance", "upper", "lower" )
  
  
  ggplot(new_data, aes(x = variable, y = Abundance, fill= OTUs)) + 
    geom_bar(stat = "identity", position=position_dodge(), width = 0.6, colour="black", size = 0.2)+ 
    theme_classic()+
    scale_y_continuous(name = "Relative abundance (%)", expand = c(0,0))+
    scale_x_discrete(name = "")+
    geom_errorbar(aes(ymin=lower, ymax=upper ), size = 0.2, width=0.4, position=position_dodge(.6))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(plot.margin = unit(c(0,01,0,0), "cm"))+
    theme( panel.grid.major.x = element_blank() ,panel.grid.major.y = element_line( size=.5, color="grey" ))+
    scale_fill_brewer(palette="Paired")

}
