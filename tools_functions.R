# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# create taxonomy matrix
#########################################
split.taxon <- function(taxamat){
  list_taxon <- as.vector(taxamat$Taxon)
  split_taxon <- strsplit(list_taxon, ";") 
  taxa_tab <- data.frame(do.call(rbind.data.frame, split_taxon))
  if(dim(taxa_tab)[2] >= 7){
    taxa_tab <- taxa_tab[,1:7]
    colnames(taxa_tab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  }
  if(dim(taxa_tab)[2] == 6){
    taxa_tab <- taxa_tab[,1:6]
    colnames(taxa_tab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  }
  
  rownames(taxa_tab) <- taxamat[,1]
  
  return(taxa_tab)
}

#########################################
# Select factors
#########################################
change.factor.order <- function(metamat, factor){
  dat <- as.factor(metamat[, factor])
  factor_var <- levels(dat)
  return(factor_var)
}

#########################################
# Sort factor levels arbitrarily
#########################################
sort.levels <- function(old.factor, level.Order) {
  if(length(level.Order) == length(levels(old.factor))) {
    reorderedFactor <- factor(old.factor, levels = levels(old.factor)[level.Order])
  }
  
  if(length(level.Order) < length(levels(old.factor))) {
    levelOrderAll <- c(level.Order, (1:length(levels(old.factor)))[-level.Order])
    reorderedFactor <- factor(old.factor, levels = levels(old.factor)[levelOrderAll])
  }
  
  return(reorderedFactor)
}

#########################################
# Select factors based on physeq object 
#########################################
change.factor.order <- function(physeq, factor.tab){
  metamat <- sample_data(physeq)
  dat <- as.factor(metamat[[factor.tab]])
  factor_var <- levels(dat)
  return(factor_var)
}

#############################################
# order factors in meta data -alpha div only-
#############################################
order.levels.df <- function(df, factor, orderd_factor1, changefac){
  if (changefac == T){
    Ids <- as.factor(df[,c(factor)])
    level_position <- match(orderd_factor1, levels(Ids))
    Ids <- sort.levels(Ids, level_position)
    df2<- cbind(df[, -which(names(df) %in% factor)], Ids)
    colnames(df2) <- colnames(df) 
    return(as.data.frame(df2))}
  else return(df)
  
}

###############################################
# Combine Relative abundance and Standard error
###############################################
merge.tab.se.fun <- function(tab, se.tab, add.se){
  
  if(add.se == T){
    var <- colnames(tab) 
    #se.tab <- se.tab[,var]
    merged.tab <-list()
    for(i in 1:length(var)){
      merged.tab[[i]] <- paste(tab[,var[i]], se.tab[,var[i]], sep = "\u00B1")
      }
    merged.tab <- as.data.frame(merged.tab, row.names = rownames(tab))
    colnames(merged.tab) <- colnames(tab)
    tab.out <- merged.tab
  }
  else
    tab.out <- tab
  
  return(tab.out)  
}

###############################################
# Filter out otus with low abundance
###############################################
otus.abund.filt <- function(physeq, abund.threshold){
  if(!taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  ##Store input data
  otu_file <- phyloseq::otu_table(physeq)
  otu_file <- data.frame(otu_file@.Data)
  otu_file <- data.frame(otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),])
  otu_file <- otu_file[sort(rownames(otu_file)), ]
  
  rel.otu.tab.all <- decostand(otu_file, method = "total", MARGIN = 2)*100
  rel.otu.tab.all <- rel.otu.tab.all[order(rownames(rel.otu.tab.all)),]
  
  taxa.sum <- rowMeans(rel.otu.tab.all)
  taxa.sum <- taxa.sum[taxa.sum >= abund.threshold]
  
  physeq.filt <- prune_taxa(names(taxa.sum),physeq)
  
  return(physeq.filt)
  
}

###############################################
# filter by frequency input
###############################################
  otus.freq.filter <- function(physeq, min.freq=0.001){
    otutab <- phyloseq::otu_table(physeq)
    freq.tab <- as.data.frame(otutab@.Data)
    total_counts <- colSums(otutab)
    for (col in colnames(otutab)) {
      freq.tab[, col] <- otutab[, col] / total_counts[col]
    }
    #min.freq <- 0.001
    filtered_freq.tab <- freq.tab %>% filter_all(any_vars(is.numeric(.) & . >= min.freq))
    
    if(nrow(filtered_freq.tab) != 0){
      filtered.otutab <- otutab[rownames(filtered_freq.tab),]
      phyloseq::otu_table(physeq) <- phyloseq::otu_table(filtered.otutab, taxa_are_rows = phyloseq::taxa_are_rows(physeq))
      return(physeq)
    }else{
      return(NULL)
    }
  }
  
###############################################
# Normalize input
###############################################
# from joey711/shiny-phyloseq/
  gm_mean <- function(x, na.rm = TRUE) {
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm = na.rm) / length(x))
  }

  trans_clr <- function(x, base = exp(1)) {
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
} 

  otus.norm <- function(physeq, method){
    if(method == "otunorm"){
      otutab <- phyloseq::otu_table(physeq)
      size <- min(phyloseq::sample_sums(physeq))
      # Apply normalization to each column
      otu.norm <- as.data.frame(apply(otutab, MARGIN = 2, FUN =  function(x) round((x * 394) / sum(x))))
      phyloseq::otu_table(physeq) <- phyloseq::otu_table(otu.norm, taxa_are_rows = phyloseq::taxa_are_rows(physeq))
      physeq.norm = physeq
    }
    
  if(method == "Rarefy"){
    size <- min(phyloseq::sample_sums(physeq))
    physeq_rarefied <- phyloseq::rarefy_even_depth(physeq, rngseed=1, replace=T, sample.size = size)
    physeq.norm = physeq_rarefied
    
  }
  if(method == "TSS"){
    otutab <- phyloseq::otu_table(physeq)
    size <- colSums(otutab)
    otu.norm <- sweep(otutab, MARGIN = 2, STATS = size, FUN = "/")
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(otu.norm, taxa_are_rows = phyloseq::taxa_are_rows(physeq))
    physeq.norm = physeq
  }
  if(method == "CLR"){
    otu <- as(phyloseq::otu_table(physeq), "matrix")
    otu.norm <- apply(otu, 2, trans_clr)
    
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(otu.norm, taxa_are_rows = phyloseq::taxa_are_rows(physeq))
    physeq.norm = physeq
  }
  if(method == "sqrt"){
    otutab <- phyloseq::otu_table(physeq)
    #rb.otutab <- decostand(otutab, method = "total", MARGIN = 2)
    #otu.norm <- decostand(rb.otutab, method = "hellinger", MARGIN = 2)
    otu.norm <- sqrt(otutab)
    
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(otu.norm, taxa_are_rows = phyloseq::taxa_are_rows(physeq))
    physeq.norm = physeq
  }
  if(method == "log"){
    otutab <- phyloseq::otu_table(physeq)
    otu.norm <- log2(otutab + 1)
    
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(otu.norm, taxa_are_rows = phyloseq::taxa_are_rows(physeq))
    physeq.norm = physeq
  }
  if(method == "None"){
    physeq.norm = physeq    
  }
  
  return(physeq.norm)
}
