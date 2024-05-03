# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
# Import requested functions
source("tools_functions.R")

core.tab <- function(physeq, prevalence, rb.threshold){
  #phydata.core.all <- microbiome::core(physeq, detection = rb/100, prevalence = prevalence/100)
  
  #core.taxa <- taxa(phydata.core.all)
  #core.taxa.mat <- phyloseq::tax_table(phydata.core.all)
  #core.taxa.df <- as.data.frame(core.taxa.mat)
  #core.taxa.df$Genus2 <- paste(rownames(core.taxa.df), core.taxa.df$Genus, sep ="_" )
  
  #core.abundance <- sample_sums(phydata.core.all)
  
  #core.rel<- core.abundance*100/sample_sums(physeq)
  #df <- mean(core.rel)
  
  
  otu_file <- phyloseq::otu_table(physeq)
  otu_file <- as.data.frame(otu_file@.Data)
  otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
  
  otu.file.rel <- decostand(otu_file, method = "total", MARGIN = 2)*100
  
  #otu.rel.core <- otu.file.rel[core.taxa,]
  otu.rel.core.temp <- otu.file.rel %>% 
                    mutate(mean = rowMeans(.)) %>% 
                    filter(mean >= rb.threshold) %>%
                    select(-mean)
  
  
  
  otu.rel.core.temp$Prevalence <- (apply(otu.rel.core.temp, MARGIN = 1, function(x) sum(x>0)*100))/dim(otu.rel.core.temp)[2]
  otu.rel.core <- otu.rel.core.temp %>% 
    filter(Prevalence >= prevalence) %>%
    select(-Prevalence)
  
  
  #df <- data.frame("Genus"=core.taxa.df[rownames(otu.rel.core),"Genus"], "RB"=rowMeans(otu.rel.core)*100, "Prevalence"=Prevalence  )
  df <- data.frame(otu.rel.core)

  return(df)
}

hashbind <- function(..., hash = "#") {
  lst <- list(...)
  Nchar <- max(rapply(lst, function(y) nchar(as.character(y)))) + 2
  do.call(
    rbind, 
    lapply(lst, function(x) {
      rbind(x, substr(paste(rep(hash, Nchar), collapse = ""), 1, Nchar))
    }))
}

subset.samples.safe = function(physeq, var, value) {
  keeps = rownames(phyloseq::sample_data(physeq)[phyloseq::sample_data(physeq)[[var]] == value,])
  phyloseq::prune_samples(keeps, physeq)
}

core.microbiome <- function(physeq, core.factor, prevalence, rb.threshold){
  
  taxa_file <- phyloseq::tax_table(physeq)
  taxa_file <- data.frame(taxa_file@.Data)
  taxa_file <- taxa_file[sort(rownames(taxa_file)), ]
  
  meta_file <- phyloseq::sample_data(physeq)
  meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),],row.names=row.names(meta_file))
  meta_file <- data.frame(meta_file[order(row.names(meta_file)),], stringsAsFactors = T)
  meta_file[, core.factor] <- as.factor(meta_file[, core.factor])
  
  groups <- levels(meta_file[, core.factor])
  
  list.core <- list()
  list.core.taxa <- list()
  library(mattsUtils)
  for(group in 1:length(groups)){
    #physeq.filtered <- subset_samples(physeq = physeq , Generation == groups[group])
    physeq.filtered <- subset.samples.safe(physeq, core.factor, groups[group])
    core.temp <- core.tab(physeq = physeq.filtered, prevalence=prevalence, rb.threshold = rb.threshold)
    list.core[[group]] <- core.temp
    list.core.taxa[[group]] <- rownames(core.temp)
  
  }
  names(list.core)<- groups
  names(list.core.taxa)<- groups
  
  library(qdapTools)
  core.combined <- t(mtabulate(list.core.taxa))
  core.combined.taxa <- data.frame(taxa_file[rownames(core.combined),], core.combined, check.names = F)
  
  output = list("core.combined.taxa" = core.combined.taxa, "list.core.taxa" = list.core.taxa)
  
  return(output)
}