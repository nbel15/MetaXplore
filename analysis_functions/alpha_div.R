# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#############################
# import requested functions
#############################
source("tools_functions.R")

# Observed OTUs = number of taxa observed in a sample
  obs.otus <- function(df){
    otus.obs <- length(which(df > 0))
    return(otus.obs)
  }

#Calculate Alpha diversity indices
  alpha_div <- function(physeq,method){
    #==check for validity of selected methods
    method<- match.arg(method,c("richness", "fisher", "simpson", "shannon", "evenness","pd", "goods", "ace", "chao", "obs.otus", "simpson_eff", "shannon_eff"), several.ok = TRUE)
    
    abund_table <- otu_table(physeq)
    df <- NULL
    if("richness"%in%method){
      R<- vegan::rarefy(abund_table,min(rowSums(abund_table)))
      df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))
      if(is.null(df)){
        df<-df_R}
      else {
        df<-rbind(df,df_R)}
    }
    if("fisher"%in%method){
      alpha <- vegan::fisher.alpha(abund_table)
      df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))
      if(is.null(df)){
        df<-df_alpha}
      else {
        df<-rbind(df,df_alpha)}
    }
    if("simpson"%in%method){
      simp <- vegan::diversity(abund_table, "simpson")
      df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))
      if(is.null(df)){
        df<-df_simp}
      else {
        df<-rbind(df,df_simp)}
    }
    if("shannon"%in%method){
      H<- vegan::diversity(abund_table)
      df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))
      if(is.null(df)){
        df<-df_H}
      else {
        df<-rbind(df,df_H)}
    }
    if("evenness"%in%method){
      H<-vegan::diversity(abund_table)
      S <- specnumber(abund_table)
      J <- H/log(S)
      df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))
      if(is.null(df)){
        df<-df_J}
      else {
        df<-rbind(df,df_J)}
    }
    if("pd"%in%method){
      otu_tree <- phyloseq::phy_tree(physeq)
      PD <- pd(as.matrix(abund_table@.Data), otu_tree, include.root = TRUE)
      df_PD<-data.frame(sample = rownames(PD),value = PD$PD, measure = rep("PD",dim(PD)[1]))
      if(is.null(df)){
        df<-df_PD}
      else {
        df<-rbind(df,df_PD)
      }
    }
    if("goods"%in%method){
      GC<-QsRutils::goods(abund_table)
      df_C<-data.frame(sample=rownames(GC),value=GC$goods,measure=rep("Goods Coverage",dim(GC)[1]))
      if(is.null(df)){
        df<-df_C}
      else {
        df<-rbind(df,df_C)}
    }
    if("ace"%in%method){
      ACE.tab <- apply(abund_table, MARGIN = 1, FUN = fossil::ACE)
      df_ACE<-data.frame(sample=names(ACE.tab),value=ACE.tab,measure=rep("ACE",length(ACE.tab)))
      if(is.null(df)){
        df<-df_ACE}
      else {
        df<-rbind(df,df_ACE)}
    }
    if("chao"%in%method){
      Chao.tab <- apply(abund_table, MARGIN = 1, FUN = fossil::chao1)
      df_Chao<-data.frame(sample=names(Chao.tab),value=Chao.tab,measure=rep("Chao1",length(Chao.tab)))
      if(is.null(df)){
        df<-df_Chao}
      else {
        df<-rbind(df,df_Chao)}
    }
    if("obs.otus"%in%method){
      obs.otus.tab <- apply(abund_table, MARGIN = 1, FUN = obs.otus)
      df_obs.otus<-data.frame(sample=names(obs.otus.tab),value=obs.otus.tab,measure=rep("Observed OTUs",length(obs.otus.tab)))
      if(is.null(df)){
        df<-df_obs.otus}
      else {
        df<-rbind(df,df_obs.otus)}
    }
    if("simpson_eff"%in%method){
      simp_eff <- vegan::diversity(abund_table, "invsimpson")
      #simp_eff <- 1/(1-simp_eff)
      df_simp_eff<-data.frame(sample=names(simp_eff),value=simp_eff,measure=rep("Simpson_eff",length(simp_eff)))
      if(is.null(df)){
        df<-df_simp_eff}
      else {
        df<-rbind(df,df_simp_eff)}
    }
    if("shannon_eff"%in%method){
      H<- vegan::diversity(abund_table)
      H <- exp(H)
      df_H_eff<-data.frame(sample=names(H),value=H,measure=rep("Shannon_eff",length(H)))
      if(is.null(df)){
        df<-df_H_eff}
      else {
        df<-rbind(df,df_H_eff)}
    }
    
    return(df)
  }

#Generate alpha diversity table
  alpha.div.tab<- function(physeq, method, grouping_column){
    #enforce orientation
    if(taxa_are_rows(physeq)){
      physeq <- t(physeq)
    }
    abund_table <- otu_table(physeq)
    meta_table <- sample_data(physeq)
    
    #get diversity measure using selected methods
    div.df <- alpha_div(physeq,method)
    
    #=add grouping information to alpha diversity measures
    df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
    df
  }

#Generate alpha diversity boxplot
# Modified from alpha diversity in microbiomSeq package (umerijaz/microbiomeSeq)  
  alpha.div.plot <- function(physeq, method, grouping_column, pValueCutoff=0.05, neword=NULL, changefac=T){
    #enforce orientation
    if(taxa_are_rows(physeq)){
      physeq <- t(physeq)
    }
    abund_table <- otu_table(physeq)
    meta_table <- sample_data(physeq)
    
    #get diversity measure using selected methods
    div.df <- alpha_div(physeq, method)
    
    #Exclude Goods coverage from the plot
    div.df <- subset(div.df, measure!="Goods Coverage")
  
    
    #=add grouping information to alpha diversity measures
    df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
    
    #change the facto's order of the targeted column
      if (changefac == T){
      df <- order.levels.df(df = df,factor = grouping_column, orderd_factor1 = neword, changefac = changefac)
      }
      
    #perform anova of diversity measure between groups if the group contain more than 1 element
    df_pw <-NULL
    if(length(unique(as.character(meta_table[[grouping_column]])))>1){
      anova_res <- perform_anova2(df,meta_table,grouping_column,pValueCutoff)
      df_pw <- data.frame(anova_res$df_pw) #get pairwise p-values
      #df_pw$measure <- factor(df_pw$measure, levels = c("Richness", "Simpson", "Shannon", "Pielou's evenness", "PD"))
      }
      
    #Order indices
     #df$measure <- factor(df$measure, levels = c("Richness", "Simpson", "Shannon", "Pielou's evenness", "PD"))
    
      #Draw the boxplots
      p <- ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)+
            geom_boxplot()+
            geom_jitter(position = position_jitter(height = 0, width=0))+
            theme_bw()+
            facet_wrap(~ measure, scales="free",nrow=1)+
            ylab("Observed Values")+
            theme(strip.background = element_rect(fill = "white"), panel.spacing = unit(0.005, "cm", data = NULL))+
            xlab("")
      
      #This loop will generate the lines and significance
      if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
        df[,grouping_column] <- as.factor(df[,grouping_column])
         for(i in 1:dim(df_pw)[1]){
          p <- p + geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
          p <- p + geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
          #p<-p + geom_text(inherit.aes=F,aes(x=x,y=y,label=p.value),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=(as.numeric(as.character(df_pw[i,"y"])))+((as.numeric(as.character(df_pw[i,"y"])))*0.025),measure=as.character(df_pw[i,"measure"]),p.value=as.character(df_pw[i,"p"])))
          
        }
      }
      return(p)
  }

#Perform statistical ANOVA test
# Modified from alpha diversity in microbiomSeq package (umerijaz/microbiomeSeq)  
  perform_anova2 <- function(df, meta_table, grouping_column, pValueCutoff){
    temp <- data.frame(df, meta_table[, grouping_column])
    colnames(temp) <- c(colnames(df), "group")
    dt<-data.table::data.table(temp)
    
    #specifying a p-value cutoff for the ggplot2 strips
    pval<-dt[, list(pvalue = sprintf("%.2g",tryCatch(summary(aov(value ~ group))[[1]][["Pr(>F)"]][1],error=function(e) NULL))),by=list(measure)]
  
    #Filter out pvals that we are not significant
    pval<-pval[!pval$pvalue == "",]
    pval<-pval[as.numeric(pval$pvalue) <= 1,]
    
    #using sapply to generate significances for pval$pvalue using the cut function.
    pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, pValueCutoff, Inf),label=c("***", "**", "*", "")))})
    
    #Update df$measure to change the measure names if the grouping_column has more than three classes
    if(length(unique(as.character(meta_table[,grouping_column])))>2){
      df$measure<-as.character(df$measure)
      if(dim(pval)[1]>0){
        for(i in seq(1:dim(pval)[1])){
          df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
        }
      }
      df$measure<-as.factor(df$measure)
    }
    #Get all possible pairwise combination of values in the grouping_column
    s<-combn(unique(as.character(df[,grouping_column])),2)
    
    #df_pw will store the pair-wise p-values
    df_pw<-NULL
    for(k in unique(as.character(df$measure))){
      #We need to calculate the coordinate to draw pair-wise significance lines
      #for this we calculate bas as the maximum value
      bas<-max(df[(df$measure==k),"value"])
      #Calculate increments as 10% of the maximum values
      inc<-0.1*bas
      #Give an initial increment
      bas<-bas+inc
      for(l in 1:dim(s)[2]){
        #Do a pair-wise anova
        tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
        #Ignore if anova fails
        if(!is.na(as.numeric(tmp[length(tmp)]))){
          #Only retain those pairs where the p-values are significant
          if(as.numeric(tmp[length(tmp)])<pValueCutoff){
            if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}
            #Generate the next position
            bas<-bas+inc
          }
        }
      }
    }
    if(!is.null(df_pw)){
      #df_pw2 <- data.frame() 
      df_pw <- rbind(df_pw )
      colnames(df_pw)<-c("measure","from","to","y","p")
      
      
    }
    out <- list("df_pw"=df_pw, "df"=df)
    return(out)
  }

#Generate alpha diversity table front-end (with Standard Error)
  alpha.div.tab.se <- function(physeq, method, grouping_column, add.se){
    #calculate alpha diversity indices
    alpha <- alpha.div.tab(physeq = physeq, 
                           method = method,
                           grouping_column = grouping_column
    )
    
    #rename calculated table
    colnames(alpha) <- c("SamplesID", "Value", "Indices", "Factor") 
    
    #Aggregate by mean alpha diversity indices based on the selected factor per index
    alpha.agg.mean <- aggregate(alpha$Value, by = list(alpha$Indices, alpha$Factor), FUN = mean)
    colnames(alpha.agg.mean) <- c("indices", "Factor", "indices_mean")
    alpha.agg.mean$indices_mean <- round(alpha.agg.mean$indices_mean, 2)
    
    #Spread aggregated indices (factor in rows and indices in column)
    alpha.mean.sp <- spread(alpha.agg.mean, "indices", "indices_mean")
    
    if(add.se == T){
      
      meta_table <- sample_data(physeq)
      #Calculate standard error based on standard deviation values
      alpha.agg.sd <- aggregate(alpha$Value, by = list(alpha$Indices, alpha$Factor), FUN = sd)
      alpha.agg.se <- alpha.agg.sd
      alpha.agg.se$x <-alpha.agg.se$x/sqrt(max(sort(table(meta_table[,grouping_column]))))
      colnames(alpha.agg.se) <- c("indices", "Factor", "indices_se")
      alpha.agg.se$indices_se <- round(alpha.agg.se$indices_se, 2)
      
      #Spread aggregated standard error (factor in rows and indices in column)
      alpha.se.sp <- spread(alpha.agg.se, "indices", "indices_se")
      
      #Select Standard error columns
      tab.indices <- as.data.frame(alpha.mean.sp[,2:length(alpha.mean.sp)])
      colnames(tab.indices) <- colnames(alpha.mean.sp)[2:length(alpha.mean.sp)]
      
      #Select indices columns
      se.indices <- as.data.frame(alpha.se.sp[,2:length(alpha.se.sp)])
      colnames(se.indices) <- colnames(alpha.se.sp)[2:length(alpha.se.sp)]
      
      #Merge aggregated indices and standard error
      alpha.tab <- merge.tab.se.fun(tab = tab.indices,se.tab = se.indices, add.se = add.se)
      #alpha.tab <- merge.tab.se.fun(tab = alpha.mean.sp[,2:length(alpha.mean.sp)], se.tab = alpha.se.sp[,2:length(alpha.se.sp)], add.se = T)
    } else {
      alpha.tab <- as.data.frame(alpha.mean.sp[,2:length(alpha.mean.sp)])
      colnames(alpha.tab) <- colnames(alpha.mean.sp)[2:length(alpha.mean.sp)]
      
    }
    
    rownames(alpha.tab) <- alpha.mean.sp$Factor
    
    return(alpha.tab)
  }
