# MetaXplore Version 1.0
# Last modified on 05/03/2024
# Author: Naima BEL MOKHTAR (Email: naima1503@gmail.com)
########################################################

#########################################
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(MicrobiotaProcess)
library(tidytree)
library(ggstar)
library(forcats)


# Calculate MPSE abundance
biomarker.abundance <- function(physeq, group.fac){
  # convert phyloseq object to mpse
  mpse <- physeq %>% as.MPSE() 
  
  #calculate abundance
  mpse %<>%
    mp_cal_abundance( # for each samples
      .abundance = RareAbundance
    ) %>%
    mp_cal_abundance( # for each groups 
      .abundance=RareAbundance,
      .group=!!rlang::sym(group.fac)
    )
  return(mpse)
}

# Perform diffential analysis 
biomarker.diff <- function(physeq, group.fac="Age2", lda.threshold=2, level.taxa="Genus", strict=T){
 #biomarker discovery based on tidy-like framework
  mpse <- biomarker.abundance(physeq, group.fac)
  
  mpse %<>%
    mp_diff_analysis(
      .abundance = RelRareAbundanceBySample,
      .group=!!rlang::sym(group.fac),
      first.test.alpha = 0.01,
      tip.level = !!rlang::sym(level.taxa), relative = T, strict = strict, ldascore = lda.threshold
    )
  return(mpse)
}

biomarker.plot.fun <- function(physeq, group.fac, rb.factor, lda.threshold, level.taxa, strict){

  mpse4 <- biomarker.diff(physeq, group.fac, lda.threshold, level.taxa, strict)
  # The result is stored to the taxatree slot.
  taxa.tree <- mpse4 %>% 
    mp_extract_tree(type="taxatree")
  
  
  # display the high light of phylum clade.
  highlight.level = "Phylum"
  
  # extract number of samples
  nsample <- mpse4 %>% mp_extract_sample() %>% nrow()
  
  #rb.fac <- ifelse(rb.fact == "Factor", group.fac, "Sample")
  #nfactor <- length(table(physeq@sam_data[[rb.fact]]))
  
  # extract selected fiels in taxa.tree object
  field.da.nm <- tidytree::get.fields(taxa.tree)
  # extract selected factor name and highled the different labels based on the selected factor
  sign.field <- field.da.nm[grepl('^Sign_', field.da.nm)][1]
  
  # Check if any sign result exist
  if(is.na(sign.field))
  {
    group.nm <- group.fac #Age2
    flag <- grepl(paste0("By", group.nm), field.da.nm)
  }else{
    group.nm <- gsub('Sign_', "", sign.field) #Age2
    flag <- grepl(paste0("By", group.nm), field.da.nm)
  }
  
  if (rb.factor == F){  #nsample != nfactor
    if (!any(flag)){
      stop_wrap("The relative abundance of each group will be displayed, but the
                       relative abundance of each group is not calculated, please run 
                       the mp_cal_abundance specified group argument before !")
    }
    abun.col <- field.da.nm[flag]
  }else{
    abun.col <- field.da.nm[grepl("BySample", field.da.nm)]
  }
  abun.col <- abun.col[1]
  #abun.col <- field.da.nm[flag] #RareAbundanceByAge2
  
  x.abun.col <- taxa.tree %>% 
    dplyr::select(!!rlang::sym(abun.col)) %>%
    tidyr::unnest(!!rlang::sym(abun.col)) %>%
    colnames()
  if (any(grepl('BySample$',  abun.col))){
    if (any(grepl('^Rel', x.abun.col))){
      x.abun.col <- paste0('Rel', abun.col)
    }else{
      x.abun.col <- gsub('BySample$', '', abun.col)
    }
  }else{
    if (any(grepl("^Rel", x.abun.col))){
      x.abun.col <- paste0('Rel', abun.col)
    }else{
      x.abun.col <- abun.col
    }
  }
  #x.abun.col <- paste0('Rel', abun.col) #RelRareAbundanceByAge2
  
  
  if (rb.factor == F){  #nsample != nfactor
    mapping <- aes(x=!!rlang::sym(group.nm), size = !!rlang::sym(x.abun.col), 
                   fill = !!rlang::sym(group.nm), subset = !!rlang::sym(x.abun.col) > 0) 
    n.pwidth <- mpse4 %>%
      mp_extract_sample() %>%
      dplyr::pull(!!rlang::sym(group.nm)) %>%
      unique() %>% length()    
  }else{
    mapping <- aes(x=forcats::fct_reorder(!!rlang::sym("Sample"), !!rlang::sym(group.nm), .fun=min),
                   size = !!rlang::sym(x.abun.col), 
                   fill = !!rlang::sym(group.nm), 
                   subset = !!rlang::sym(x.abun.col) > 0
    )
    n.pwidth <- ncol(mpse4)
  }
  
  # Since taxa.tree is treedata object, it can be visualized by ggtree and ggtreeExtra
  p1 <- ggtree(taxa.tree, layout="rectangular", size = 0.3) +
    geom_point(data = td_filter(!isTip), fill="white", size=1, shape=21)
  
  p2 <- p1 +
    geom_hilight(
      data = td_filter(nodeClass == highlight.level),
      mapping = aes(node = node, fill = label)#, type = "encircle"
    )
  
  # display the relative abundance of features(OTU)
  p3 <- p2 +
    ggnewscale::new_scale("fill") +
    geom_fruit(
      data = td_unnest(abun.col),
      geom = geom_star,
      mapping = mapping, #aes(x = !!rlang::sym(group.nm),size = !!rlang::sym(x.abun.col),fill = !!rlang::sym(group.nm),subset = !!rlang::sym(x.abun.col) > 0),
      starshape = 13,
      starstroke = 1,
      offset = 0.02,
      pwidth = n.pwidth*0.02, 
      #axis.params = aes(axis ="x", text.size= 3, text.angle = 90, hjust=1),
      grid.params = list(linetype=2)
    ) +
    scale_size_continuous(
      name="Relative Abundance (%)",
      range = c(.3, 3)) +
    scale_fill_viridis(discrete = TRUE, option = "D", direction = -1)
  #+    scale_fill_manual(values=c("#1B9E77", "#D95F02", "blue", "red"))
  
  # display the tip labels of taxa tree
  p4 <- p3 + 
    geom_tiplab(size=5, offset = max(p3$data$xmaxtmp, na.rm=TRUE) - 0.95*max(p3$data$x, na.rm=TRUE)) #7.2
  # display the LDA of significant OTU.
  title.height <- 4.4e-06 * sum(p4$data$isTip)
  
  if(!is.na(sign.field))
  { 
  p5 <- p4 +
    ggnewscale::new_scale("fill") +
    geom_fruit(
      #data = td_filter(!is.na(!!rlang::sym(sign.field))),
      geom = geom_col,
      mapping = aes(
        x = LDAmean,
        fill = !!rlang::sym(sign.field)
      ),
      orientation = "y",
      offset = 0.4,
      pwidth = 0.2,
      axis.params = list(axis = "x",
                         title = "Log10(LDA)",
                         title.height = title.height,
                         title.size = 4,
                         text.size = 3.5,
                         vjust = 1),
      grid.params = list(linetype = 2)
    )
  
  # display the significant (FDR) taxonomy after kruskal.test (default)
  p.final <- p5 +
    ggnewscale::new_scale("size") +
    geom_point(
      data=td_filter(!is.na(!!rlang::sym(sign.field))),
      mapping = aes(size = -log10(fdr),
                    fill = !!rlang::sym(sign.field),
      ),
      shape = 21, 
    ) +
    scale_size_continuous(range=c(1, 4)) +
    scale_fill_viridis(discrete = TRUE, option = "D", direction = -1)+
    theme(legend.text=element_text(size=14))
  
  #theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))
  }else {
    p.final <- p4 + ggplot2::xlim(0, 8)
}
  
  return(p.final)
  
}