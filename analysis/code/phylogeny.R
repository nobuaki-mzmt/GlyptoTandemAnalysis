# Pylogenetic analysis
{
  options(warn = 0)
  rm(list = ls())
  
  library(dplyr)
  
  require(phytools)
  require(stringr)
  
  library(ggplot2)
  library(viridis)
  library(scales)
}

# ---------------------------------------------------------------------------- #
# prep data
# ---------------------------------------------------------------------------- #
{
  # tree data
  {
    tree.file <- "data_raw/tree/run5_400M_burnin20_mcc_median.tree"
    
    tree <- read.nexus(tree.file)
    labels <- tree$tip.label
    labels <- str_replace(labels, "Cubi_tenu", "Cubitermes_tenuiceps")
    labels <- str_replace(labels, "Tern_pall", "Ternicubitermes_pall")
    labels <- str_replace(labels, "Tern_pall", "Ternicubitermes_pall")
    # as for Anacanthotermes, the tree is from sp. 
    # The mating system information is available for macrocephalus, so we use macrocephalus for reconstruction.
    labels <- str_replace(labels, "Anacanthotermes_sp", "Anacanthotermes_macrocephalus")
  }

  # tandem data (remove some taxa without phylogeny information)
  {
    d_tandem <- read.csv("data_raw/tree/tandem_info_Mizumoto-etal-2022-PNAS.csv")[,1:6]
    d_mate   <- read.csv("data_raw/tree/mating_system.csv")[,1:7]
    d_tandem <- full_join(d_tandem, d_mate, by = c("Group", "Genus", "Species", "Tandem"))
    
    d_tandem <- d_tandem[d_tandem$Species != "convulsionarius",]
    d_tandem <- d_tandem[d_tandem$Genus != "Odontotermes" | (
                           d_tandem$Species != "distans" &
                           d_tandem$Species != "brunneus" &
                           d_tandem$Species != "assmuthi")   ,]
    d_tandem <- d_tandem[d_tandem$Genus != "Nasutitermes" | (
        d_tandem$Species != "nigriceps" &
        d_tandem$Species != "ephratae" &
        d_tandem$Species != "costalis")   ,]
    d_tandem <- d_tandem[d_tandem$Species != "dimorphus",]
    d_tandem <- d_tandem[d_tandem$Species != "unicolor",]
    d_tandem <- d_tandem[d_tandem$Species != "edentatus",]
    d_tandem <- d_tandem[d_tandem$Species != "suspensus",]
    d_tandem <- d_tandem[d_tandem$Species != "bequarerti",]
    d_tandem <- d_tandem[d_tandem$Species != "ochraceus",]
    d_tandem <- d_tandem[d_tandem$Species != "wheeleri",]
    d_tandem <- d_tandem[d_tandem$Species != "atlanticus",]
  }

  # matching data and tree
  {
    label_num_list <- numeric(length(d_tandem[,1]))
    options(warn=2)
    for(i_d in 1:length(d_tandem[,1])){
      
      genus <- d_tandem[i_d,]$Genus
      species <- d_tandem[i_d,]$Species
      
      label_num <- which(str_detect(labels, genus))
      if(length(label_num) > 1){
        label_num <- which(str_detect(labels, genus) & str_detect(labels, species))
        if(length(label_num) > 1){
          print(paste(genus, species, "multiple tips exist.", "use first one."))
          label_num <- label_num[1]
        }
      } else if(length(label_num) == 0){
        label_num <- NA
      }
      
      label_num_list[i_d] <- label_num
    }
    
    d_tandem$label <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
    d_tandem$label_num <- label_num_list
    d_tandem <- d_tandem[!is.na(d_tandem$label_num),]
    
    labels[d_tandem$label_num] <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
    tree$tip.label <- labels
    
    tree_dropped <- drop.tip(tree, labels[!(1:length(labels) %in% d_tandem$label_num)])
    tandem_tree <- ladderize(tree_dropped, right = TRUE)
    
  }
  
  # get vector for ace
  {
    leader <- as.factor(d_tandem$Leader2)
    incipient <- as.factor(d_tandem$Incipient)
    mature <- as.factor(d_tandem$Mature)
    
    names(leader) = names(incipient) = names(mature) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  }
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# dependence between female and male leaders
# ---------------------------------------------------------------------------- #
{
  female_leader <- leader == "female" | leader == "both"
  male_leader   <- leader == "male"   | leader == "both"
  names(female_leader) = names(male_leader) <- names(leader)
  res <- fitPagel(tandem_tree, female_leader, male_leader)
  res
}

# ---------------------------------------------------------------------------- #
# plot / analysis
# ---------------------------------------------------------------------------- #
{
  # Ancestral state reconstruction of tandem (four state model)
  save_plot = T
  {
    fitARD <- fitMk(tandem_tree, leader, model="ARD", pi="fitzjohn")
    fitER  <- fitMk(tandem_tree, leader, model="ER", pi="fitzjohn")
    fitSYM <- fitMk(tandem_tree, leader, model="SYM", pi="fitzjohn")
  
    fithrmER1 <- fitHRM(tandem_tree, leader, model="ER", pi="fitzjohn",
                       ncat=c(1,1,2,1), parallel=TRUE, ncores=20, niter=20)
    
    fithrmER2 <- fitHRM(tandem_tree, leader, model="ER", pi="fitzjohn",
                        ncat=c(2,2,2,2), parallel=TRUE, ncores=20, niter=20)
    
    anova(fitARD, fitER, fitSYM, fithrmER1, fithrmER2)
    plot(fithrmER1)
    
    ace1 <- ancr(fithrmER1, tips=TRUE)
    
    # plot
    leader_col   <- c("gray50", "#9370DB", "#D55E00", "#F63C00", "#009E73")
    plot(ace1, args.nodelabels=list(piecol=leader_col))
    if(save_plot){
      pdf("output/leader_ase.pdf", width=6, height=8)
      plot(ace1, args.nodelabels=list(piecol=leader_col))
      dev.off()
    }
    
    # ancestral state data
    save(ace1, file = "data_fmt/ace_tandem.rda")
    load("data_fmt/ace_tandem.rda")
    plot(tandem_tree)
    nodelabels() 
    round(ace1$ace, 2)
    
    # Kalotermitidae
    #0 both female female* male
    #0.32 0.28   0.29    0.00 0.11
    
    # Glyptotermes
    #0.38 0.42   0.11    0.00 0.08
    
    # G. nakajimai
    #0.93 0.04   0.02    0.00 0.02
    
  }  
  
  # Ancestral state reconstruction of mate system
  {
    # incipient
    {
      valid_tips <- names(incipient[!is.na(incipient)])
      pruned_tree <- drop.tip(tandem_tree, setdiff(tandem_tree$tip.label, valid_tips))
      incipient_pruned <- incipient[valid_tips]
      
      fitARD <- fitMk(pruned_tree, incipient_pruned, model="ARD", pi="fitzjohn")
      fitER  <- fitMk(pruned_tree, incipient_pruned, model="ER" , pi="fitzjohn")
      fitSYM <- fitMk(pruned_tree, incipient_pruned, model="SYM", pi="fitzjohn")
      
      anova(fitARD, fitER, fitSYM)
      
      ace1 <- ancr(fitER, tips=TRUE)
      
      plot(ace1, args.nodelabels=list(piecol=viridis(3)))
      
      # get value
      plot(pruned_tree)
      nodelabels() 
      round(ace1$ace, 2)
      
      # Kalotermitidae
      # 1 0 0
      
      # Glyptotermes
      # 1 0 0 
    }

    # mature
    {
      valid_tips <- names(mature[!is.na(mature)])
      pruned_tree <- drop.tip(tandem_tree, setdiff(tandem_tree$tip.label, valid_tips))
      mature_pruned <- mature[valid_tips]
      
      fitARD <- fitMk(pruned_tree, mature_pruned, model="ARD", pi="fitzjohn")
      fitER  <- fitMk(pruned_tree, mature_pruned, model="ER" , pi="fitzjohn")
      fitSYM <- fitMk(pruned_tree, mature_pruned, model="SYM", pi="fitzjohn")
      
      anova(fitARD, fitER, fitSYM)
      
      ace1 <- ancr(fitER, tips=TRUE)
      
      plot(ace1, args.nodelabels=list(piecol=viridis(3)))
      
      # get value
      plot(pruned_tree)
      nodelabels() 
      round(ace1$ace, 2)
      
      # Kalotermitidae
      # 0.99 0 0.01
      
      # Glyptotermes
      # 0.03 0 0.97
      
      # G. nakajimai
      # 0 0 1
    }
    
    
  }  
  
  # plot tip labels
  save_plot <- TRUE
  {
    partheno <- as.factor(d_tandem$Parthenogenesis)
    names(partheno) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
    
    is_tip <- tandem_tree$edge[,2] <= length(tandem_tree$tip.label)
    ordered_tips <- tandem_tree$edge[is_tip, 2]
    ordered_tips_label <- tandem_tree$tip.label[ordered_tips]
    ordered_leader    <- as.numeric(rev(leader[ordered_tips_label]))
    ordered_incipient <- as.numeric(rev(incipient[ordered_tips_label]))
    ordered_mature    <- as.numeric(rev(mature[ordered_tips_label]))
    ordered_partheno  <- as.numeric(rev(partheno[ordered_tips_label]))
    
    cols = viridis(4, end=0.5, alpha=0.5)
    col.show <- alpha(c(cols, "white"), 0.5)
    
    leader_col   <- c("gray50", "#9370DB", "#D55E00", "#009E73")
    mate_col     <- viridis(3, direction = -1)
    partheno_col <- c("#D72638", "#F49AC2", "#FFFFFF")  
    data.tile <- c(matrix(c(leader_col[ordered_leader],  
                            mate_col[ordered_incipient], 
                            mate_col[ordered_mature],
                            partheno_col[ordered_partheno]),
                          4, byrow=T))
    
    show_col(data.tile, ncol = 4, labels =F)
    
    if(save_plot){
      pdf("output/tip_lebel.pdf", width=5, height=6)
      show_col(data.tile, ncol = 4, labels =F)
      dev.off()
    }
    
  }
  
  # plot tree
  save_plot <- TRUE
  {
    plotTree(tandem_tree, fsize=0.5, ftype="i", offset=2 )
    if(save_plot){
      pdf("output/tree.pdf", width=5, height=6)
      plotTree(tandem_tree, fsize=0.5, ftype="i", offset=2 )
      dev.off()
    }
  }
}
  
