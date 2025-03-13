# Pylogenetic analysis

{
  require(phytools)
  require(stringr)
  library(viridis)
}

# ---------------------------------------------------------------------------- #
# prep data
# ---------------------------------------------------------------------------- #
{
  # tree data
  {
    #tree.file <- "data_raw/tree/run1_burn20_mcc_median.tree"
    tree.file <- "data_raw/tree/run5_400M_burnin20_mcc_median.tree"
    
    tree <- read.nexus(tree.file)
    labels <- tree$tip.label
    labels <- str_replace(labels, "Cubi_tenu", "Cubitermes_tenuiceps")
    labels <- str_replace(labels, "Tern_pall", "Ternicubitermes_pall")
  }

  # tandem data (remove some taxa without phylogeny information)
  {
    d_tandem <- read.csv("data_raw/tree/tandem_info_Mizumoto-etal-2022-PNAS.csv")
    
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
    d_tandem <- d_tandem[d_tandem$Species != "macrocephalus",]
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
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# plot / analysis
# ---------------------------------------------------------------------------- #
{
  # 
  {
  tandem <- d_tandem$Tandem
  names(tandem) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  
  
  fitER <- fitMk(tandem_tree, tandem, model="ER", pi="fitzjohn")
  fitARD <- fitMk(tandem_tree, tandem, model="ARD", pi="fitzjohn")

  fithrmER1<-fitHRM(tandem_tree,tandem,model="ER", parallel=T, ncat=c(2,2), pi="fitzjohn")
  fithrmARD1<-fitHRM(tandem_tree,tandem,model="ARD", parallel=T, ncat=c(2,2), pi="fitzjohn")
  
  fithrmER2<-fitHRM(tandem_tree,tandem,model="ER", parallel=T, ncat=c(1,2), pi="fitzjohn")
  fithrmARD2<-fitHRM(tandem_tree,tandem,model="ARD", parallel=T, ncat=c(1,2), pi="fitzjohn")
  
  fit_aov <- anova(fitER, fitARD, fithrmER1, fithrmARD1, fithrmER2, fithrmARD2)
  
  cols = viridis(3)
  anc_parity <- ancr(fithrmARD2,tips=TRUE)
  plot(anc_parity,args.nodelabels=list(piecol=cols),
       args.tiplabels=list(piecol=cols))
  
  plot(fithrmARD2)
  
  }
  
  
  
  {
    leader <- d_tandem$Leader
    names(leader) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
    leader <- as.factor(leader)
    
    X <- to.matrix(leader, levels(leader))
    labs <- row.names(X[apply(X, 1, sum) == 0,])
    X[apply(X, 1, sum) == 0, 2:4] <- 1
    
    fit1 <- fitMk(tandem_tree, X, model="ARD", pi="fitzjohn")
    fit2 <- fitMk(tandem_tree, leader, model="ER", pi="fitzjohn")
    fit3 <- fitMk(tandem_tree, X, model="SYM", pi="fitzjohn")
  
    anova(fit1, fit2, fit3)
    
    ncores<-min(c(10,parallel::detectCores()))
    
    fithrmARD2<-fitHRM(tandem_tree, leader, model="ARD", pi="fitzjohn",
                       ncat=c(1,1,2,1))
    
    anova(fit1, fit2, fit3, fit4)
    
    ace1 <- ancr(fit3,tips=TRUE)
    plot(ace1)
    tips<-sapply(labs,function(x,y) which(y==x),y=tandem_tree$tip.label)
    add.arrow(tandem_tree,tips,arrl=3,offset=2,lwd=2,
              col=palette()[4])
  }  
  
  
  {
    options(warn=0)
    leader <- d_tandem$Leader
    names(leader) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
    leader <- as.factor(leader)
    
    tandem_tree2 <- drop.tip(tandem_tree,  d_tandem$label[ is.na(d_tandem$Leader) ])
    leader <- d_tandem$Leader[ !is.na(d_tandem$Leader) ]
    names(leader) <- d_tandem$label[ !is.na(d_tandem$Leader) ]
    leader <- as.factor(leader)
    
    fit1 <- fitMk(tandem_tree2, leader, model="ARD", pi="fitzjohn")
    fit2 <- fitMk(tandem_tree2, leader, model="ER", pi="fitzjohn")
    fit3 <- fitMk(tandem_tree2, leader, model="SYM", pi="fitzjohn")
    
    fithrmER1<-fitHRM(tandem_tree2, leader, model="ER", pi="fitzjohn",
                      ncat=c(1,1,2,1), parallel=TRUE, ncores=20, niter=20)
    fithrmER2<-fitHRM(tandem_tree2, leader, model="ER", pi="fitzjohn",
                       ncat=c(2,2,2,2), parallel=TRUE, ncores=20, niter=20)
    
    anova(fit1, fit2, fit3, fithrmER1, fithrmER2)
    
    plot(fithrmER1)
    
    ace1 <- ancr(fithrmER1,tips=TRUE)
    plot(ace1)
  }

  {
    options(warn=0)
    leader <- d_tandem$Leader2
    names(leader) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
    leader <- as.factor(leader)
    
    fit1 <- fitMk(tandem_tree, leader, model="ARD", pi="fitzjohn")
    fit2 <- fitMk(tandem_tree, leader, model="ER", pi="fitzjohn")
    fit3 <- fitMk(tandem_tree, leader, model="SYM", pi="fitzjohn")
    
    fithrmER1<-fitHRM(tandem_tree, leader, model="ER", pi="fitzjohn",
                      ncat=c(1,1,2,1), parallel=TRUE, ncores=20, niter=20)
    fithrmER2<-fitHRM(tandem_tree, leader, model="ER", pi="fitzjohn",
                      ncat=c(2,2,2,2), parallel=TRUE, ncores=20, niter=20)
    
    fithrmSYM1<-fitHRM(tandem_tree, leader, model="SYM", pi="fitzjohn",
                      ncat=c(1,1,2,1), parallel=TRUE, ncores=20, niter=20)
    fithrmSYM2<-fitHRM(tandem_tree, leader, model="SYM", pi="fitzjohn",
                      ncat=c(2,2,2,2), parallel=TRUE, ncores=20, niter=20)
    
    
    anova(fit1, fit2, fit3, fithrmER1, fithrmER2, fithrmSYM1, fithrmSYM2)
    
    plot(fithrmER1)
    
    ace1 <- ancr(fithrmER1,tips=TRUE)
    plot(ace1)
  }
  
  
  {
    
    
    f_leader <- (leader == "female" | leader == "both")
    names(f_leader) <- d_tandem$label[ !is.na(d_tandem$Leader) ]
    m_leader <- (leader == "male" | leader == "both")
    names(m_leader) <- d_tandem$label[ !is.na(d_tandem$Leader) ]
    
    fit_fm <- fitPagel(tandem_tree2, f_leader, m_leader)
    plot(fit_cor, lwd.by.rate=TRUE)
    fit_f <- fitPagel(tandem_tree2, f_leader, m_leader, dep.var="x")
    fit_m <- fitPagel(tandem_tree2, f_leader, m_leader, dep.var="y")
    
    aic<-setNames(c(fit_fm$independent.AIC,
                    fit_m$dependent.AIC,
                    fit_f$dependent.AIC,
                    fit_fm$dependent.AIC),
                  c("independent","dependent m",
                    "dependent f","dependent f&m"))
    aic.w(aic)
    
    
    fitER <- fitMk(tandem_tree2, f_leader, model="ER", pi="fitzjohn")
    fitARD <- fitMk(tandem_tree2, f_leader, model="ARD", pi="fitzjohn")
    
    fithrmER2<-fitHRM(tandem_tree2, f_leader,model="ER", parallel=T, ncat=c(1,2))
    fithrmARD2<-fitHRM(tandem_tree2, f_leader,model="ARD", parallel=T, ncat=c(1,2))
    
    fit_aov <- anova(fitER, fitARD, fithrmER2, fithrmARD2)
    
    cols = viridis(3)
    anc_parity <- ancr(fitARD,tips=TRUE)
    plot(anc_parity,args.nodelabels=list(piecol=cols),
         args.tiplabels=list(piecol=cols))
    
    plot(fithrmARD2)
    
    
    
    fitER <- fitMk(tandem_tree2, m_leader, model="ER", pi="fitzjohn")
    fitARD <- fitMk(tandem_tree2, m_leader, model="ARD", pi="fitzjohn")
    
    fithrmER2<-fitHRM(tandem_tree2, m_leader,model="ER", parallel=T, ncat=c(1,2))
    fithrmARD2<-fitHRM(tandem_tree2, m_leader,model="ARD", parallel=T, ncat=c(1,2))
    
    fit_aov <- anova(fitER, fitARD, fithrmER2, fithrmARD2)
    
    cols = viridis(3)
    anc_parity <- ancr(fitARD,tips=TRUE)
    plot(anc_parity,args.nodelabels=list(piecol=cols),
         args.tiplabels=list(piecol=cols))
    
    plot(fithrmARD2)
    
    
  }
  
  
  
  
  
  #X["Glyptotermes_satsumensis",] <- c(0, 1-0.513776, 0.513776)
  #X["Glyptotermes_fuscus",] <- c(0, 1-0.2784785, 0.2784785)
  
  X[X[,1] == 1, ] <- 1
  X <- X[,-1]
  X[apply(X, 1, sum) == 0, ] <- 1
  
  labs <- rownames(X)[apply(X, 1, sum) == 3]
  
  
  fit1 <- fitMk(tandem_tree, X, model="ARD", pi="fitzjohn",
              lik.func="pruning", logscale=TRUE)
  ace1 <- ancr(fit1,tips=TRUE)
  plot(ace1,args.plotTree=list(direction="upwards"))
  tips<-sapply(labs,function(x,y) which(y==x),y=tandem_tree$tip.label)
  add.arrow(tandem_tree,tips,arrl=3,offset=2,lwd=2,
            col=palette()[4])
  
  ordered_model<-matrix(c(
    0,1,2,
    3,0,4,
    5,6,0),3,3,
    byrow=TRUE,
    dimnames=list(0:2,0:2))
  ordered_model
  fitER <- fitMk(tandem_tree, tandem, model="ER", pi="fitzjohn")
  fitARD <- fitMk(tandem_tree, tandem, model="ARD", pi="fitzjohn")
  
  fit_aov <- anova(fitER, fitARD)
  tandem_ancr<-ancr(fit_aov)
  library(viridis)
  cols = viridis(2)[1:2]
  plot(tandem_ancr,args.nodelabels=list(piecol=cols),
       args.tiplabels=list(piecol=cols),direction="upwards")
  
  fitER <- fitgammaMk(tandem_tree, tandem, rand_start=TRUE, model="ER")
  fitARD <- fitgammaMk(tandem_tree, tandem, rand_start=TRUE, model="ARD")
  
  fithrmER<-fitHRM(tandem_tree,tandem,model="ER", parallel=T)
  fithrmARD<-fitHRM(tandem_tree,tandem,model="ARD", parallel=T)
  plot(fithrmARD, spacer=0.4, asp=0.5, offset=0.05, signif=4)
  
  cols = viridis(4)
  anc_parity <- ancr(fithrmER,tips=TRUE)
  plot(anc_parity,args.nodelabels=list(piecol=cols),
       args.tiplabels=list(piecol=cols),direction="upwards")
  
  anc_parity <- ancr(fithrmARD, tips=TRUE)
  plot(anc_parity,args.nodelabels=list(piecol=cols),
       args.tiplabels=list(piecol=cols))
  
  
  
  
  
  
  
  cols = c(0,1)
  labels <- tandem_tree$tip.label
  tandem <- tandem[labels]
  plotTree(tandem_tree, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(tandem,sort(unique(tandem))), piecol = cols, cex = 0.3)
  
  fitER <- ace(tandem, tandem_tree, model="ER", type="discrete")
  fitARD <- ace(tandem, tandem_tree, model="ARD", type="discrete")
  
  AIC(fitER, k = 1)
  AIC(fitARD, k = 2)
  
  nodelabels(node = 1:tandem_tree$Nnode+Ntip(tandem_tree),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)
  nodelabels(node = 1:tandem_tree$Nnode+Ntip(tandem_tree),
             pie = fitARD$lik.anc, piecol = cols, cex=0.5)
  
  
  library(corHMM)
  rate.mat <- getRateCatMat(ntraits = 1, nstates = 2, nratecats = 2)
  
  # Fit the model
  fit <- corHMM(tandem_tree, tandem, rate.cat = 2, node.states = "marginal")
  
  # Print results
  print(fit)
  
  
  
  
  tandem <- d_tandem$Leader
  names(tandem) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  tandem[is.na(tandem)] <- "female"
  
  cols = 1:3
  labels <- tree2$tip.label
  tandem <- tandem[labels]
  plotTree(tree2, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(tandem, sort(unique(tandem))), piecol = cols, cex = 0.3)
  
  fitER <- ace(tandem, tree2, model="ER", type="discrete")
  nodelabels(node = 1:tree2$Nnode+Ntip(tree2),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)
  
  
  
  
  tandem <- d_tandem$Female
  names(tandem) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  tandem[is.na(tandem)] <- 1
  
  cols = 1:3
  labels <- tree2$tip.label
  tandem <- tandem[labels]
  plotTree(tree2, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(tandem, sort(unique(tandem))), piecol = cols, cex = 0.3)
  
  fitER <- ace(tandem, tree2, model="ER", type="discrete")
  nodelabels(node = 1:tree2$Nnode+Ntip(tree2),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)
  
  
  
  tandem <- d_tandem$Male
  names(tandem) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  tandem[is.na(tandem)] <- 1
  
  cols = 1:3
  labels <- tree2$tip.label
  tandem <- tandem[labels]
  plotTree(tree2, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(tandem, sort(unique(tandem))), piecol = cols, cex = 0.3)
  
  fitER <- ace(tandem, tree2, model="ER", type="discrete")
  nodelabels(node = 1:tree2$Nnode+Ntip(tree2),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)
  
  
  
  Female <- d_tandem$Female.active
  names(Female) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  tree3 <- drop.tip(tree2, names(Female)[is.na(Female)])
  Female <- Female[!is.na(Female)]
  Male <- d_tandem$Male.active
  names(Male) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  Male <- Male[!is.na(Male)]
  
  #fit <- fitPagel(tree3, Male, Female)
  
  plot(fit)
  
  cols = 1:2
  
  labels <- tree3$tip.label
  Female <- Female[labels]
  Male <- Male[labels]
  plotTree(tree3, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(Female, sort(unique(Female))), piecol = cols, cex = 0.3)
  fitER <- ace(Female, tree3, model="ER", type="discrete")
  nodelabels(node = 1:tree3$Nnode+Ntip(tree3),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)
  
  plotTree(tree3, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(Male, sort(unique(Male))), piecol = cols, cex = 0.3)
  fitER <- ace(Male, tree3, model="ER", type="discrete")
  nodelabels(node = 1:tree3$Nnode+Ntip(tree3),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)
  
  
  
  
  
  
  tandem <- d_tandem$Leader
  names(tandem) <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  tandem[is.na(tandem)] <- "female"
  
  dtemp <- subset(d_tandem, Tandem == 1)
  
  d_tandem$
  
  
  Female <- d_tandem$Female.active
  Male <- d_tandem$Male.active
  
  Female[d_tandem$Tandem == 1]
  Male[d_tandem$Tandem == 1]
  
  cols = 1:3
  labels <- tree2$tip.label
  tandem <- tandem[labels]
  plotTree(tree2, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(tandem, sort(unique(tandem))), piecol = cols, cex = 0.3)
  
  fitER <- ace(tandem, tree2, model="ER", type="discrete")
  nodelabels(node = 1:tree2$Nnode+Ntip(tree2),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)
  
  
  
  
  
  
  
  
  Female <- d_tandem$Female.active
  otus <- paste(d_tandem$Genus, d_tandem$Species, sep="_")
  names(Female) <- otus
  tree3 <- drop.tip(tree2, names(Female)[is.na(Female) | d_tandem$Tandem == 0])
  Female <- Female[!is.na(Female) & d_tandem$Tandem == 1]
  Male <- d_tandem$Male.active
  names(Male) <- otus
  Male <- Male[!is.na(Male) & d_tandem$Tandem == 1]
  
  
  cols = 1:2
  
  labels <- tree3$tip.label
  Female <- Female[labels]
  Male <- Male[labels]
  plotTree(tree3, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(Female, sort(unique(Female))), piecol = cols, cex = 0.3)
  fitER <- ace(Female, tree3, model="ER", type="discrete")
  nodelabels(node = 1:tree3$Nnode+Ntip(tree3),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)
  
  plotTree(tree3, fsize=0.8, ftype="i")
  tiplabels(pie = to.matrix(Male, sort(unique(Male))), piecol = cols, cex = 0.3)
  fitER <- ace(Male, tree3, model="ER", type="discrete")
  nodelabels(node = 1:tree3$Nnode+Ntip(tree3),
             pie = fitER$lik.anc, piecol = cols, cex=0.5)

}



tandem_tree2 <- drop.tip(tandem_tree,  d_tandem$label[ is.na(d_tandem$Leader) ])
leader <- d_tandem$Leader[ !is.na(d_tandem$Leader) ]

f_leader <- (leader == "female" | leader == "both")
names(f_leader) <- d_tandem$label[ !is.na(d_tandem$Leader) ]

fitER <- fitMk(tandem_tree2, f_leader, model="ER", pi="fitzjohn")
fitARD <- fitMk(tandem_tree2, f_leader, model="ARD", pi="fitzjohn")

fithrmER2<-fitHRM(tandem_tree2, f_leader,model="ER", parallel=T, ncat=c(1,2))
fithrmARD2<-fitHRM(tandem_tree2, f_leader,model="ARD", parallel=T, ncat=c(1,2))

fit_aov <- anova(fitER, fitARD, fithrmER2, fithrmARD2)

cols = viridis(3)
anc_parity <- ancr(fitARD,tips=TRUE)
plot(anc_parity,args.nodelabels=list(piecol=cols),
     args.tiplabels=list(piecol=cols))

plot(fithrmARD2)


m_leader <- (leader == "male" | leader == "both")
names(m_leader) <- d_tandem$label[ !is.na(d_tandem$Leader) ]

fitER <- fitMk(tandem_tree2, m_leader, model="ER", pi="fitzjohn")
fitARD <- fitMk(tandem_tree2, m_leader, model="ARD", pi="fitzjohn")

fithrmER2<-fitHRM(tandem_tree2, m_leader,model="ER", parallel=T, ncat=c(1,2))
fithrmARD2<-fitHRM(tandem_tree2, m_leader,model="ARD", parallel=T, ncat=c(1,2))

fit_aov <- anova(fitER, fitARD, fithrmER2, fithrmARD2)

cols = viridis(3)
anc_parity <- ancr(fitARD,tips=TRUE)
plot(anc_parity,args.nodelabels=list(piecol=cols),
     args.tiplabels=list(piecol=cols))

plot(fithrmARD2)

