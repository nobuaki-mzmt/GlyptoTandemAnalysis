# Glyptotermes tandem
# plot all results + stats
# N Mizumoto

{
  rm(list = ls())
  options(warn = 0)
  
  library(arrow)
  library(stringr)
  library(data.table)
  
  library(dplyr)
  library(tidyr)
  library(forcats)
  
  library(ggplot2)
  library(viridis)
  library(ggridges)
  
  library(circular)
  library(NPCirc)
  library(CircMLE)
  
  library(compiler)
  
  library(survival)
  library(survminer)
  library(car)
  library(coxme)
  
  library(lme4)
  library(car)
  library(multcomp)
  
  species_list <- c("Gly_fus", "Gly_sat", "Gly_nak_sexual", "Gly_nak_asexual")
}

# ---------------------------------------------------------------------------- #
# data prep
# ---------------------------------------------------------------------------- #
{
  get_body_length <- function(df, species){
    df <- subset(df, df$frame %% 30 == 0)
    
    tbl0 <- sqrt((df$fHead_x - df$fCenter_x)^2) + 
      sqrt((df$fCenter_x - df$fTip_x)^2)
    tbl1 <- sqrt((df$mHead_x - df$mCenter_x)^2) + 
      sqrt((df$mCenter_x - df$mTip_x)^2)
    tbl0 = tapply(tbl0, df$video, mean) 
    tbl1 = tapply(tbl1, df$video, mean)
    
    df_temp <- data.frame(
      tbl = (tbl0+tbl1)/2,
      tbl0, tbl1,
      video = unique(df$video),
      species
    )
    
    return(df_temp)
  }
  
  relative_data <- function(df, relative2 = "female", species){
    df <- subset(df, df$frame %% 30 == 0)
    
    if(relative2 == "female"){
      dir_vec_x = df$fHead_x - df$fTip_x
      dir_vec_y = df$fHead_y - df$fTip_y
      rel_vec_x = df$mCenter_x - df$fCenter_x
      rel_vec_y = df$mCenter_y - df$fCenter_y
    } else if (relative2 == "male") {
      dir_vec_x = df$mHead_x - df$mTip_x
      dir_vec_y = df$mHead_y - df$mTip_y
      rel_vec_x = df$fCenter_x - df$mCenter_x
      rel_vec_y = df$fCenter_y - df$mCenter_y
    } else if (relative2 == "both") {
      dir_vec_x = c(df$mHead_x - df$mTip_x, df$fHead_x - df$fTip_x)
      dir_vec_y = c(df$mHead_y - df$mTip_y, df$fHead_y - df$fTip_y)
      rel_vec_x = c(df$mCenter_x - df$fCenter_x, df$fCenter_x - df$mCenter_x)
      rel_vec_y = c(df$mCenter_y - df$fCenter_y, df$fCenter_y - df$mCenter_y)
    }
    
    sta_ang <- atan2(rel_vec_y, rel_vec_x) - atan2(dir_vec_y, dir_vec_x) + pi/2
    dis <- sqrt(rel_vec_x^2 + rel_vec_y^2)
    
    df_stat_temp <- data.frame(
      angle = sta_ang,
      dis = dis,
      video = df$video,
      relative2 = relative2,
      species = species
    )
    
    return(df_stat_temp)
  }
  
  tandem_data <- function(df, species){
    # angle female to male
    {
      dir_vec_x = df$fHead_x - df$fTip_x
      dir_vec_y = df$fHead_y - df$fTip_y
      rel_vec_x = df$mCenter_x - df$fCenter_x
      rel_vec_y = df$mCenter_y - df$fCenter_y
      ang_ftom <- atan2(rel_vec_y, rel_vec_x) - atan2(dir_vec_y, dir_vec_x) + pi/2
      ang_ftom <- (ang_ftom + pi) %% (2 * pi) - pi
    }
    # angle male to female
    {
      dir_vec_x = df$mHead_x - df$mTip_x
      dir_vec_y = df$mHead_y - df$mTip_y
      rel_vec_x = df$fCenter_x -  df$mCenter_x
      rel_vec_y = df$fCenter_y -  df$mCenter_y
      ang_mtof <- atan2(rel_vec_y, rel_vec_x) - atan2(dir_vec_y, dir_vec_x) + pi/2
      ang_mtof <- (ang_mtof + pi) %% (2 * pi) - pi
    }
    dis <- sqrt(rel_vec_x^2 + rel_vec_y^2)
    
    df_temp <- data.frame(
      dis, ang_ftom, ang_mtof,
      frame = df$frame,
      video = df$video,
      mx = df$mCenter_x, fx = df$fCenter_x,
      my = df$mCenter_y, fy = df$fCenter_y,
      species = species
    )
    return(df_temp)
  }
  
  df_relative = df_bodylength = df_tandem_analysis <- NULL
  for(i_s in species_list){
    print(i_s)
    data_name <- paste0("data_fmt/", i_s, "_df.feather")
    df <- arrow::read_feather(data_name)
    # the pair male did not move for 1 h and female was next to it.
    df <- subset(df, video != "Gly_nak-sex_NM2344_9") 
    if(i_s != "Gly_nak_asexual"){
      df_relative <- rbind(df_relative,
                           relative_data(df, "female", i_s),
                           relative_data(df, "male", i_s)
      )
    } else {
      df_relative <- rbind(df_relative, relative_data(df, "both", i_s))
    }
    
    df_bodylength <- rbind(df_bodylength, get_body_length(df, i_s))
    df_tandem_analysis <- rbind(df_tandem_analysis, tandem_data(df, i_s))
    
  }
  
  save(df_bodylength, file = "data_fmt/df_bodylength.rda")
  save(df_relative, file = "data_fmt/df_relative.rda")
  save(df_tandem_analysis, file = "data_fmt/df_tandem_analysis.rda")
  
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# plot relative positioning
# ---------------------------------------------------------------------------- #
{
  load("data_fmt/df_relative.rda")
  load("data_fmt/df_bodylength.rda")
  # data prep
  {
    df_relative$x = cos(df_relative$angle)*df_relative$dis
    df_relative$y = sin(df_relative$angle)*df_relative$dis
    
    df_relative$relative2[df_relative$relative2 == "both"] <- "female"
    plot.range <- 15
    
    df_relative <- df_relative %>%
      group_by(species, relative2) %>%
      mutate(weight = 1 / n())  
    
    df_relative_normalized <- df_relative %>%
      left_join(df_bodylength[,c("video","tbl")], by = "video") %>% 
      mutate(
        x = x / tbl,  
        y = y / tbl,  
        dis = dis / tbl
      )
    df_relative_normalized$species <- factor(df_relative_normalized$species, 
                                             levels = c("Gly_sat", "Gly_fus",
                                                        "Gly_nak_sexual", "Gly_nak_asexual"))
  }
  
  ## 2D density plot
  {
    plot.range <- 2.5
    ggplot(subset(df_relative_normalized, abs(x) < plot.range & abs(y) < plot.range), 
           aes(x = x, y = y, weight = weight)) +
      geom_bin2d(bins = 25) +
      scale_fill_viridis(option = "viridis", 
                         guide = guide_colorbar(barwidth = 0.5, barheight = 3)) +
      scale_x_continuous(expand = c(0, 0.01),
                         limits = c(-plot.range, plot.range)) +
      scale_y_continuous(expand = c(0, 0.01), 
                         limits = c(-plot.range, plot.range), position = "right") + 
      theme_classic() +
      coord_fixed() + 
      xlab("Distance left−right (body length)")+
      ylab("Distance front−back (body length)")+
      theme(
        legend.position = c(0.8, 0.15),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text.y.right = element_text(size = 7), 
        axis.text.x.bottom = element_text(size = 7),
        strip.text.y.left = element_text(angle = 0, size = 7),
        strip.text.x = element_text(size = 7)
      ) + 
      facet_grid(species ~ relative2, switch = "y")
    ggsave("output/2d_hist.pdf", width = 4, height = 7)
  }

  ## Angles
  # plot
  {
    df_relative_normalized$angle <- (df_relative_normalized$angle + pi) %% (2 * pi) - pi
    ggplot(subset(df_relative_normalized, relative2 == "female" & dis < 2), 
           aes(x = angle, y = forcats::fct_rev(species))) +
      geom_density_ridges(stat = "binline", bins = 40, scale = 0.7) +
      labs(y = "", x = "Relaive angle from female (rad)") + 
      scale_x_continuous(
        breaks = c(-pi, -pi/2, 0, pi/2, pi),
        labels = c("-π", "-π/2", "0", "π/2", "π")
      ) + 
      geom_vline(xintercept = c(-pi/2, pi/2),  linetype = 2) +
      theme_classic()
    ggsave("output/angle_hist.pdf", width = 3, height = 3)
  }
  
  # stat
  {
    df_temp <- subset(df_relative_normalized, relative2 == "female" & dis < 2)
    
    min(table(df_temp$species))
    
    angle_gly_fus <- circular(subset(df_temp, species =="Gly_fus")$angle)
    angle_gly_sat <- circular(subset(df_temp, species =="Gly_sat")$angle)
    angle_gly_nak_sexual <- circular(subset(df_temp, species =="Gly_nak_sexual")$angle)
    angle_gly_nak_asexual <- circular(subset(df_temp, species =="Gly_nak_asexual")$angle)
    
    # HR tests T function 
    HermansRasson2Tunc <- function(sample){
      n <- length(sample)
      total <- 0
      for (i in 1:n){
        for (j in 1:n){ total <- total + abs(abs(sample[i]-sample[j])-pi)-(pi/2)
        total <- total - (2.895*(abs(sin(sample[i]-sample[j]))-(2/pi)))}}
      T <- total/n
      return(T)}
    HermansRasson2T <- cmpfun(HermansRasson2Tunc)
    
    subsample_HermansRasson2T <- function(data, n_min = 1000, iterations = 100, name = "start") {
      print(name)
      HermansRasson2T_stat <- numeric(iterations)
      for (i in 1:iterations) {
        print(paste(i, "/", iterations))
        subsampled_data <- data[sample(1:length(data), size = n_min, replace = FALSE)]
        HermansRasson2T_stat[i] <- HermansRasson2T(subsampled_data)
      }
      return(HermansRasson2T_stat)
    }
    
    if(F){
      # this part takes time for computation. load pre-runned results in default
      gly_fus_HRT <- subsample_HermansRasson2T(angle_gly_fus, name = "Gly fus")
      gly_sat_HRT <- subsample_HermansRasson2T(angle_gly_sat)
      gly_nak_s_HRT <- subsample_HermansRasson2T(angle_gly_nak_sexual)
      gly_nak_a_HRT <- subsample_HermansRasson2T(angle_gly_nak_asexual)
      HRT_data <- data.frame(
        Statistic = c(gly_fus_HRT, gly_sat_HRT, gly_nak_s_HRT, gly_nak_a_HRT),
        Dataset = rep(
          c("gly_fus", "gly_sat", "gly_nak_sexual", "gly_nak_asexual"),
          each = length(gly_nak_a_HRT)
        )
      )
    } else {
      load("data_fmt/HRT_data.rda")
    }
    
    ggplot(HRT_data, aes(x = Statistic, fill = Dataset)) +
      geom_histogram(alpha = 0.6, bins = 100, position = "identity") +
      theme_minimal() +
      scale_fill_viridis(discrete =T)+
      labs(x = "HRT Statistic",
           y = "Frequency") +
      theme(legend.position = "none") +
      theme_classic()
    ggsave("output/HRT_stat.pdf", width = 4, height = 3)
  }
  
  ## Distance
  # plot
  {
    ggplot(subset(df_relative_normalized, relative2 == "female" & dis < 10), 
           aes(x = dis, y = forcats::fct_rev(species))) +
      geom_density_ridges(stat = "binline", bins = 40, scale = 0.7) +
      labs(y = "", x = "Relaive angle from female (rad)") + 
      scale_x_continuous(
        breaks = c(0, 2, 10),
        labels = c("0", "2", "10")
      ) +
      geom_vline(xintercept = 2) + 
      theme_classic()
    
    ggsave("output/dis_hist.pdf", width = 3, height = 2)
    
    
    tbl_mean <- tapply(df_bodylength$tbl, df_bodylength$species, mean) 
    res <- NULL
    for(i_s in species_list){
      print(i_s)
      data_name <- paste0("data_fmt/", i_s, "_df.feather")
      df <- arrow::read_feather(data_name)
      # the pair male did not move for 1 h and female was next to it.
      df <- subset(df, video != "Gly_nak-sex_NM2344_9") 
      
      df <- df[df$frame %% 30 == 0,]
      videos <- unique(df$video)
      
      for(i in 1:1000){
        video_com <- sample(videos, 2, replace = F)
        df1 <- subset(df, video == video_com[1])
        df2 <- subset(df, video == video_com[2])
        
        df_row_lim <-min(dim(df1)[1], dim(df2)[1])
        df1 <- df1[1:df_row_lim,]
        df2 <- df2[1:df_row_lim,]
        
        dis <- sqrt( (df1$fCenter_x - df2$mCenter_x)^2 + (df1$fCenter_y - df2$mCenter_y)^2 )
        res <- c(res, mean(dis < tbl_mean[i_s]))
      }
      
    }
    df_res <- data.frame(res, species = rep(species_list, each = 1000))
    random_dis <- tapply(df_res$res, df_res$species, mean)
    
    library(multimode)
    
    mode_res <- modetest(sample(subset(df_relative_normalized, species == "Gly_fus")$dis, 1000))
    print(mode_res)
    
    mean_2body <- tapply(df_relative_normalized$dis<2, df_relative_normalized$video, mean)
    videos <- names(mean_2body)
    colony <- str_split_fixed(videos, "_", 4)[,3]
    species <- str_split_fixed(videos, "_", 4)[,2]
    df_dis <- data.frame(
      mean_2body, videos, colony, species
    )
    ggplot(df_dis, aes(x=species, y = mean_2body)) +geom_point()
    y = df_dis$mean_2body
    df_dis$logit_prop = log((y)/(1-y))
    df_dis$colony <- factor(df_dis$colony)
    df_dis$species[df_dis$species == "nak-asex"] <- "asex"
    df_dis$species[df_dis$species == "nak-sex"] <- "sex"
    r <- lmer(logit_prop ~ species + (1|colony), data = df_dis)
    Anova(r)
  }
  
  if(F){
    HermansRasson2T_stat <- numeric(length(videos))
    for(i_v in 1:length(videos)){
      df_temp2 <- subset(df_temp, video == videos[i_v])
      print(paste(i_v, "/", length(videos)))
      ang_data <- circular(df_temp2$angle)
      if(length(ang_data) > 100){
        for(i in 1:100){
          ang_data <- sample(ang_data, 100, replace = F)
          HermansRasson2T_stat[i_v] <- HermansRasson2T_stat[i_v] + HermansRasson2T(ang_data)
        }
        HermansRasson2T_stat[i_v] <- HermansRasson2T_stat[i_v] / 100
      } else {
        HermansRasson2T_stat[i_v] <- HermansRasson2T(ang_data)
      }
    }
    df_dis$HermansRasson2T_stat <- HermansRasson2T_stat
    ggplot(df_dis, aes(x=species, y = HermansRasson2T_stat)) +geom_point()
    df_dis$colony <- factor(df_dis$colony)
    df_dis$species[df_dis$species == "nak-asex"] <- "asex"
    df_dis$species[df_dis$species == "nak-sex"] <- "sex"
    r <- lmer(HermansRasson2T_stat ~ species + (1|colony), data = df_dis)
    Anova(r)
  }
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# tandem movement
# ---------------------------------------------------------------------------- #
{
  analyze_FPS = 5
  minimum_tandem_sec <- 5
  minimum_sep_sec <- 2
  
  tandem.smoothing <- function(vec, min.sec){ 
    if(sum(vec)>0){
      timing <- which(vec)[c(T, diff(which(vec))>1)]
      end    <- which(vec)[c(diff(which(vec))>1,T)]
      for(fi in 1:length(timing)){
        if(length( vec[timing[fi]:end[fi]]) < min.sec ){
          vec[timing[fi]:end[fi]] <- F
        }
      }
    }
    return(vec)
  }
  
  # tandem proportion
  {
    # data prep
    {
      load("data_fmt/df_bodylength.rda")
      load("data_fmt/df_tandem_analysis.rda")
      
      df_tandem_analysis <- df_tandem_analysis[df_tandem_analysis$frame %% (30/analyze_FPS) == 0,]
      
      df_tandem_analysis <- df_tandem_analysis %>%
        left_join(df_bodylength[,c("video","tbl")], by = "video") %>%  
        mutate(
          dis = dis / tbl
        ) 
      
      # defining tandem running
      interactions <- df_tandem_analysis$dis <= 2
      !tandem.smoothing(!interactions, analyze_FPS * minimum_tandem_sec)
      interactions <- tandem.smoothing(interactions, analyze_FPS * minimum_sep_sec)
      
      df_tandem_analysis <- df_tandem_analysis %>%
        mutate(
          status = case_when(
            !interactions ~ "sep",  
            interactions & ang_ftom < 0 & ang_mtof > 0 ~ "fleader", 
            interactions & ang_ftom > 0 & ang_mtof < 0 ~ "mleader", 
            interactions ~ "interact",
            TRUE ~ NA_character_ 
          )
        )
      
      df_tandem_analysis$status <- factor(df_tandem_analysis$status, 
                                          levels = c("interact", "fleader", "mleader", "sep"))
      df_tandem_analysis[df_tandem_analysis$species == "Gly_nak_asexual" &
                           df_tandem_analysis$status == "mleader",]$status <- "fleader"
      
      df_proportions <- df_tandem_analysis %>%
        group_by(species, video, status) %>%
        tally() %>% ungroup() %>%
        group_by(video) %>%
        mutate(proportion = n / sum(n)) %>% 
        mutate(total = sum(n)) %>% 
        ungroup()
      
      df_proportions$sep_proportion <- 
        rep(df_proportions[df_proportions$status == "sep",]$proportion, 
            rle(df_proportions$video)$length)
      
      df_proportions$status <- factor(df_proportions$status, 
                                      levels = c("interact", "fleader", "mleader", "sep"))
      df_proportions$species <- str_replace(df_proportions$species, "_fus", "_2fus")
      df_proportions$species <- str_replace(df_proportions$species, "_sat", "_1sat")
      df_proportions$species <- str_replace(df_proportions$species, "_nak_sexual", "_3nak")
      df_proportions$species <- str_replace(df_proportions$species, "_nak_asexual", "_4nak")
      df_proportions <- df_proportions %>%
        arrange(species, sep_proportion) %>%
        mutate(video = factor(video, levels = unique(video)))
    }
    
    # plot for each video
    {
      ggplot(df_proportions, aes(x = video, y = proportion, fill = status)) +
        geom_bar(stat = "identity", position = "stack", linewidth = .1, color = NA) +
        #scale_fill_viridis(discrete = TRUE, option = "G") +
        scale_fill_manual(values = c( "gray50", "red2", "blue4","grey95")) + 
        theme_classic() + 
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "top") +
        scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        labs(x = "", y = "Proportion of time")
      
      ggsave("output/tandem_prop_pair.pdf", width = 5, height = 2)
      
      
    }
    
    df_proportions$interact_prop <- 1 - df_proportions$sep_proportion
    
    df_proportions$colony <- str_split_fixed(df_proportions$video, "_", 4)[,3]
    
    df_subprop <- subset(df_proportions, status == "fleader" | status == "mleader")
    
    df_prop_stat <- data.frame(
      video = sort(unique(df_tandem_analysis$video)),
      sep = tapply(df_tandem_analysis$status == "sep", df_tandem_analysis$video, sum),
      fleader = tapply(df_tandem_analysis$status == "fleader", df_tandem_analysis$video, sum),
      mleader = tapply(df_tandem_analysis$status == "mleader", df_tandem_analysis$video, sum),
      interact = tapply(df_tandem_analysis$status == "interact", df_tandem_analysis$video, sum),
      frame = tapply(df_tandem_analysis$frame, df_tandem_analysis$video, max)/6+1,
      tandem_prop = tapply(df_tandem_analysis$status == "fleader" | 
                             df_tandem_analysis$status == "mleader", df_tandem_analysis$video, mean)>.4
    )
    row.names(df_prop_stat) <- NULL
    df_prop_stat$colony <- str_split_fixed(df_prop_stat$video, "_", 4)[,3]
    df_prop_stat$species <- str_split_fixed(df_prop_stat$video, "_", 4)[,2]
    
    r <- glmer(cbind(fleader+mleader, frame ) ~ 
                 species + (1|colony) + (1|colony:video), 
               family = "binomial", data = df_prop_stat)
    Anova(r)
    
    df_prop_fuscus <- df_proportions %>% 
      filter(species == "Gly_2fus") %>%
      filter(status == "fleader" | status == "mleader")
    r <- glmer(cbind(n, total) ~ status + (1|colony) + (1|colony:video), 
               family = "binomial", data = df_prop_fuscus)
    Anova(r)
    summary(r)$coefficient
    
    df_prop_sat <- df_proportions %>% 
      filter(species == "Gly_1sat") %>%
      filter(status == "fleader" | status == "mleader")
    
    r <- glmer(cbind(n, total) ~ status +(1|colony:video), 
               family = "binomial", data = df_prop_sat)
    Anova(r)
    summary(r)$coefficient
    
    
    df_prop_nak <- df_proportions %>% 
      filter(species == "Gly_3nak") %>%
      filter(status == "fleader" | status == "mleader")
    
    r <- glmer(cbind(n, total) ~ status +(1|colony:video), 
               family = "binomial", data = df_prop_nak)
    Anova(r)
    summary(r)$coefficient
    
    table(df_prop_stat[, c("tandem_prop", "species")])
  }
  
  # tandem stability
  {
    # data prep
    {
      load("data_fmt/df_bodylength.rda")
      load("data_fmt/df_tandem_analysis.rda")
      df <- df_tandem_analysis
      
      df <- df %>%
        left_join(df_bodylength[,c("video","tbl")], by = "video") %>%  
        mutate(
          dis = dis / tbl
        ) 
      
      df <- df[df$frame %% (30/analyze_FPS) == 0,]
      
      videos <- unique(df$video)
      
      df_tandem = df_sep <- NULL
      for(i_v in 1:length(videos)){
        print(videos[i_v])
        df_temp <- subset(df, video == videos[i_v])
        
        interaction <- df_temp$dis < 2
        interaction <- !tandem.smoothing(!interaction, analyze_FPS * minimum_tandem_sec)
        interaction <- tandem.smoothing(interaction, analyze_FPS * minimum_sep_sec)
        
        rle_interaction <- rle(interaction)
        
        f_dis <- c(NA, sqrt( diff(df_temp$fx)^2 + diff(df_temp$fy)^2 ))
        m_dis <- c(NA, sqrt( diff(df_temp$mx)^2 + diff(df_temp$my)^2 ))
        
        start <- 1
        for (i in seq_along(rle_interaction$lengths)) {
          run_length <- rle_interaction$lengths[i]
          end <- start + run_length - 1
          
          if(rle_interaction$values[i]){
            
            f_lead_prop <- mean(df_temp[start:end,]$ang_ftom < 0 & df_temp[start:end,]$ang_mtof > 0)
            m_lead_prop <- mean(df_temp[start:end,]$ang_ftom > 0 & df_temp[start:end,]$ang_mtof < 0)
            
            print(paste0("tandem id: ", i, "/", length(rle_interaction$values), 
                         "; start: ", start/analyze_FPS,
                         "; flead: ", round(f_lead_prop,2),
                         "; mlead: ", round(m_lead_prop,2)))
            if(f_lead_prop > .5){
              leader <- "female"
            } else if (m_lead_prop > .5){
              leader <- "male"
            } else {
              leader <- "none"
            }
          } else {
            leader <- "sep"
          }
          
          df_tandem <- rbind(df_tandem, data.frame(
            video = videos[i_v],
            leader,
            duration = run_length,
            status = sum(i != max(seq_along(rle_interaction$lengths))), # 0: censored
            f_dis = sum(f_dis[start:end], na.rm=T),
            m_dis = sum(m_dis[start:end], na.rm=T),
            species = df_temp$species[1]
          ))
          
          start <- start + run_length
        }
      }
      df_tandem[df_tandem$species == "Gly_nak_asexual" &
                  df_tandem$leader == "male", "leader"] <- "female"
      
    }
    
    #
    max(subset(df_tandem, species == "Gly_nak_sexual" & leader == "male")$m_dis)
    max(subset(df_tandem, species == "Gly_nak_sexual" & leader == "female")$f_dis)
    max(subset(df_tandem, species == "Gly_nak_asexual" & leader == "female")$f_dis)
    long_tandem <- rbind(df_tandem[df_tandem$f_dis > 1000 & df_tandem$leader == "female",],
      df_tandem[df_tandem$m_dis > 1000 & df_tandem$leader == "male",])
    unique(long_tandem$video)
    table((long_tandem$species))
    
    #
    {
      r<-survfit(Surv(f_dis, status) ~ species, type = "kaplan-meier", 
                 data = subset(df_tandem, leader == "female"))
      p <- ggsurvplot(fit = r, data = df_tandem,
                 pval = F, pval.method = TRUE,
                 risk.table = F, conf.int = FALSE,
                 ncensor.plot = T, size = 1,
                 xlim = c(1, 3000) )
      p$plot + 
        scale_color_viridis(discrete = T) +
        scale_y_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
        scale_x_continuous(breaks = c(0, 1500, 3000), labels = c("0", "1500", "3000")) +
        #scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1k", "10k")) +
        labs(x = "Leader moved ditance (mm)",
             y = "Tandem probability") + 
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_blank())
      
      ggsave("output/female_leader_tandem_stability.pdf", width = 4, height= 3)
      df_tandem$species <- factor(df_tandem$species)
      m <- coxme(Surv(f_dis, status) ~ species + (1|video), 
                 data = subset(df_tandem, leader == "female"))
      Anova(m)
      summary(glht(m, linfct = mcp(species = "Tukey")))
    }  
    
    # male leader
    {
      r<-survfit(Surv(m_dis, status) ~ species, type = "kaplan-meier", 
                 data = subset(df_tandem, leader == "male"))
      p <- ggsurvplot(fit = r, data = df_tandem,
                 pval = F, pval.method = TRUE,
                 risk.table = F, conf.int = FALSE,
                 ncensor.plot = FALSE, size = 1,
                 xlab="Time (min)", ggtheme = theme_bw()  + 
                   theme(aspect.ratio = 0.75),
                 xlim=c(0, 3000))
      p$plot + 
        scale_color_viridis(discrete = T) +
        scale_y_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
        scale_x_continuous(breaks = c(0, 1500, 3000), labels = c("0", "1500", "3000")) +
        #scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1k", "10k")) +
        labs(x = "Leader moved ditance (mm)",
             y = "Tandem probability") + 
        theme_classic() +
        theme(legend.position = c(0.8,0.8),
              legend.title = element_blank())
      ggsave("output/male_leader_tandem_stability.pdf", width = 4, height= 3)
      
      df_temp <- subset(df_tandem, leader == "male")
      df_temp$species <- factor(df_temp$species)
      m <- coxme(Surv(m_dis, status) ~ species + (1|video), 
                 data = df_temp)
      Anova(m)
      summary(glht(m, linfct = mcp(species = "Tukey")))
    }
    
    df_fuscus <- df_tandem %>% 
      filter(species == "Gly_fus") %>%
      filter(leader == "female" | leader == "male") %>%
      mutate(leader_dis = case_when(
        leader == "female" ~ f_dis,
        leader == "male" ~ m_dis,
        TRUE ~ NA_real_
      )) %>%
      dplyr::select(-species, -f_dis, -m_dis)
    
    m <- coxme(Surv(leader_dis, status) ~ leader + (1|video), 
               data = df_fuscus)
    Anova(m)
    
    df_fuscus <- df_tandem %>% 
      filter(species == "Gly_sat") %>%
      filter(leader == "female" | leader == "male") %>%
      mutate(leader_dis = case_when(
        leader == "female" ~ f_dis,
        leader == "male" ~ m_dis,
        TRUE ~ NA_real_
      )) %>%
      dplyr::select(-species, -f_dis, -m_dis)
    
    m <- coxme(Surv(leader_dis, status) ~ leader + (1|video), 
               data = df_fuscus)
    Anova(m)
    
    
  }
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# breeding structure
# ---------------------------------------------------------------------------- #
{
  df <- fread("data_raw/Glyptotermes_reproductive_number.csv")
  df$species <- factor(df$species, levels = c("satsumensis", "fuscus", "nakajimai_sex", "nakajimai_asex"))
  # show some stats
  
  # data count
  table(df[, c("species", "development")])
  # data lacking primary reproductives
  subset(df, p_total < 2)
  
  #
  df_stat <- subset(df, p_total > 1)
  
  df_stat$monogamous <- df_stat$p_total == 2
  
  df_summary <- df_stat %>%
    group_by(species, development, monogamous) %>%
    summarise(count = n()) %>%
    group_by(species, development) %>%
    mutate(prop = count / sum(count), total = sum(count))
  
  ggplot(df_summary, aes(x = interaction(development, species), y = prop, fill = factor(monogamous))) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.75, width = 0.5) +
    ylab("Proportion Monogamous") +
    theme_classic() +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
    scale_fill_viridis(discrete = T) +
    theme(legend.position = "none", aspect.ratio = 0.75) +
    xlab("")
  ggsave("output/mate_type_glypto.pdf", width = 3, height = 2)
  
  
  
  #ggplot(subset(df_stat, development == "mature"), 
  ggplot(df_stat, 
         aes(x = fct_rev(species), y = p_total,  col = development)) + 
    #geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
    geom_jitter(width = 0.075, height = 0, alpha = 0.75, size = 0.5)+
    scale_color_viridis(discrete = T, end = 0.5, direction = -1) +
    theme_classic() +
    coord_flip() +
    scale_y_log10(breaks = c(2, 5, 10, 20, 50))+
    labs(x = "", y = "") +
    theme(legend.position = "none", aspect.ratio = 1.8, axis.text.y = element_blank())
  ggsave("output/num_repro_glypto.pdf", width = 2.5, height = 3)
  
  
  # comp species 
  r <- glm(p_total ~ species, family = "poisson",
           data = subset(df_stat, development == "incipient"))
  Anova(r)
  multicomparison<-glht(r, linfct = mcp(species = "Tukey"))
  summary(multicomparison)
  
  r <- glm(p_total ~ species, family = "poisson",
           data = subset(df_stat, development == "mature"))
  Anova(r)
  multicomparison<-glht(r, linfct = mcp(species = "Tukey"))
  
  summary(multicomparison)
  
  # comp incipient vs mature for each species
  r <- glm(p_total ~ development, family = "poisson",
           data = subset(df_stat, species == "nakajimai_asex"))
  Anova(r)
  
  r <- glm(p_total ~ development, family = "poisson",
           data = subset(df_stat, species == "satsumensis"))
  Anova(r)
  
  r <- glm(p_total ~ development, family = "poisson",
           data = subset(df_stat, species == "fuscus"))
  Anova(r)
  
}
# ---------------------------------------------------------------------------- #

