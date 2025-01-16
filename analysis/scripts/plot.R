# Glyptotermes tandem
# N Mizumoto
{
  rm(list = ls())
  
  library(arrow)
  library(stringr)
  
  library(dplyr)
  library(tidyr)
  library(forcats)
  
  library(ggplot2)
  library(viridis)
  library(ggridges)
  
  library(circular)
  
  library(survival)
  library(survminer)
  library(car)
  library(coxme)
  
  species_list <- c("Gly_fus", "Gly_sat", "Gly_nak_sexual", "Gly_nak_asexual")
}

# ---------------------------------------------------------------------------- #
# data prep
# ---------------------------------------------------------------------------- #
{
  get_body_length <- function(df, species){
    df <- subset(df, df$frame %% 30 == 0)
    
    tbl0 <- sqrt((df$fHead_x - df$fCenter_x)^2) + sqrt((df$fCenter_x - df$fTip_x)^2)
    tbl1 <- sqrt((df$mHead_x - df$mCenter_x)^2) + sqrt((df$mCenter_x - df$mTip_x)^2)
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
    }
    # angle male to female
    {
      dir_vec_x = df$mHead_x - df$mTip_x
      dir_vec_y = df$mHead_y - df$mTip_y
      rel_vec_x = df$fCenter_x -  df$mCenter_x
      rel_vec_y = df$fCenter_y -  df$mCenter_y
      ang_mtof <- atan2(rel_vec_y, rel_vec_x) - atan2(dir_vec_y, dir_vec_x) + pi/2
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
    data_name <- paste0("data_fmt/", i_s, "_df.feather")
    df <- arrow::read_feather(data_name)
    # the pair male did not move for 1 h and female was next to it.
    df <- subset(df, video != "Gly_nak_NM2344_9") 
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
  ## 2D density plot
  {
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
        ) %>%
        select(-tbl)  
    }
    
    # plot
    {
      plot.range <- 2.5
      df_relative_normalized$species <- factor(df_relative_normalized$species, 
                                               levels = c("Gly_fus", "Gly_sat", 
                                                          "Gly_nak_sexual", "Gly_nak_asexual"))
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
    #df_temp <- df_temp[seq(1, nrow(df_angle), by = 10), ]
    
    angle_gly_fus <- circular(subset(df_temp, species =="Gly_fus")$angle)
    angle_gly_sat <- circular(subset(df_temp, species =="Gly_sat")$angle)
    angle_gly_nak_sexual <- circular(subset(df_temp, species =="Gly_nak_sexual")$angle)
    angle_gly_nak_asexual <- circular(subset(df_temp, species =="Gly_nak_asexual")$angle)
    
    watson.test(angle_gly_fus)
    watson.test(angle_gly_sat)
    watson.test(angle_gly_nak_sexual)
    watson.test(angle_gly_nak_asexual)
    
    # HR test is favorable but cannot complerete calculation in time. try with better PC
    # testdata = circular::rvonmises(20, mu = circular::circular(pi), kappa = 3)
    # HR_test(testdata, iter = 999)
    
    watson.two.test(angle_gly_nak_sexual, angle_gly_nak_asexual)
    watson.two.test(angle_gly_fus, angle_gly_sat)
    watson.two.test(angle_gly_sat, angle_gly_nak_sexual)
    watson.two.test(angle_gly_sat, angle_gly_nak_asexual)
    watson.two.test(angle_gly_sat, angle_gly_nak_sexual)
    
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
  }
  
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# tandem movement
# ---------------------------------------------------------------------------- #
{
  # data prep
  {
    load("data_fmt/df_bodylength.rda")
    load("data_fmt/df_tandem_analysis.rda")
    df_tandem_analysis <- df_tandem_analysis %>%
      left_join(df_bodylength[,c("video","tbl")], by = "video") %>%  
      mutate(
        dis = dis / tbl
      ) %>%
      select(-tbl) 
    
    # defining tandem running
    df_tandem_analysis <- df_tandem_analysis %>%
      mutate(
        status = case_when(
          dis > 2 ~ "sep",  
          dis <= 2 & ang_ftom < 0 & ang_mtof > 0 ~ "fleader", 
          dis <= 2 & ang_ftom > 0 & ang_mtof < 0 ~ "mleader", 
          dis <= 2 ~ "interact",
          TRUE ~ NA_character_ 
        )
      )
    
    df_tandem_analysis$status <- factor(df_tandem_analysis$status)
    df_tandem_analysis[df_tandem_analysis$species == "Gly_nak_asexual" &
                         df_tandem_analysis$status == "mleader",]$status <- "fleader"
    
    df_proportions <- df_tandem_analysis %>%
      group_by(species, video, status) %>%
      tally() %>% ungroup() %>%
      group_by(video) %>%
      mutate(proportion = n / sum(n)) %>% 
      ungroup()
    
    df_proportions$sep_proportion <- 
      rep(df_proportions[df_proportions$status == "sep",]$proportion, 
          rle(df_proportions$video)$length)
    
    df_proportions$status <- factor(df_proportions$status, levels = c("interact", "fleader", "mleader", "sep"))
    df_proportions$species <- str_replace(df_proportions$species, "_fus", "_1fus")
    df_proportions$species <- str_replace(df_proportions$species, "_sat", "_2sat")
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
  
  # data prep 2
  {
    df_wide <- df_proportions %>%
      select(species, video, status, proportion) %>%
      spread(key = status, value = proportion, fill = 0)
    
    df_wide$tandem <- df_wide$fleader + df_wide$mleader
    
    ggplot(df_wide, aes(x = fleader, y = mleader, col=species)) +
      geom_point(alpha = .75) +
      scale_color_viridis(discrete = T) + 
      coord_cartesian(xlim = c(0, .5), ylim = c(0, .5)) +
      theme_classic() +
      theme(aspect.ratio = 1) +
      labs(x = "Prop of female leader", y = "Prop of male leader")
    ggsave("output/tandem_proportion.pdf")
  }
  
  # tandem stability
  # data prep
  minimum_tandem_sec <- 2
  {
    load("data_fmt/df_bodylength.rda")
    load("data_fmt/df_tandem_analysis.rda")
    df_tandem_analysis <- df_tandem_analysis %>%
      left_join(df_bodylength[,c("video","tbl")], by = "video") %>%  
      mutate(
        dis = dis / tbl
      ) %>%
      select(-tbl) 
    
    df_tandem_analysis <- df_tandem_analysis %>%
      mutate(
        status = case_when(
          dis > 2 ~ "sep",  
          dis <= 2 & ang_ftom < 0 & ang_mtof > 0 ~ "fleader", 
          dis <= 2 & ang_ftom > 0 & ang_mtof < 0 ~ "mleader", 
          dis <= 2 ~ "interact",
          TRUE ~ NA_character_ 
        )
      )
    
    
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
    
    df <- df_tandem_analysis
    analyze_FPS = 5
    df <- df[df$frame %% (30/analyze_FPS) == 0,]
    
    videos <- unique(df$video)
    
    df_tandem = df_sep <- NULL
    thresh_sec <- 1
    for(i_v in 1:length(videos)){
      print(videos[i_v])
      df_temp <- subset(df, video == videos[i_v])
      tandem_status <- df_temp$status
      
      f_dis <- c(NA, sqrt( diff(df_temp$fx)^2 + diff(df_temp$fy)^2 ))
      m_dis <- c(NA, sqrt( diff(df_temp$mx)^2 + diff(df_temp$my)^2 ))
      
      f_leader = tandem_status == "fleader"
      f_leader <- !tandem.smoothing(!f_leader, analyze_FPS * minimum_tandem_sec)
      f_leader <- tandem.smoothing(f_leader, analyze_FPS * minimum_tandem_sec)
      
      m_leader = tandem_status == "mleader"
      m_leader <- !tandem.smoothing(!m_leader, analyze_FPS * minimum_tandem_sec)
      m_leader <- tandem.smoothing(m_leader, analyze_FPS * minimum_tandem_sec)
      
      if( sum(f_leader & m_leader) > 0 ){
        print("tandem duplication")
      }
      
      sep = (!f_leader & !m_leader & !(tandem_status == "interact"))
      sep <- !tandem.smoothing(!sep, analyze_FPS * minimum_tandem_sec)
      sep <- tandem.smoothing(sep, analyze_FPS * minimum_tandem_sec)
      
      rle_fleader <- rle(f_leader)
      rle_mleader <- rle(m_leader)
      rle_sep <- rle(sep)
      
      start <- 1
      for (i in seq_along(rle_fleader$lengths)) {
        run_length <- rle_fleader$lengths[i]
        
        if(rle_fleader$values[i]){
          end <- start + run_length - 1
          df_tandem <- rbind(df_tandem, data.frame(
            video = videos[i_v],
            leader = "female",
            duration = run_length,
            f_dis = sum(f_dis[start:end]),
            m_dis = sum(m_dis[start:end]),
            species = df_temp$species[1]
          ))
        }
        
        start <- start + run_length
      }
       
      start <- 1
      for (i in seq_along(rle_mleader$lengths)) {
        run_length <- rle_mleader$lengths[i]
        
        if(rle_mleader$values[i]){
          end <- start + run_length - 1
          df_tandem <- rbind(df_tandem, data.frame(
            video = videos[i_v],
            leader = "male",
            duration = run_length,
            f_dis = sum(f_dis[start:end]),
            m_dis = sum(m_dis[start:end]),
            species = df_temp$species[1]
          ))
        }
        
        start <- start + run_length
      }
      
      start <- 1 + rle_sep$lengths[1]
      for (i in seq_along(rle_sep$lengths)[-1]) {
        run_length <- rle_sep$lengths[i]
        
        if(rle_sep$values[i]){
          end <- start + run_length - 1
          df_sep <- rbind(df_sep, data.frame(
            video = videos[i_v],
            sep_event = paste0(videos[i_v], "_", i),
            tandem = tandem_status[start-1],
            time = (-(analyze_FPS*minimum_tandem_sec-1)):run_length,
            f_dis = (f_dis[(start-analyze_FPS*minimum_tandem_sec):end]),
            m_dis = (m_dis[(start-analyze_FPS*minimum_tandem_sec):end]),
            species = df_temp$species[1]
          ))
        }
        
        start <- start + run_length
      }
      
    }
    df_tandem[df_tandem$species == "Gly_nak_asexual" &
                df_tandem$leader == "male", "leader"] <- "female"
    
  }
  
  # female leader
  {
    r<-survfit(Surv(f_dis) ~ species, type = "kaplan-meier", 
               data = subset(df_tandem, leader == "female"))
    p <- ggsurvplot(fit = r, data = df_tandem,
               pval = F, pval.method = TRUE,
               risk.table = F, conf.int = FALSE,
               ncensor.plot = FALSE, size = 1,
               xlim = c(0,180) )
    p$plot + 
      scale_color_viridis(discrete = T) +
      scale_y_continuous(breaks = c(0,.5,1), labels = c("0", "0.5", "1")) +
      scale_x_continuous(breaks = c(0, 80, 160), labels = c("0", "80", "160")) +
      labs(x = "Leader moved ditance (mm)",
           y = "Tandem probability") + 
      theme_classic() +
      theme(legend.position = c(0.8,0.8),
            legend.title = element_blank())
    
    ggsave("output/female_leader_tandem_stability.pdf", width = 4, height= 3)
    m <- coxme(Surv(f_dis) ~ species + (1|video), 
               data = subset(df_tandem, leader == "female"))
    Anova(m)
  }  
  
  # male leader
  {
    r<-survfit(Surv(m_dis) ~ species, type = "kaplan-meier", 
               data = subset(df_tandem, leader == "male"))
    p <- ggsurvplot(fit = r, data = df_tandem,
               pval = F, pval.method = TRUE,
               risk.table = F, conf.int = FALSE,
               ncensor.plot = FALSE, size = 1,
               xlab="Time (min)", ggtheme = theme_bw()  + 
                 theme(aspect.ratio = 0.75),
               xlim=c(0,180))
    p$plot + 
      scale_color_viridis(discrete = T) +
      scale_y_continuous(breaks = c(0,.5,1), labels = c("0", "0.5", "1")) +
      scale_x_continuous(breaks = c(0, 80, 160), labels = c("0", "80", "160")) +
      labs(x = "Leader moved ditance (mm)",
           y = "Tandem probability") + 
      theme_classic() +
      theme(legend.position = c(0.8,0.8),
            legend.title = element_blank())
    ggsave("output/male_leader_tandem_stability.pdf", width = 4, height= 3)
    
    m <- coxme(Surv(m_dis) ~ species + (1|video), 
               data = subset(df_tandem, leader == "male"))
    Anova(m)
  }
}
# ---------------------------------------------------------------------------- #

