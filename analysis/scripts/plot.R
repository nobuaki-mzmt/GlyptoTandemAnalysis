
{
  library(arrow)
  library(stringr)
  
  library(dplyr)
  library(forcats)
  
  library(ggplot2)
  library(viridis)
  library(ggridges)
  
  species_list <- c("Gly_sat", "Gly_fus", "Gly_nak_sexual", "Gly_nak_asexual")
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
      rel_vec_y = df$fCenter_y -  df$mCenter_x
      ang_mtof <- atan2(rel_vec_y, rel_vec_x) - atan2(dir_vec_y, dir_vec_x) + pi/2
    }
    dis <- sqrt(rel_vec_x^2 + rel_vec_y^2)
    
    df_temp <- data.frame(
      dis, ang_ftom, ang_mtof,
      frame = df$frame,
      video = df$video,
      species = species
    )
    
    return(df_temp)
  }
  
  df_relative = df_bodylength = df_tandem_analysis <- NULL
  for(i_s in species_list){
    data_name <- paste0("data_fmt/", i_s, "_df.feather")
    df <- arrow::read_feather(data_name)
    df <- subset(df, video != "Gly_nak_NM2344_9") # the pair male did not move for 1 h and female was next to it.
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
  save(df_relative, df_bodylength, df_tandem_analysis, file = "data_fmt/df_analysis.rda")
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# plot relative positioning
# ---------------------------------------------------------------------------- #
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
  
  plot.range <- 2.5
  df_relative_normalized$species <- factor(df_relative_normalized$species, 
                                           levels = c("Gly_fus", "Gly_sat", 
                                                      "Gly_nak_sexual", "Gly_nak_asexual"))
  ggplot(subset(df_relative_normalized, abs(x) < plot.range & abs(y) < plot.range), 
         aes(x = x, y = y, weight = weight)) +
    geom_bin2d(bins = 40) +
    scale_fill_viridis(option = "viridis", guide = guide_colorbar(barwidth = 0.5, barheight = 3)) +
    scale_x_continuous(expand = c(0, 0.01), limits = c(-plot.range, plot.range)) +
    scale_y_continuous(expand = c(0, 0.01), limits = c(-plot.range, plot.range), position = "right") + 
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
      axis.text.y.right = element_text(size = 7), # Adjust y-axis text size
      axis.text.x.bottom = element_text(size = 7), # Adjust y-axis text size
      strip.text.y.left = element_text(angle = 0, size = 7), # Left-aligned facet labels
      strip.text.x = element_text(size = 7)
    ) + 
    facet_grid(species ~ relative2, switch = "y") # Move facet labels to the left
  ggsave("output/2d_hist.pdf", width = 3, height = 7)
  
  
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
  
  
  
  df_angle <- subset(df_relative_normalized, relative2 == "female")
  df_subset <- df_angle[seq(1, nrow(df_angle), by = 10), ]
  
  acf(df_subset$angle, lag.max = 30)
  
  library(circular)
  
  df_temp <- subset(df_relative_normalized, relative2 == "female" & dis < 2)
  
  watson.test(circular(subset(df_temp, species =="Gly_fus")$angle))
  watson.test(circular(subset(df_temp, species =="Gly_sat")$angle))
  watson.test(circular(subset(df_temp, species =="Gly_nak_sexual")$angle))
  watson.test(circular(subset(df_temp, species =="Gly_nak_asexual")$angle))
  
  # HR test is favorable but cannot complerete calculation in time. try with better PC
  # testdata = circular::rvonmises(20, mu = circular::circular(pi), kappa = 3)
  # HR_test(testdata, iter = 999)
  
  # Wheeler-Watson Test
  
  library(circular)
  
  # Prepare circular data
  
  # Initialize results data frame
  results <- data.frame(pair = character(), statistic = numeric(), 
                        p_value = numeric(), stringsAsFactors = FALSE)
  
  # Perform Watson Two-Sample Test for all pairwise combinations using a for loop
  pairs <- combn(names(angles), 2, simplify = FALSE)
  
  for (pair in pairs) {
    print(pair)
    test_result <- watson.two.test(angles[[pair[1]]], angles[[pair[2]]])
    results <- rbind(results, data.frame(
      pair = paste(pair[1], "vs", pair[2]),
      statistic = test_result$statistic,
      p_value = test_result$p.value
    ))
  }
  
  # Print results
  print(results)
  
  
  angle_gly_fus <- circular(subset(df_temp, species =="Gly_fus")$angle)
  angle_gly_sat <- circular(subset(df_temp, species =="Gly_sat")$angle)
  angle_gly_nak_sexual <- circular(subset(df_temp, species =="Gly_nak_sexual")$angle)
  angle_gly_nak_asexual <- circular(subset(df_temp, species =="Gly_nak_asexual")$angle)
  angles <- list(
    Gly_fus = circular(subset(df_temp, species == "Gly_fus")$angle),
    Gly_sat = circular(subset(df_temp, species == "Gly_sat")$angle),
    Gly_nak_sexual = circular(subset(df_temp, species == "Gly_nak_sexual")$angle),
    Gly_nak_asexual = circular(subset(df_temp, species == "Gly_nak_asexual")$angle)
  )
  
  
  angles[[1]]
  watson.two.test(angle_gly_nak_sexual, angle_gly_nak_asexual)
  watson.two.test(angle_gly_fus, angle_gly_sat)
  watson.two.test(angle_gly_sat, angle_gly_nak_sexual)
  watson.two.test(angle_gly_sat, angle_gly_nak_asexual)
  watson.two.test(angle_gly_sat, angle_gly_nak_sexual)
  
  ggplot(subset(df_relative_normalized, relative2 == "female"), 
         aes(x=dis, y =species))+
    geom_density_ridges()
  
  ggplot(subset(df_relative_normalized, relative2 == "female" & dis < 10), 
         aes(x = dis, y = forcats::fct_rev(species))) +
    geom_density_ridges(stat = "binline", bins = 40, scale = 0.7) +
    labs(y = "", x = "Relaive angle from female (rad)") + 
    scale_x_continuous(
      breaks = c(0, 2),
      labels = c("0", "2")
    ) +
    geom_vline(xintercept = 2) + 
    theme_classic()
  
  
}
# ---------------------------------------------------------------------------- #

load("data_fmt/df_analysis.rda")
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

df_tandem_analysis$status <- factor(df_tandem_analysis$status)
df_tandem_analysis[df_tandem_analysis$species == "Gly_nak_asexual" &
                     df_tandem_analysis$status == "mleader",]$status <- "fleader"


df_proportions <- df_tandem_analysis %>%
  group_by(species, video, status) %>%
  tally() %>% 
  ungroup() %>%
  group_by(video) %>%
  mutate(proportion = n / sum(n)) %>% 
  ungroup()

ggplot(df_proportions, aes(x = video, y = proportion, 
                           fill = status, col=species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_color_viridis(discrete = T)

library(dplyr)
library(tidyr)


df_proportions$sep_proportion <- 
  rep(df_proportions[df_proportions$status == "sep",]$proportion, 
      rle(df_proportions$video)$length)

df_proportions2 <- df_proportions %>%
  arrange(species, sep_proportion) %>%
  mutate(video = factor(video, levels = unique(video)))

ggplot(df_proportions2, aes(x = video, y = proportion, fill = status)) +
  geom_bar(stat = "identity", position = "stack", linewidth = .1, color = NA) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top")

df_wide <- df_proportidiscrete = df_wide <- df_proportions %>%
  select(species, video, status, proportion) %>%
  spread(key = status, value = proportion, fill = 0)

df_wide$tandem <- df_wide$fleader + df_wide$mleader

df_wide$species <- factor(df_wide$species)
ggplot(df_wide) +
  geom_jitter(aes(x = as.numeric(species), y = tandem), position = position_dodge(width = 0))+
  geom_jitter(aes(x = as.numeric(species)+.1, y = mleader), position = position_dodge(width = 0.3), col = 2) +
  geom_jitter(aes(x = as.numeric(species)+.2, y = fleader), position = position_dodge(width = 0.6), col = 3)

ggplot(df_wide, aes(x = fleader, y = mleader, col=species)) +
  geom_point()



ggplot(df_wide) +
  geom_jitter(aes(x = species, y = tandem), position = position_jitter(width = 0.1)) +
  geom_jitter(aes(x = species, y = mleader), position = position_jitter(width = 0.1), col = 2, 
              aes(x = as.numeric(species) + 0.3)) +
  geom_jitter(aes(x = species, y = fleader), position = position_jitter(width = 0.1), col = 3, 
              aes(x = as.numeric(species) + 0.6))

  


geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) 
  labs(x = "Species", y = "Proportion", fill = "Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


  
  df_combined <- df_proportions %>%
  group_by(species) %>%
  mutate(
    fleader_mleader = ifelse(status %in% c("fleader", "mleader"), proportion, 0)
  ) %>%
  ungroup()

# Sum the combined proportion for each species (fleader + mleader)
df_combined_sum <- df_combined %>%
  group_by(species) %>%
  summarise(fleader_mleader_proportion = sum(fleader_mleader)) %>%
  ungroup()

# Add the combined category to the original df_proportions
df_proportions_combined <- df_proportions %>%
  bind_rows(df_combined_sum %>% mutate(status = "fleader+mleader", proportion = fleader_mleader_proportion))

# Plot the proportions
ggplot(df_proportions_combined, aes(x = species, y = proportion, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Species", y = "Proportion", fill = "Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




df_combined <- df_proportions %>%
  group_by(species) %>%
  mutate(
    fleader_mleader = ifelse(status %in% c("fleader", "mleader"), proportion, 0)
  ) %>%
  ungroup()

# Sum the combined proportion for each species (fleader + mleader)
df_combined_sum <- df_combined %>%
  group_by(species) %>%
  summarise(fleader_mleader_proportion = sum(fleader_mleader)) %>%
  ungroup()

# Add the combined category to the original df_proportions
df_proportions_combined <- df_proportions %>%
  bind_rows(df_combined_sum %>% mutate(status = "fleader+mleader", proportion = fleader_mleader_proportion))

# Plot the proportions
ggplot(df_proportions_combined, aes(x = species, y = proportion, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Species", y = "Proportion", fill = "Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ----------------------------------------------------------------------------------------- #
# 2D Plot
# ----------------------------------------------------------------------------------------- #
{
  plot_2D_density <- function(df, relative2 = "female", plot.range = 15, xlab, ylab){
    df <- subset(df, df$frame %% 30 == 0)
    
    if(relative2 == "female"){
      dir_vec_x = df$fHead_x - df$fTip_x
      dir_vec_y = df$fHead_y - df$fTip_y
      rel_vec_x = df$fCenter_x - df$mCenter_x
      rel_vec_y = df$fCenter_y - df$mCenter_y
    } else if (relative2 == "male") {
      dir_vec_x = df$mHead_x - df$mTip_x
      dir_vec_y = df$mHead_y - df$mTip_y
      rel_vec_x = df$mCenter_x - df$fCenter_x
      rel_vec_y = df$mCenter_y - df$fCenter_y
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
      
    )
    
    df_plot_temp <- data.frame(
      x = cos(sta_ang)*dis,
      y = sin(sta_ang)*dis
    )
  
    p <- ggplot(subset(df_plot_temp, abs(x) < plot.range & abs(y) < plot.range), aes(x=x, y=y))+
      geom_bin2d(bins = 40) +
      scale_fill_viridis(option = "viridis", guide = guide_colorbar(barwidth = .5, barheight = 3)) +
      scale_x_continuous(expand = c(0, 0.01), limits = c(-plot.range, plot.range)) +
      scale_y_continuous(expand = c(0, 0.01), limits = c(-plot.range, plot.range)) +
      theme_classic() +
      coord_fixed() +
      theme(legend.position = "right", 
            legend.title = element_blank(),
            legend.text = element_text(size = 6)) +
      xlab(xlab) +
      ylab(ylab) 
    
    return(p)
  }
  {
    df <- arrow::read_feather("data_fmt/Gly_sat_df.feather")
    p1 <- plot_2D_density(df, "female", 15, "", "Distance back-front (body length)")
    p2 <- plot_2D_density(df, "male", 15, "", "")
    
    df <- arrow::read_feather("data_fmt/Gly_fus_df.feather")
    p3 <- plot_2D_density(df, "female", 10, "", "Distance back-front (body length)")
    p4 <- plot_2D_density(df, "male", 10, "", "")
    
    df <- arrow::read_feather("data_fmt/Gly_nak_sexual_df.feather")
    df <- subset(df, video != "Gly_nak_NM2344_9") # the pair male did not move for 1 h and female was next to it.
    p5 <- plot_2D_density(df, "female", 10, "", "Distance back-front (body length)")
    p6 <- plot_2D_density(df, "male", 10, "Distance left-right (body length)", "")
    
    df <- arrow::read_feather("data_fmt/Gly_nak_asexual_df.feather")
    p7 <- plot_2D_density(df, "both", 10, "Distance left-right (body length)", "Distance back-front (body length)")
    
    gridExtra::grid.arrange(p3, p4, p1, p2, p5, p6, p7, nrow = 4, ncol = 2)
    ggsave("output/2Dhist.pdf")
  }
}
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# pair distance (smaller than random distribution?)
# ----------------------------------------------------------------------------------------- #
dis_return <- function(df, species){
  df <- subset(df, df$frame %% 30 == 0)
  rel_vec_x = df$fCenter_x - df$mCenter_x
  rel_vec_y = df$fCenter_y - df$mCenter_y
  dis <- sqrt(rel_vec_x^2 + rel_vec_y^2)
  return(data.frame(dis, species))
}

df_dis <-NULL
df <- arrow::read_feather("data_fmt/Gly_sat_df.feather")
df_dis <- rbind(df_dis, dis_return(df, "Gly_sat"))

df <- arrow::read_feather("data_fmt/Gly_fus_df.feather")
df_dis <- rbind(df_dis, dis_return(df, "Gly_fus"))

df <- arrow::read_feather("data_fmt/Gly_nak_sexual_df.feather")
df <- subset(df, video != "Gly_nak_NM2344_9") # the pair male did not move for 1 h and female was next to it.
df_dis <- rbind(df_dis, dis_return(df, "Gly_nak_sexual"))

df <- arrow::read_feather("data_fmt/Gly_nak_asexual_df.feather")
df_dis <- rbind(df_dis, dis_return(df, "Gly_nak_asexual"))

library(ggridges )
ggplot(df_dis, aes(x=dis, y = species)) + 
  geom_density_ridges()+
  xlim(c(0,20))
# ----------------------------------------------------------------------------------------- #



df <- arrow::read_feather("data_fmt/Gly_sat_df.feather")
df <- arrow::read_feather("data_fmt/Gly_fus_df.feather")

{
  
  video_names <- unique(df$video)
  
  dir_vec_x = df$fHead_x - df$fTip_x
  dir_vec_y = df$fHead_y - df$fTip_y
  rel_vec_x = df$fCenter_x - df$mCenter_x
  rel_vec_y = df$fCenter_y - df$mCenter_y
  sta_ang_ftom <- atan2(rel_vec_y, rel_vec_x) - atan2(dir_vec_y, dir_vec_x) + pi/2
  dis <- sqrt(rel_vec_x^2 + rel_vec_y^2)
  
  dir_vec_x = df$mHead_x - df$mTip_x
  dir_vec_y = df$mHead_y - df$mTip_y
  rel_vec_x = df$mCenter_x - df$fCenter_x
  rel_vec_y = df$mCenter_y - df$fCenter_y
  sta_ang_mtof <- atan2(rel_vec_y, rel_vec_x) - atan2(dir_vec_y, dir_vec_x) + pi/2
  
  sta_ang_ftom[sta_ang_ftom > 2*pi] <- sta_ang_ftom[sta_ang_ftom > 2*pi] - 2*pi
  sta_ang_ftom[sta_ang_ftom < -2*pi] <- sta_ang_ftom[sta_ang_ftom < -2*pi] + 2*pi
  
  sta_ang_mtof[sta_ang_mtof > 2*pi] <- sta_ang_mtof[sta_ang_mtof > 2*pi] - 2*pi
  sta_ang_mtof[sta_ang_mtof < -2*pi] <- sta_ang_mtof[sta_ang_mtof < -2*pi] + 2*pi
  
  sta_ang_ftom[sta_ang_ftom > pi] <- sta_ang_ftom[sta_ang_ftom > pi] - 2*pi
  sta_ang_ftom[sta_ang_ftom < -pi] <- sta_ang_ftom[sta_ang_ftom < -pi] + 2*pi
  
  sta_ang_mtof[sta_ang_mtof > pi] <- sta_ang_mtof[sta_ang_mtof > pi] - 2*pi
  sta_ang_mtof[sta_ang_mtof < -pi] <- sta_ang_mtof[sta_ang_mtof < -pi] + 2*pi
  
  df$leader <- "none"
  df$leader[sta_ang_ftom > 0 & sta_ang_mtof < 0 &  dis < mean(tbl)*2] <- "male"
  df$leader[sta_ang_ftom < 0 & sta_ang_mtof > 0 &  dis < mean(tbl)*2] <- "female"
  
  videos <- unique(df$video)
  
  df_tandem <- NULL
  for(i_v in 1:length(videos)){
    df_temp <- subset(df, video == videos[i_v])
    rle_leader <- rle(df_temp$leader)
    
    fleader_duration <- rle_leader$lengths[rle_leader$values == "female"]
    mleader_duration <- rle_leader$lengths[rle_leader$values == "male"]
    
    leader_sums <- list()
    start <- 1
    
    f_dis <- sqrt( diff(df_temp$fCenter_x)^2 + diff(df_temp$fCenter_y)^2 )
    m_dis <- sqrt( diff(df_temp$mCenter_x)^2 + diff(df_temp$mCenter_y)^2 )
    
    for (i in seq_along(rle_leader$lengths)) {
      leader_type <- rle_leader$values[i]
      run_length <- rle_leader$lengths[i]
      
      if (leader_type != "none") { 
        end <- start + run_length - 1
        f_dis_temp <- sum(f_dis[start:end])
        m_dis_temp <- sum(m_dis[start:end])
        df_tandem <- rbind(df_tandem, data.frame(
          video = videos[i_v],
          leader = leader_type,
          duration = run_length,
          f_dis = f_dis_temp,
          m_dis = m_dis_temp
        ))
      }
      start <- start + run_length
    }
    
    
    
  }
}

df_tandem



  
  ggplot(subset(df_tandem, duration/30 > 0 & f_dis > 0), 
         aes(x = duration/30, y = f_dis, col=leader))+
    geom_point() +
    facet_grid(~leader) +
    scale_x_log10() + scale_y_log10()
  
  plot(tapply(df$leader == "female", df$video, sum), tapply(df$leader == "male", df$video, sum))

  library(survival)
  library(survminer)
  library(car)
  library(coxme)
  m <- coxme(Surv(duration) ~ leader + (1|video), data = df_tandem)
  Anova(m)
  r<-survfit(Surv((duration/30)) ~leader, type = "kaplan-meier", data=df_tandem)
  ggsurvplot(fit = r, data = df_tandem,
             pval = F, pval.method = TRUE,
             risk.table = F, conf.int = FALSE,
             ncensor.plot = FALSE, size = 1, linetype = 1:3,
             xlab="Time (min)", ggtheme = theme_bw()  + theme(aspect.ratio = 0.75))
  
  



df_tandem <- subset(df_r, fhead_mhead_dis < 2)

pca_res <- prcomp(df_tandem[,2:10])
summary(pca_res)

df_tandem <- cbind(df_tandem, pca_res$x[,1:2])

library(ggplot2)
library(viridis)
ggplot(df_tandem,aes(x=PC1, y=PC2)) + 
  stat_density_2d(geom = "polygon",  aes(alpha = ..level..))+
  scale_fill_viridis(discrete = T, direction = 1, end=0.5)

  coord_fixed(ylim = c(-1,1), xlim = c(-2,2), expand = F) +
  theme(aspect.ratio = 1)+
  theme_classic()

c1 <- kmeans(df_tandem[,c("PC1","PC2")], center=6)

df_tandem$cluster = c1$cluster

ggplot(df_tandem,aes(x=PC1, y=PC2)) + 
  stat_density_2d(geom = "polygon",  aes(alpha = ..level.., fill=as.factor(cluster)))+
  scale_fill_viridis(discrete = T, trans = "log10")

ggplot(df_tandem[sample(1:dim(df_tandem)[1],10000),], aes(x=PC1, y=PC2, col=as.factor(cluster))) + 
  geom_point(alpha=.2)+
  scale_color_viridis(discrete = T)




df_r$fhead_mtip_dis
dis_change <- c(NA, diff(df_r$fhead_mtip_dis))
dis_change[df_r$frame==0] <- NA
df_r$dis_change = c(dis_change[2:length(dis_change)], NA)
a <- tapply(df_r$dis_change, round(df_r$fhead_mtip_dis,2), mean, na.rm=T)
plot(as.numeric(names(a)), a, ylim=c(-0.02,0.02), xlim = c(0,4), type="o")
abline(h=0)

dis_change <- c(NA, diff(df_r$ftip_mhead_dis))
dis_change[df_r$frame==0] <- NA
df_r$dis_change <- c(dis_change[2:length(dis_change)], NA)
a <- tapply(df_r$dis_change, round(df_r$ftip_mhead_dis,2), mean, na.rm=T)
plot(as.numeric(names(a)), a, ylim=c(-0.02,0.02), xlim = c(0,4), type="o")
abline(h=0)
