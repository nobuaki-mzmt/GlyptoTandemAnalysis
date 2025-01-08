library(arrow)
library(stringr)
df <- arrow::read_feather("data_fmt/Gly_sat_df.feather")

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
    
    gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, nrow = 4)
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
