d <- read.delim("clipboard", header =T)

library(ggplot2)
library(ggridges)

ggplot(d, aes(y = pq + pk, x = species, col = development, fill = development))+
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2, alpha=0.5)+
  facet_grid(development~.)

ggplot(subset(d, species == "nakajimai_asex"),
       aes(x = pq + pk + sq + sk, y = species, col = development, fill = development))+
  #geom_density_ridges(stat = "binline", bins = 40, scale = 0.7, alpha = .4)
  geom_density_ridges(alpha = .4)

table(d_g$p_total)
d_g$p_total[d_g$p_total > 10] <- 11
ggplot(d_g,
       aes(x = p_total, y = species, col = development, fill = development))+
        geom_density_ridges(stat = "binline", bins = 40, scale = 0.5, alpha = .4) +
  facet_grid(~development)

ggplot(subset(d_g, development == "mature" & species != "satsumensis"),
       aes(x = p_total, y = species, col = development, fill = development))+
  geom_density_ridges(stat = "binline", bins = 80, scale = 0.5, alpha = .4) +
  geom_vline(xintercept = 2)



d_g <- subset(d, genus == "Glyptotermes")
d_g$r_total <- d_g$p_total + d_g$sk + d_g$sq
d_g$r_total[d_g$r_total > 11] <- 11
ggplot(d_g,
       aes(x = r_total, y = 1, col = species, fill = species))+
  geom_density_ridges(stat = "binline", bins = 40, scale = 0.5, alpha = .4) +
  scale_fill_viridis(discrete =T )+
  facet_grid(~development)

d_g <- subset(d, genus == "Glyptotermes")
ggplot(d_g, aes(y = r_total, x = development, fill = species))+
  scale_fill_viridis(discrete =T )+
  geom_dotplot(binaxis = "y", binwidth = 1, alpha=0.9, dotsize = .2, stackgroups = TRUE)

d_g$p_total[d_g$p_total > 11] <- 11
ggplot(subset(d_g, p_total >= 2), 
       aes(y = p_total, x = development, fill = species))+
  scale_fill_viridis(discrete =T )+
  geom_dotplot(binaxis = "y", binwidth = 1, alpha=0.9, dotsize = .2, 
               stackgroups = TRUE, binpositions = "all") +
  scale_y_continuous(limits = c(2,11), breaks = c(2,4,6,8,10,11), 
                     labels = c("2", "4", "6", "8", "10", ">11")) +
  labs(x = "", y = "Number of primary reproductives") + 
  theme_classic()

d_temp <- subset(d_g, p_total >= 2)
tapply(d_temp$p_total, d_temp[,c("species", "development")], mean)

qd_g <- subset(d, genus == "Glyptotermes")
ggplot(d_g, aes(y = pq + pk, x = species, col = development, fill = development))+
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2, alpha=0.5)+
  facet_grid(development~.)

d_g <- subset(d, genus == "Glyptotermes")
ggplot(d_g, aes(y = pq + pk, x = species, col = development, fill = development))+
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2, alpha=0.5)+
  facet_grid(development~.)

library(lme4)
library(car)
library(multcomp)
d_g$species <- factor(d_g$species)
r <- glm(p_total ~ species, family = "poisson", data = subset(d_g, development == "mature"))
Anova(r)
multicomparison<-glht(r, linfct = mcp(species = "Tukey"))
summary(multicomparison)

r <- glm(p_total ~ species, family = "poisson", data = subset(d_g, development == "incipient"))
Anova(r)
multicomparison<-glht(r, linfct = mcp(species = "Tukey"))
summary(multicomparison)

r <- glm(p_total ~ development, family = "poisson", data = subset(d_g, species == "fuscus"))
Anova(r)

r <- glm(p_total ~ development, family = "poisson", data = subset(d_g, species == "nakajimai_asex"))
Anova(r)
