library(ggplot2)
library(vegan)
library(plyr)
library(tidyverse)
library(lme4)
library(emmeans)
library(car)
library(patchwork)

theme<-theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
             axis.title.y = element_text(face="bold", size = 18, margin = margin(r = 10)), 
             axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16),
             axis.title.x = element_text(size = 18, face = "bold",color = "black"),
             plot.margin = unit(c(0.1,0.1,0,0.1),"cm") ,
             strip.text = element_text(face = "bold", size = 16), 
             legend.position = "none")

dist_cols <- colorRampPalette(c("#9C6EB0", "cyan3"))(6)
dist_cols

#load in data
dna_data<- read.csv("DNA_amounts_distanceData.csv")
dna_data$distance<- as.factor(dna_data$distance)

#filter out NA columns
dna_data <- dna_data %>% 
  filter(!is.na(DNA_amount))

#run model for just latrines
#make data with out veg patch
dna_lat<- dna_data %>% 
  filter(type=='latrine')

m_dna_dis<- lmer(DNA_amount~distance+(1|latrine), data=dna_lat)
summary(m_dna_dis)
Anova(m_dna_dis)
emmeans(m_dna_dis, pairwise~distance)
# 1 was different from everything else but everything else was very similar 
# to each other

#model with inside and outside factor
dna_lat$inside<- as.factor(dna_lat$inside)
m_dna_ins<- lmer(DNA_amount ~inside+(1|latrine), data=dna_lat)
summary(m_dna_ins)
Anova(m_dna_ins)

#model with ins, ous, and veg patch
dna_data$type_ins<- as.factor(dna_data$type_ins)
m_dna_veg<- lmer(DNA_amount~type_ins+(1|latrine), data=dna_data)
summary(m_dna_veg)
Anova(m_dna_veg)
emmeans(m_dna_veg, pairwise~type_ins)
#inside and veg patch were similar, and each were different from outside latrine

#make plots for each

#plot for latrine distance
#make mean and s values
#we have to make a wide data frame so it has mean and SE for each 
dna_latSum<-dna_lat %>%
  group_by(distance) %>%
  summarise_at(vars(DNA_amount), list(mean_dis=mean, sd_dis=sd)) %>% 
  mutate(se_dis=sd_dis/sqrt(36)) %>% 
  mutate(distance=recode_factor(distance, '1'='Latrine center','2'='Latrine edge','3'='0.5m', '4'='1m', '5'='2m','6'='4m')) %>% 
  as.data.frame()

dna_lat <- dna_lat %>% 
  mutate(dis_name=recode_factor(distance, '1'='Latrine center','2'='Latrine edge','3'='0.5m', '4'='1m', '5'='2m','6'='4m')) 

#background colors
bg <- data.frame(
  xmin = c(0.5, 2.5),
  xmax = c(2.5, 6.5),
  ymin = -Inf,
  ymax = Inf,
  fill = c("#9C6EB0", "cyan3"))
labels_bg <- data.frame(
  x = c(1.5, 4.5),
  y = Inf,
  label = c("Inside Latrine", "Outside Latrine"))

#plot it
dna_dists <- ggplot(dna_lat, aes(x = dis_name, y = DNA_amount)) +
  geom_rect(data = bg,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.1) +
  scale_fill_identity() +
  ggnewscale::new_scale_fill() +
  geom_boxplot(aes(fill = dis_name),
    outliers = FALSE) +
  geom_point(color = "black",
    alpha = 0.5,
    size = 3) +
  geom_text(
    data = labels_bg,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 4,
    size = 5,
    fontface = "bold",
    color = c("#5B3A7A", "#1A8FA3")) +
  scale_fill_manual(values = dist_cols, guide = "none") +
  labs(x = NULL, y = "DNA (µg/g)", tag = "(g)") +
  theme_bw() +
  theme +
  theme(plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.98),
    legend.position = "none")
dna_dists


#plot inside outside
dna_insSum<-dna_lat %>%
  group_by(inside) %>%
  summarise_at(vars(DNA_amount), list(mean_ins=mean, sd_ins=sd)) %>% 
  mutate(se_ins=sd_ins/sqrt(36)) %>% 
  as.data.frame()

ggplot(dna_insSum, aes(x=inside, y=mean_ins, color=inside))+
  geom_errorbar(aes(ymin=mean_ins-se_ins, ymax=mean_ins+se_ins))+
  geom_point(data=dna_lat, aes(x=inside, y=DNA_amount, color=inside), size=3, alpha=.5)+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))+
  theme_bw()

#plot with veg patch
dna_vegSum<-dna_data %>%
  group_by(type_ins) %>%
  summarise_at(vars(DNA_amount), list(mean=mean, sd=sd)) %>% 
  mutate(se=sd/sqrt(36)) %>% 
  as.data.frame()
dna_vegSum$type_ins <- as.factor(dna_vegSum$type_ins)

dna_type <- dna_type <- ggplot(dna_data, aes(type_ins, DNA_amount)) +
  geom_boxplot(
    aes(fill = type_ins),
    outliers = FALSE) +
  geom_point(
    aes(color = type_ins),
    alpha = 0.5,
    size = 3) +
  labs(
    x = NULL,
    y = "DNA (µg/g)",
    tag = "(h)") +
  scale_fill_manual(
    values = c("#9C6EB0","#4ACD66","#4DD7CE"),
    guide = "none") +
  scale_color_manual(
    values = c("black", "black", "black")) +
  scale_x_discrete(
    labels = c(
      "ins_latrine" = "Inside Latrine",
      "ins_veg_patch" = "Vegetation Patch",
      "out_latrine" = "Outside Latrine") ) +
  theme_bw() +
  theme +
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.98),
    legend.position = "none")
dna_type

#put the two together for manuscript 
dna_dec <- dna_dists + dna_type
dna_dec
#export 1800 by 600

#multi combined plot 
dna_and_nuts / dna_dec
#export 2000 x 1600