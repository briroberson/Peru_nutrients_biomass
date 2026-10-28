library(ggplot2)
library(vegan)
library(tidyverse)
library(lme4)
library(emmeans)
library(car)
library(patchwork)

theme<-theme(plot.title = element_text(size = 16, hjust = 0.5),
             axis.title.y = element_text(face="bold", size = 18), 
             axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16),
             axis.title.x = element_text(size = 18, face = "bold",color = "black"),
             plot.margin = unit(c(0.1,0.1,0,0.1),"cm") ,
             strip.text = element_text(face = "bold", size = 16), 
             legend.position = "none")


#load in the metadata
dna_metadata<- read.csv("formatted_metadata.csv")

#load in the dna amounts data
dna_amount<- read.csv("DNA_amounts_BigData.csv")
#need to make sample id a character
dna_amount$SampleID<- as.character(dna_amount$SampleID)
str(dna_amount)

#load in waypoints
waypoints<- read.csv("waypoints.csv")

#format latrine names
waypoints$latrineF<-gsub("^.{0,4}", "", waypoints$latrine)
waypoints$latrine<- paste("L", waypoints$latrineF, sep='') 

#merge them together and drop some duplicate columns 
dna_data<- dna_metadata %>% 
  left_join(dna_amount, by='SampleID') %>% 
  dplyr::select(-location, -treatment.y, -season, -year, -replicate.y)

#drop NAs AND the 3 latrines that were probably swapped         
dna_data <- dna_data %>% 
  filter(!is.na(DNA_ug_per_g)) %>% 
  filter(!latrine_trt_month %in% c("L83_latrine_wet", 'L83_control_wet', "L94_latrine_wet", 'L94_control_wet', "L62_control_wet", 'L62_latrine_wet'))

#merge waypoints
dna_data<- dna_data %>% 
  left_join(waypoints, by='latrine') 

############## run models

###########subset data for lia and rgm wet
dna_wet<- dna_data %>% 
  filter(month.collected=='wet')

m_soilAge<-lmer(DNA_ug_per_g~treatment.x*soilAge+treatment.x*scale(elevation)+(1|latrine_trt_month)+(1|latrine), data=dna_wet)
summary(m_soilAge)
Anova(m_soilAge, type='III')
emmeans(m_soilAge, pairwise~treatment.x*soilAge)

#check model assumptions
qqnorm(residuals(m_soilAge)) #checking normality


#plot the DNA amounts to examine
ggplot(dna_wet, aes(treatment.x, DNA_ug_per_g)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment.x), outliers=F) +
  geom_jitter(aes(color=treatment.x))+
  labs(x = NULL, y = "DNA ug/g", title = "DNA Amounts by Soil Age") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_manual(values=c('cyan', 'purple1'))+
  theme_bw() +
  facet_wrap(~soilAge)

####### models for RGM wet data
dna_rgmW<- dna_data %>% 
  filter(soilAge=='rgm' & month.collected=='wet')

m_rgmW<- lmer(DNA_ug_per_g~treatment.x*scale(elevation)+(1|latrine_trt_month)+(1|latrine), data=dna_rgmW)
Anova(m_rgmW, type='III')
summary(m_rgmW)

#check model assumptions
qqnorm(residuals(m_rgmW)) #checking normality


#plot DNA by season
ggplot(dna_rgmW, aes(treatment.x, DNA_ug_per_g)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment.x), outliers=F) +
  geom_jitter(aes(color=treatment.x))+
  labs(x = NULL, y = "DNA ug/g", title = "DNA Amounts for RGM Wet") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_manual(values=c('cyan', 'purple1'))+
  theme_bw() 

#plot elevation just to see
ggplot(dna_rgmW, aes(x=elevation, y=DNA_ug_per_g))+
  geom_point()+
  geom_smooth(method='lm')

####### models for LIA wet data
dna_lia<- dna_data %>% 
  filter(soilAge=='lia')

m_lia<- lmer(DNA_ug_per_g~treatment.x*scale(elevation)+(1|latrine_trt_month)+(1|latrine), data=dna_lia)
Anova(m_lia, type='III')
summary(m_lia)

#check model assumptions
qqnorm(residuals(m_lia)) #checking normality


#plot DNA by season
ggplot(dna_lia, aes(treatment.x, DNA_ug_per_g)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment.x), outliers=F) +
  geom_jitter(aes(color=treatment.x))+
  labs(x = NULL, y = "DNA ug/g", title = "DNA Amounts for LIA Wet") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_manual(values=c('cyan', 'purple1'))+
  theme_bw() 


####### models for RGM dry data
dna_rgmD<- dna_data %>% 
  filter(soilAge=='rgm' & month.collected=='dry')

m_season<- lmer(DNA_ug_per_g~treatment.x*scale(elevation)+(1|latrine_trt_month)+(1|latrine), data=dna_rgmD)
Anova(m_season, type='III')
summary(m_season)

#check model assumptions
qqnorm(residuals(m_season)) #checking normality


#plot DNA by season
ggplot(dna_rgmD, aes(treatment.x, DNA_ug_per_g)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment.x), outliers=F) +
  geom_jitter(aes(color=treatment.x))+
  labs(x = NULL, y = "DNA ug/g", title = "DNA Amounts for RGM Dry") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_manual(values=c('cyan', 'purple1'))+
  theme_bw() 



#extract significance for plotting
tests <- list(
  wet_LIA = m_lia, 
  wet_RGM = m_rgmW,
  dry_RGM = m_season)

p_to_star <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns")}

DNA_sig_table <- data.frame(
  variable = names(tests),
  p.value = sapply(tests, function(x) {
    a <- car::Anova(x, type = "III")
    a["treatment.x", "Pr(>Chisq)"]
  }))

DNA_sig_table$stars <- p_to_star(DNA_sig_table$p.value)

DNA_sig_table$panel <- dplyr::case_when(
  DNA_sig_table$variable == "wet_LIA" ~ "lia",
  DNA_sig_table$variable == "wet_RGM" ~ "rgm",
  DNA_sig_table$variable == "dry_RGM" ~ "dry_RGM")

DNA_sig_table$panel <- factor(
  DNA_sig_table$panel,
  levels = c("lia", "rgm", "dry_RGM"))





#new combined manuscript plot 

dna_all <- bind_rows(
  dna_wet %>% 
    mutate(panel = soilAge),  
  dna_rgmD %>% 
    mutate(panel = "dry_RGM")) 

dna_all$panel <- factor(
  dna_all$panel,
  levels = c("lia", "rgm", "dry_RGM"))

dna_plot <- ggplot(dna_all, aes(treatment.x, DNA_ug_per_g)) +
  geom_boxplot(aes(fill = treatment.x), outliers = FALSE) +
  geom_point(aes(color = treatment.x), alpha = 0.5) +
  labs(x = NULL, y = "DNA (µg/g)") +
  scale_fill_manual(
    values = c("cyan3", "#9C6EB0"),
    guide = "none") +
  scale_color_manual(
    values = c("black", "black")) +
  scale_x_discrete(labels = c(
    control = "Reference",
    latrine = "Latrine")) +
  facet_wrap( ~panel, nrow = 1,
              labeller = labeller(panel = c(
                lia = " Wet season LIA",
                rgm = "Wet season RGM",
                dry_RGM = "Dry season RGM"))) +
  geom_text(data = DNA_sig_table,
            aes(x = 1.5, y = Inf, label = stars),
            inherit.aes = FALSE,
            vjust = 1.2,
            size = 6) + 
  theme_bw() + theme +
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.98)) 

dna_plot

#1200x600


#combined with nutrient plots
design <- c(
  area(t = 1, l = 1, b = 1, r = 1),  # both_c
  area(t = 1, l = 2, b = 1, r = 2),  # both_cn
  area(t = 2, l = 1, b = 2, r = 1),  # both_N
  area(t = 2, l = 2, b = 2, r = 2),  # both_d13C
  area(t = 3, l = 1, b = 3, r = 1),  # both_d15N
  area(t = 1, l = 3, b = 3, r = 3)   # dna_plots w own column
  )

final_plot <-
  both_c +
  both_N +
  both_d13C +
  both_d15N +
  both_cn +
  dna_plot +
  plot_layout(design = design)

dna_and_nuts <- final_plot +
  plot_layout(
    design = design,
    widths = c(1, 1, 2))
dna_and_nuts
#export 2000 x 1000


#run DEC dna analysis for the rest











