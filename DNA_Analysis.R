library(ggplot2)
library(vegan)
library(tidyverse)
library(lme4)
library(emmeans)
library(car)

#load in the metadata
dna_metadata<- read.csv("F:\\Research\\16S_Soil\\formatted_metadata.csv")
  
#load in the dna amounts data
dna_amount<- read.csv("F:\\Research\\16S_Soil\\DNA_amounts.csv")
#need to make sample id a character
dna_amount$SampleID<- as.character(dna_amount$SampleID)
str(dna_amount)

#load in waypoints
waypoints<- read.csv("F:\\Research\\waypoints.csv")

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
ggplot(dna_rgm, aes(treatment.x, DNA_ug_per_g)) +
  geom_boxplot(alpha = 0.5, aes(fill=treatment.x), outliers=F) +
  geom_jitter(aes(color=treatment.x))+
  labs(x = NULL, y = "DNA ug/g", title = "DNA Amounts for RGM Dry") +
  scale_fill_manual(values=c('cyan3','purple3'))+
  scale_color_manual(values=c('cyan', 'purple1'))+
  theme_bw() 
