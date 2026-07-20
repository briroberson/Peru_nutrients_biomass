library(ggplot2)
library(vegan)
library(plyr)
library(tidyverse)
library(lme4)
library(emmeans)
library(car)

#load in data
dna_data<- read.csv("F:\\Research\\Dec24_16S\\DNA_extraction.csv")
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
  mutate(distance=recode_factor(distance, '1'='in middle','2'='on edge','3'='.5m out', '4'='1m out', '5'='2m out','6'='4m out')) %>% 
  as.data.frame()

dna_lat <- dna_lat %>% 
  mutate(dis_name=recode_factor(distance, '1'='in middle','2'='on edge','3'='.5m out', '4'='1m out', '5'='2m out','6'='4m out')) 

#plot it
ggplot(dna_latSum, aes(x=distance, y=mean_dis, color=distance))+
  geom_errorbar(aes(ymin=mean_dis-se_dis, ymax=mean_dis+se_dis))+
  geom_point(data=dna_lat, aes(x=dis_name, y=DNA_amount, color=dis_name), size=3, alpha=.5)+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))+
  theme_bw()


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

ggplot(dna_vegSum, aes(x=type_ins, y=mean, color=type_ins))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se))+
  geom_point(data=dna_data, aes(x=type_ins, y=DNA_amount, color=type_ins), size=3, alpha=.5)+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))+
  theme_bw()

