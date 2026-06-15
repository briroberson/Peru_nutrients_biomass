library(ggplot2)
library(vegan)
library(tidyverse)
library(doBy)

#load in nutrient data
nutrient<- read.csv("Peru_Nutrient.csv")


########## LIA----

#filter for just lia and remove na
nut_LIA<- nutrient[c(1:4,6:8,11:13),]

#make data wide
wide_LIA<-pivot_wider(names_from='Type', values_from=c('mass_N','mol_N','Mass_C','mol_C','d15N','d13C','CN_Ratio'), data=nut_LIA)

#normality test on the difference between columns
shapiro.test(wide_LIA$mol_N_control-wide_LIA$mol_N_latrine)
shapiro.test(wide_LIA$mol_C_control-wide_LIA$mol_C_latrine)
shapiro.test(wide_LIA$d15N_control-wide_LIA$d15N_latrine)
shapiro.test(wide_LIA$d13C_control-wide_LIA$d13C_latrine)
shapiro.test(wide_LIA$CN_Ratio_control-wide_LIA$CN_Ratio_latrine)

#run paired ttest for control and latrine
t.test(Pair(mol_N_control,mol_N_latrine)~1, data=wide_LIA)
t.test(Pair(mol_C_control,mol_C_latrine)~1, data=wide_LIA)
t.test(Pair(d15N_control,d15N_latrine)~1, data=wide_LIA) #not different
t.test(Pair(d13C_control,d13C_latrine)~1, data=wide_LIA)
t.test(Pair(CN_Ratio_control,CN_Ratio_latrine)~1, data=wide_LIA) #not different

#plots
molN_stats <- summaryBy(mol_N ~ Type , nut_LIA, FUN = c(mean))

ggplot(nut_LIA, aes(x=Type, y=mol_N, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=molN_stats, aes(x = Type, y = mol_N.mean, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=mol_N), geom="errorbar", fun.data="mean_se",
                 fun.args=list(mult=1), size=1, width=.4) +
 scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

molC_stats <- summaryBy(mol_C ~ Type , nut_LIA, FUN = c(mean))

ggplot(nut_LIA, aes(x=Type, y=mol_C, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=molC_stats, aes(x = Type, y = mol_C.mean, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=mol_C), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

d13C_stats <- summaryBy(d13C ~ Type , nut_LIA, FUN = c(mean))

ggplot(nut_LIA, aes(x=Type, y=d13C, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=d13C_stats, aes(x = Type, y = d13C.mean, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=d13C), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

############ RGM----
#filter for just rgm
nut_RGM<- nutrient[14:27,1:10]

#make data wide
wide_RGM<-pivot_wider(names_from='Type', values_from=c('mass_N','mol_N','Mass_C','mol_C','d15N','d13C','CN_Ratio'), data=nut_RGM)

#normality test on the difference between columns
shapiro.test(wide_RGM$mol_N_control-wide_RGM$mol_N_latrine)
shapiro.test(wide_RGM$mol_C_control-wide_RGM$mol_C_latrine)
shapiro.test(wide_RGM$d15N_control-wide_RGM$d15N_latrine)
shapiro.test(wide_RGM$d13C_control-wide_RGM$d13C_latrine) #not normal
shapiro.test(wide_RGM$CN_Ratio_control-wide_RGM$CN_Ratio_latrine)

qqnorm(wide_RGM$d13C_control-wide_RGM$d13C_latrine)

#run paired ttest for control and latrine
t.test(Pair(mol_N_control,mol_N_latrine)~1, data=wide_RGM)
t.test(Pair(mol_C_control,mol_C_latrine)~1, data=wide_RGM)
t.test(Pair(d15N_control,d15N_latrine)~1, data=wide_RGM) #not different
t.test(Pair(d13C_control,d13C_latrine)~1, data=wide_RGM) # p= .05124, but not normal
t.test(Pair(CN_Ratio_control,CN_Ratio_latrine)~1, data=wide_RGM) 

RGM_d13C_c<- wide_RGM$d13C_control
RGM_d13C_l<- wide_RGM$d13C_latrine
wilcox.test(RGM_d13C_c, RGM_d13C_l, paired=T) # not different

#plots
RGM_molN_stats <- summaryBy(mol_N ~ Type , nut_RGM, FUN = c(mean))

ggplot(nut_RGM, aes(x=Type, y=mol_N, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=RGM_molN_stats, aes(x = Type, y = mol_N.mean, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=mol_N), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

RGM_molC_stats <- summaryBy(mol_C ~ Type , nut_RGM, FUN = c(mean))

ggplot(nut_RGM, aes(x=Type, y=mol_C, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=RGM_molC_stats, aes(x = Type, y = mol_C.mean, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=mol_C), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

RGM_CN_stats <- summaryBy(CN_Ratio ~ Type , nut_RGM, FUN = c(mean))

ggplot(nut_RGM, aes(x=Type, y=CN_Ratio, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=RGM_CN_stats, aes(x = Type, y = CN_Ratio.mean, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=CN_Ratio), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))



#####plot both soil types----

#mol C
RGM_C_mean<- data.frame(soilAge=c('rgm','rgm'), Type=c('control', 'latrine'), mol_C=c(231.7143, 26869))
LIA_C_mean<-data.frame(soilAge=c('lia','lia'), Type=c('control', 'latrine'), mol_C=c(6805.8, 28181.8))

ggplot(nutrient, aes(x=Type, y=mol_C, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=RGM_C_mean, aes(x = Type, y = mol_C, color=Type), size=5, alpha=.5) +
  geom_point(data=LIA_C_mean, aes(x = Type, y = mol_C, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=mol_C), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  facet_wrap(~soilAge)+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

#mol_N
RGM_N_mean<- data.frame(soilAge=c('rgm','rgm'), Type=c('control', 'latrine'), mol_N=c(21.37143, 1976.61429))
LIA_N_mean<-data.frame(soilAge=c('lia','lia'), Type=c('control', 'latrine'), mol_N=c(415.52, 1970.44))

nutrient[c(1:4,6:8,11:27),] %>% 
  ggplot(., aes(x=Type, y=mol_N, color=Type, group=Type))+
    geom_point(size=3)+
    geom_point(data=RGM_N_mean, aes(x = Type, y = mol_N, color=Type), size=5, alpha=.5) +
    geom_point(data=LIA_N_mean, aes(x = Type, y = mol_N, color=Type), size=5, alpha=.5) +
    stat_summary(aes(x=Type, y=mol_N), geom="errorbar", fun.data="mean_se",
                 fun.args=list(mult=1), size=1, width=.4) +
    scale_color_manual(values=c('cyan3','purple3'))+ 
    theme_bw()+
    facet_wrap(~soilAge)+
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title.y = element_text(face="bold", size = 18), 
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 18, face = "bold",color = "black"),
      axis.text.x=element_text(size=16),
      plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))


#CN ratio
RGM_CN_mean<- data.frame(soilAge=c('rgm','rgm'), Type=c('control', 'latrine'), CN_Ratio=c(mean(wide_RGM$CN_Ratio_control), mean(wide_RGM$CN_Ratio_latrine)))
LIA_CN_mean<-data.frame(soilAge=c('lia','lia'), Type=c('control', 'latrine'), CN_Ratio=c(mean(wide_LIA$CN_Ratio_control), mean(wide_LIA$CN_Ratio_latrine)))

ggplot(nutrient, aes(x=Type, y=CN_Ratio, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=RGM_CN_mean, aes(x = Type, y = CN_Ratio, color=Type), size=5, alpha=.5) +
  geom_point(data=LIA_CN_mean, aes(x = Type, y = CN_Ratio, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=CN_Ratio), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  facet_wrap(~soilAge)+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

#d15N ratio
RGM_d15N_mean<- data.frame(soilAge=c('rgm','rgm'), Type=c('control', 'latrine'), d15N=c(mean(wide_RGM$d15N_control), mean(wide_RGM$d15N_latrine)))
LIA_d15N_mean<-data.frame(soilAge=c('lia','lia'), Type=c('control', 'latrine'), d15N=c(mean(wide_LIA$d15N_control), mean(wide_LIA$d15N_latrine)))

ggplot(nutrient, aes(x=Type, y=d15N, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=RGM_d15N_mean, aes(x = Type, y = d15N, color=Type), size=5, alpha=.5) +
  geom_point(data=LIA_d15N_mean, aes(x = Type, y = d15N, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=d15N), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  facet_wrap(~soilAge)+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

#d13C ratio
RGM_d13C_mean<- data.frame(soilAge=c('rgm','rgm'), Type=c('control', 'latrine'), d13C=c(mean(wide_RGM$d13C_control), mean(wide_RGM$d13C_latrine)))
LIA_d13C_mean<-data.frame(soilAge=c('lia','lia'), Type=c('control', 'latrine'), d13C=c(mean(wide_LIA$d13C_control), mean(wide_LIA$d13C_latrine)))

ggplot(nutrient, aes(x=Type, y=d13C, color=Type, group=Type))+
  geom_point(size=3)+
  geom_point(data=RGM_d13C_mean, aes(x = Type, y = d13C, color=Type), size=5, alpha=.5) +
  geom_point(data=LIA_d13C_mean, aes(x = Type, y = d13C, color=Type), size=5, alpha=.5) +
  stat_summary(aes(x=Type, y=d13C), geom="errorbar", fun.data="mean_se",
               fun.args=list(mult=1), size=1, width=.4) +
  scale_color_manual(values=c('cyan3','purple3'))+ 
  theme_bw()+
  facet_wrap(~soilAge)+
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    axis.text.x=element_text(size=16),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

############ latrines----
#filter for just latrines
nut_lat<- nutrient %>% 
  filter(Type=='latrine') %>% 
  filter(!is.na(mol_C))

#make data wide
wide_lat<-pivot_wider(names_from='soilAge', values_from=c('mass_N','mol_N','Mass_C','mol_C','d15N','d13C','CN_Ratio'), data=nut_lat)

#normality test on the difference between columns
shapiro.test(nut_lat$mol_N)
shapiro.test(nut_lat$mol_C)
shapiro.test(nut_lat$d15N)
shapiro.test(nut_lat$d13C) # not normal
shapiro.test(nut_lat$CN_Ratio)

qqnorm(nut_lat$d13C)

#run paired ttest for control and latrine
t.test(mol_N~soilAge, data=nut_lat)
t.test(mol_C~soilAge, data=nut_lat)
t.test(d15N~soilAge, data=nut_lat)
t.test(CN_Ratio~soilAge, data=nut_lat)
#none of these are different

wilcox.test(d13C~soilAge, data=nut_lat) # not different
