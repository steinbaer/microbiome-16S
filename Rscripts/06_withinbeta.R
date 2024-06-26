###Within beta div####
#Beta-Diversities within individual (pad)
#Author: Ruth Steinberg

library(phyloseq)
library(Biostrings)
library(ggplot2)
library(readr)
library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)
library(tidyverse)
library(mgcv)
library(Maaslin2)
library(vegan)
library(plyr)
library(lme4)

# Load dist matrix and Metadata--------------------------------------------------------
load("R/Rfiles/20221026alldata.RData")
load("R/Rfiles/distance.RData")

# All individuals ---------------------------------------------------------
bray_dist_all<-as.matrix(dist_Bray) %>% 
  as_tibble(rownames = "sample") %>% #tibble
  pivot_longer(-sample) %>% #long
  filter(sample<name) %>% 
  separate(sample,into = c("ID_a","number_a"), sep = "-", remove = FALSE) %>% #get IDs
  separate(name,into = c("ID_b","number_b"), sep = "-", remove = FALSE) %>% 
  select(-number_a, -number_b)#get IDs

meta_all<-data_new%>% 
  select(id,CF,ns_week,symptom_freq, sev_freq, ab_freq, breastfeeding,bf_atswab, season, ag_sectio_n,anf_sibshome,ab_atswab)

meta_all_2<-data_new %>% 
  select(id,CF,ns_week,symptom_freq, sev_freq, ab_freq)

bray_dist_all_time<-left_join(bray_dist_all,meta_all,by=c("sample"="id"))%>%
  dplyr::rename(CF_a=CF) %>% 
  dplyr::rename(ns_week_a=ns_week) %>% 
  dplyr::rename(symptom_freq_a=symptom_freq) %>% 
  dplyr::rename(sev_freq_a=sev_freq) %>% 
  dplyr::rename(ab_freq_a=ab_freq)

bray_dist_all_time<-left_join(bray_dist_all_time,meta_all_2,by=c("name"="id")) %>%
  dplyr::rename(CF_b=CF) %>% 
  dplyr::rename(ns_week_b=ns_week) %>% 
  dplyr::rename(symptom_freq_b=symptom_freq) %>% 
  dplyr::rename(sev_freq_b=sev_freq) %>% 
  dplyr::rename(ab_freq_b=ab_freq)

###Filter
bray_dist_within_all<-bray_dist_all_time %>% 
  mutate(diff=abs(ns_week_a-ns_week_b)) %>% 
  filter(ID_a==ID_b,diff<3) 
table(bray_dist_within_all$CF_a)
# CF Healthy 
# 707     442 

bray_dist_within_naive<-bray_dist_within_all %>% 
  filter(sample%in%data_new_abnaive$Subject) %>% 
  filter(name%in%data_new_abnaive$Subject)

bray_dist_naive_CF<-bray_dist_within_naive %>% 
  filter(CF_a=="CF")

table(bray_dist_within_naive$CF_a)
# CF Healthy 
# 397     426 
bray_dist_within_naive_tab<-bray_dist_within_naive %>% select(CF_a,ID_a) %>% distinct()
table(bray_dist_within_naive_tab$CF_a)
# CF Healthy 
# 42      30


# Calculate gamm and lmer for beta within ---------------------------------
gamm.within_all<-gamm(value ~ CF_a +
                        s(ns_week_a, bs="cr",fx=FALSE)+
                        as.factor(season)+
                        as.factor(bf_atswab)+
                        as.factor(ab_atswab)+
                        as.factor(anf_sibshome),
                      family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_dist_within_all,na.action=na.omit)
summary(gamm.within_all$gam)

bray_dist_within_all$fitBeta<-predict(gamm.within_all$gam)
intervals(gamm.within_all$lme)

lmer.within_all<-lmer(median~CF_a+bf_atswab+season+anf_sibshome+(ns_week_a|ID_a),data = bray_dist_within_all)
summary(lmer.within_all)
ggplot(bray_dist_within_all,aes(x=ns_week_a,y=median,colour=CF_a))+geom_point(size=1,alpha=0.5)+geom_smooth(method = "gam")+theme_bw()

p<-bray_dist_within_all %>% 
  mutate(control=ifelse(CF_a=="Healthy","Controls","CF")) %>% 
  ggplot(aes(x=ns_week_a,y=fitBeta,colour=control))+
  geom_point(alpha=0.3)+
  theme_bw()+geom_smooth()+
  labs(x="Age (weeks)",y="Within-subject dissimilarity",col="")+
  ylim(0.48,0.72)
p  


# Subset to CF and healthy ------------------------------------------------
bray_dist_within_CF<-bray_dist_within_all %>% 
  filter(CF_a=="CF")
summary(bray_dist_within_CF$sev_freq_a)
# normal  often 
# 464    243
tab<-bray_dist_within_CF %>% select(ID_a,CF_a,sev_freq_a) %>% 
  distinct()
table(tab$CF_a,tab$sev_freq_a)

# CF and severe frequency -------------------------------------------------
modCF<-lmer(median~sev_freq_a+(1|ID_a),data = bray_dist_within_CF)
summary(modCF)

gamm.mod_betaCF_sev<-gamm(value ~ sev_freq_a+
                          s(ns_week_a, bs="cr",fx=FALSE),
                          as.factor(season)+
                          as.factor(ab_atswab)+
                          as.factor(bf_atswab)+
                          as.factor(anf_sibshome),
                          family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_dist_within_CF,na.action=na.omit)
summary(gamm.mod_betaCF_sev$gam)

intervals(gamm.mod_betaCF_sev$lme)
bray_dist_within_CF$fitBeta<-predict(gamm.mod_betaCF_sev$gam)
a<-ggplot(bray_dist_within_CF,aes(y=fitBeta,x=ns_week_a,col=sev_freq_a))+
  geom_point(size=1,alpha=0.5)+
  theme_bw()+geom_smooth()+
  labs(x="Age (weeks)",y="Consecutive Distance",col="")+
  ylim(0.5,0.75)
  
##### What about AB Naive CF?
table(bray_dist_naive_CF$sev_freq_a)
# normal  often 
# 282    115 
tab<-bray_dist_naive_CF %>% select(sev_freq_a,CF_a,ID_a) %>% 
  distinct()
table(tab$CF_a,tab$sev_freq_a)
gamm.mod_betaCF_sevnaive<-gamm(value ~ sev_freq_a+
                            s(ns_week_a, bs="cr",fx=FALSE)+
                            as.factor(season)+
                            as.factor(bf_atswab)+
                            as.factor(anf_sibshome),
                          family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_dist_naive_CF,na.action=na.omit)
summary(gamm.mod_betaCF_sev$gam)
  
intervals(gamm.mod_betaCF_sevnaive$lme,which = "fixed")

bray_dist_naive_CF$fitBeta<-predict(gamm.mod_betaCF_sevnaive$gam)
b<-ggplot(bray_dist_naive_CF,aes(y=fitBeta,x=ns_week_a,col=sev_freq_a))+
  geom_point(size=1,alpha=0.5)+
  theme_bw()+geom_smooth()+
  labs(x="Age (weeks)",y="Consecutive Distance",col="")+
  ylim(0.5,0.75)
b

#### What about before RTI?
data_beforesev_CF<-data_new_beforesev %>% filter(CF=="CF")
bray_dist_CF_before<-bray_dist_within_all %>% 
  filter(sample%in%data_beforesev_CF$Subject)
tab<-bray_dist_CF_before %>% select(sev_freq_a,CF_a,ID_a) %>% 
  distinct()
table(tab$CF_a,tab$sev_freq_a)
# often 14 normal 33
table(bray_dist_CF_before$sev_freq_a)

gamm.mod_betaCF_befsev<-gamm(value ~ sev_freq_a+
                                 s(ns_week_a, bs="cr",fx=FALSE)+
                                 as.factor(season)+
                                 as.factor(bf_atswab)+
                                 as.factor(ab_atswab)+
                                 as.factor(anf_sibshome),
                               family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_dist_CF_before,na.action=na.omit)
summary(gamm.mod_betaCF_befsev$gam) 

intervals(gamm.mod_betaCF_befsev$lme)
bray_dist_CF_before$fitBeta<-predict(gamm.mod_betaCF_befsev$gam)
bray_dist_CF_before %>% filter(sev_freq_a=="normal")%>% mutate(mean=mean(value)) %>% 
  select(mean) %>% distinct()

c<-ggplot(bray_dist_CF_before,aes(x=ns_week_a,y=value,col=sev_freq_a))+
  geom_point(size=1,alpha=0.5)+
  theme_bw()+geom_smooth()+
  labs(x="Age (weeks)",y="Consecutive Distance",col="")
c
cowplot::plot_grid(a,b,nrow = 2,labels = "AUTO")
ggsave("R/Results/Distance_ordination/withinbeta_sev_freq.png",height = 10,width = 12)

# Antibiotic naive individuals ---------------------------------------------------------
# Calculate gamm and lmer for beta within ab_naive ---------------------------------
gamm.within_naive<-gamm(value ~ CF_a +
                        s(ns_week_a, bs="cr",fx=FALSE)+
                        as.factor(season)+
                        as.factor(bf_atswab)+
                        as.factor(ag_sectio_n)+
                        as.factor(anf_sibshome),
                      family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_dist_within_naive,na.action=na.omit)
summary(gamm.within_naive$gam)

bray_dist_within_naive$fitBeta<-predict(gamm.within_naive$gam)
intervals(gamm.within_naive$lme,which="fixed")

lmer.within_naive<-lmer(median~CF_a+bf_atswab+season+anf_sibshome+(ns_week_a|ID_a),data = bray_dist_within_naive)
summary(lmer.within_naive)
ggplot(bray_dist_within_naive,aes(x=ns_week_a,y=median,colour=CF_a))+geom_point(size=1,alpha=0.5)+geom_smooth(method = "gam")+theme_bw()
q<-bray_dist_within_naive %>% 
  mutate(control=ifelse(CF_a=="Healthy","Controls","CF")) %>% 
  ggplot(aes(x=ns_week_a,y=fitBeta,colour=control))+
  geom_point(size=1,alpha=0.5)+
  theme_bw()+geom_smooth()+
  labs(x="Age (weeks)",y="Within-subject dissimilarity",col="")+
  ylim(0.48,0.72)
q
r<-cowplot::plot_grid(p,q,nrow = 2,labels="auto")
r

# Antibiotic individuals before RTI (dis_abnosev) ---------------------------------------------------------
bray_abnosev<-as.matrix(dis_abnosev) %>% 
  as_tibble(rownames = "sample") %>% #tibble
  pivot_longer(-sample) %>% #long
  filter(sample<name) %>% 
  separate(sample,into = c("ID_a","number_a"), sep = "-", remove = FALSE) %>% #get IDs
  separate(name,into = c("ID_b","number_b"), sep = "-", remove = FALSE) %>% 
  select(-number_a, -number_b) #get IDs

meta_abnosev_within<-meta_abnosev %>% 
  select(id,CF,ns_week,symptom_freq, sev_freq, breastfeeding,bf_atswab, season, ag_sectio_n,anf_sibshome,ab_atswab)

meta_abnosev_2<-meta_abnosev %>% 
  select(id,CF,ns_week,symptom_freq, sev_freq)

bray_dist_abnosev<-right_join(bray_abnosev,meta_abnosev_within,by=c("sample"="id"))%>%
  dplyr::rename(CF_a=CF) %>% 
  dplyr::rename(ns_week_a=ns_week) %>% 
  dplyr::rename(symptom_freq_a=symptom_freq) %>% 
  dplyr::rename(sev_freq_a=sev_freq) 

bray_dist_abnosev<-left_join(bray_dist_abnosev,meta_abnosev_2,by=c("name"="id")) %>%
  dplyr::rename(CF_b=CF) %>% 
  dplyr::rename(ns_week_b=ns_week) %>% 
  dplyr::rename(symptom_freq_b=symptom_freq) %>% 
  dplyr::rename(sev_freq_b=sev_freq)

###Filter
bray_dist_within_abnosev<-bray_dist_abnosev %>% 
  mutate(diff=abs(ns_week_a-ns_week_b)) %>% 
  filter(ID_a==ID_b,diff<3) %>%
  group_by(ID_a,ns_week_a,CF_a,sev_freq_a,diff)%>% 
  dplyr::mutate(median=median(value)) %>% 
  ungroup()


# Calculate gamm and lmer for beta within ab_nosev ---------------------------------
gamm.within_abnosev<-gamm(median ~ CF_a +
                          s(ns_week_a, bs="cr",fx=FALSE)+
                          as.factor(season)+
                          as.factor(bf_atswab)+
                          as.factor(ag_sectio_n)+
                          as.factor(anf_sibshome),
                        family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_dist_within_abnosev,na.action=na.omit)
summary(gamm.within_abnosev$gam)

# before RTI (dis_beforesev) ---------------------------------------------------------
ps_beforesev<-subset_samples(ps, id%in%data_new_beforesev$Subject)
ps_beforesev<-transform_sample_counts(ps_beforesev,function (x) x/sum(x))
dis_beforesev<-phyloseq::distance(ps_beforesev,method="bray")

bray_beforesev<-as.matrix(dis_beforesev) %>% 
  as_tibble(rownames = "sample") %>% #tibble
  pivot_longer(-sample) %>% #long
  filter(sample<name) %>% 
  separate(sample,into = c("ID_a","number_a"), sep = "-", remove = FALSE) %>% #get IDs
  separate(name,into = c("ID_b","number_b"), sep = "-", remove = FALSE) %>% 
  select(-number_a, -number_b) #get IDs

meta_beforesev<-get_variable(ps_beforesev)
meta_beforesev_within<-meta_beforesev %>% 
  select(id,CF,ns_week,symptom_freq, sev_freq, breastfeeding,bf_atswab, season, ag_sectio_n,anf_sibshome,ab_atswab)

meta_beforesev_2<-meta_beforesev %>% 
  select(id,CF,ns_week,symptom_freq, sev_freq)

bray_dist_beforesev<-right_join(bray_beforesev,meta_beforesev_within,by=c("sample"="id"))%>%
  dplyr::rename(CF_a=CF) %>% 
  dplyr::rename(ns_week_a=ns_week) %>% 
  dplyr::rename(symptom_freq_a=symptom_freq) %>% 
  dplyr::rename(sev_freq_a=sev_freq) 

bray_dist_beforesev<-left_join(bray_dist_beforesev,meta_beforesev_2,by=c("name"="id")) %>%
  dplyr::rename(CF_b=CF) %>% 
  dplyr::rename(ns_week_b=ns_week) %>% 
  dplyr::rename(symptom_freq_b=symptom_freq) %>% 
  dplyr::rename(sev_freq_b=sev_freq)

###Filter
bray_dist_within_beforesev<-bray_dist_beforesev %>% 
  mutate(diff=abs(ns_week_a-ns_week_b)) %>% 
  filter(ID_a==ID_b,diff<3) %>%
  group_by(ID_a,ns_week_a,CF_a,sev_freq_a,diff)%>% 
  dplyr::mutate(median=median(value)) %>% 
  ungroup()
bray_dist_within_beforesev_tab<-bray_dist_within_beforesev %>% select(ID_a,CF_a) %>% distinct()
table(bray_dist_within_beforesev$CF_a)
# CF Healthy 
# 305     263 
table(bray_dist_within_beforesev_tab$CF_a)
# CF Healthy 
# 43      27

bray_dist_within_all_tab<-bray_dist_within_all %>% select(ID_a,CF_a) %>% distinct()
table(bray_dist_within_all_tab$CF_a)
table(bray_dist_within_all$CF_a)

# Calculate gamm and lmer for beta within ab_nosev ---------------------------------
gamm.within_beforesev<-gamm(median ~ CF_a +
                            s(ns_week_a, bs="cr",fx=FALSE)+
                            as.factor(season)+
                            as.factor(bf_atswab)+
                            as.factor(ag_sectio_n)+
                            as.factor(anf_sibshome),
                          family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_dist_within_beforesev,na.action=na.omit)
summary(gamm.within_beforesev$gam)

intervals(gamm.within_beforesev$lme,which = "fixed")

# after RTI (dis_noabsev) ---------------------------------------------------------

bray_aftersev<-as.matrix(dis) %>% 
  as_tibble(rownames = "sample") %>% #tibble
  pivot_longer(-sample) %>% #long
  filter(sample<name) %>% 
  separate(sample,into = c("ID_a","number_a"), sep = "-", remove = FALSE) %>% #get IDs
  separate(name,into = c("ID_b","number_b"), sep = "-", remove = FALSE) %>% 
  select(-number_a, -number_b) #get IDs

meta_aftersev_within<-data_new_sevsymptoms %>% 
filter(era_sevsymptom=="After") %>% 
  select(id,CF,ns_week,symptom_freq, sev_freq, breastfeeding,bf_atswab, season, ag_sectio_n,anf_sibshome,ab_atswab)

meta_aftersev_2<-data_new_sevsymptoms %>% 
  filter(era_sevsymptom=="After") %>% 
  select(id,CF,ns_week,symptom_freq, sev_freq)

bray_dist_aftersev<-right_join(bray_aftersev,meta_aftersev_within,by=c("sample"="id"))%>%
  dplyr::rename(CF_a=CF) %>% 
  dplyr::rename(ns_week_a=ns_week) %>% 
  dplyr::rename(symptom_freq_a=symptom_freq) %>% 
  dplyr::rename(sev_freq_a=sev_freq) 

bray_dist_aftersev<-left_join(bray_dist_aftersev,meta_aftersev_2,by=c("name"="id")) %>%
  dplyr::rename(CF_b=CF) %>% 
  dplyr::rename(ns_week_b=ns_week) %>% 
  dplyr::rename(symptom_freq_b=symptom_freq) %>% 
  dplyr::rename(sev_freq_b=sev_freq)

###Filter
bray_dist_within_aftersev<-bray_dist_aftersev %>% 
  mutate(diff=abs(ns_week_a-ns_week_b)) %>% 
  filter(ID_a==ID_b,diff<3) %>%
  group_by(ID_a,ns_week_a,CF_a,sev_freq_a,diff)%>% 
  dplyr::mutate(median=median(value)) %>% 
  ungroup()
bray_dist_within_aftersev_tab<-bray_dist_within_aftersev %>% select(ID_a,CF_a) %>% distinct()
table(bray_dist_within_aftersev$CF_a)
# CF Healthy 
# 310     134
table(bray_dist_within_aftersev_tab$CF_a)
# CF Healthy 
# 38      20


# Calculate gamm and lmer for beta within after RTI ---------------------------------
gamm.within_aftersev<-gamm(median ~ CF_a +
                              s(ns_week_a, bs="cr",fx=FALSE)+
                              as.factor(season)+
                              as.factor(bf_atswab)+
                              as.factor(ag_sectio_n)+
                              as.factor(anf_sibshome),
                            family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_dist_within_aftersev,na.action=na.omit)
summary(gamm.within_aftersev$gam)
intervals(gamm.within_aftersev$lme,which = "fixed")

# Subset to CF after AB and Healthy ---------------------------------------
data_CF_ab_after<-data_CF_ab %>% 
  filter(era_ab=="After")
bray_dist_within_CF_ab<-bray_dist_within_all %>% 
  filter(sample%in%data_CF_ab_after$Subject)
bray_dist_within_healthy<-bray_dist_within_all %>% 
  filter(CF_a=="Healthy")

bray_within_healthy_CFab<-rbind(bray_dist_within_CF_ab,bray_dist_within_healthy)
table(bray_within_healthy_CFab$CF_a)
# CF Healthy 
# 287     442
bray_within_healthy_CFab_tab<-bray_within_healthy_CFab %>% select(CF_a,ID_a) %>% distinct()
table(bray_within_healthy_CFab_tab$CF_a)
# CF Healthy 
# 30      30 
gamm.within_CFab_healthy<-gamm(median ~ CF_a +
                                 s(ns_week_a, bs="cr",fx=FALSE)+
                                 as.factor(season)+
                                 as.factor(bf_atswab)+
                                 as.factor(ab_atswab)+
                                 as.factor(anf_sibshome),
                               family=gaussian,correlation=corAR1(form=~1|ID_a),data=bray_within_healthy_CFab,na.action=na.omit)
summary(gamm.within_CFab_healthy$gam)
intervals(gamm.within_CFab_healthy$lme)
