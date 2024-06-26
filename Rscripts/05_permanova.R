# PERMANOVA Testing, beta diversity ---------------------------------------------------------
# Author: Ruth Steinberg
# Last change:08.03.2023

library(BiocManager)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(vegan)
library(ggforce)

# Load data  ---------------------------------------------------------------
load("R/Rfiles/20221026alldata.RData")
 
# CF vs Healthy -----------------------------------------------------------

ps_RA<-transform_sample_counts(ps,function (x) x/sum(x))
dis_all<-phyloseq::distance(ps_RA,method = "bray")
meta_all<-get_variable(ps_RA)

meta_all$CF<-as.factor(meta_all$CF)
meta_all$ab_atswab<-as.factor(meta_all$ab_atswab)
meta_all$bf_atswab<-as.factor(meta_all$bf_atswab)
meta_all$anf_sibshome<-as.factor(meta_all$anf_sibshome)
meta_all$season<-as.factor(meta_all$season)
meta_all$ag_sectio_n<-as.factor(meta_all$ag_sectio_n)

set.seed(13081992)
#Adonis2
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_all,uid)
adonis2_res<-adonis2(dis_all~CF+ab_atswab+bf_atswab+season+anf_sibshome,by="margin", permutations = perm, data=meta_all) #Terms added marginally
adonis2_res

# CF vs Healthy before sev symptoms, ab naive -----------------------------
ps_beforenaive<-subset_samples(ps, id%in%data_new_naivebeforesev$Subject)
ps_beforenaive<-transform_sample_counts(ps_beforenaive,function (x) x/sum(x))
dis_beforenaive<-phyloseq::distance(ps_beforenaive,method="bray")

meta_beforenaive<-get_variable(ps_beforenaive)
meta_beforenaive$CF<-as.factor(meta_beforenaive$CF)
meta_beforenaive$ab_atswab<-as.factor(meta_beforenaive$ab_atswab)
meta_beforenaive$bf_atswab<-as.factor(meta_beforenaive$bf_atswab)
meta_beforenaive$anf_sibshome<-as.factor(meta_beforenaive$anf_sibshome)
meta_beforenaive$season<-as.factor(meta_beforenaive$season)
meta_beforenaive$ag_sectio_n<-as.factor(meta_beforenaive$ag_sectio_n)

###PERMANOVA
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_beforenaive,uid)
adonis2(dis_beforenaive~CF+bf_atswab+season+anf_sibshome+ag_sectio_n+ns_week,by="margin", permutations = perm, data=meta_beforenaive)#

# Before first antibiotics -------------------------
ps_naive<-subset_samples(ps, id%in%data_new_abnaive$Subject)
ps_naive<-transform_sample_counts(ps_naive,function (x) x/sum(x))
dis_naive<-phyloseq::distance(ps_naive,method="bray")

meta_naive<-get_variable(ps_naive)
meta_naive$CF<-as.factor(meta_naive$CF)
meta_naive$ab_atswab<-as.factor(meta_naive$ab_atswab)
meta_naive$bf_atswab<-as.factor(meta_naive$bf_atswab)
meta_naive$anf_sibshome<-as.factor(meta_naive$anf_sibshome)
meta_naive$season<-as.factor(meta_naive$season)
meta_naive$ag_sectio_n<-as.factor(meta_naive$ag_sectio_n)

table(meta_naive$CF) #539 CF, 559 Healthy
meta_naive_tab<-meta_naive %>% select(CF,uid) %>% distinct()
table(meta_naive_tab$CF) #46 CF; 30 Healthy

###PERMANOVA
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_naive,uid)
adonis2(dis_naive~CF+bf_atswab+season+anf_sibshome+ag_sectio_n+ns_week,by="margin", permutations = perm, data=meta_naive)#gives adonis2 same output?

# Before first RTI -------------------------
ps_beforesev<-subset_samples(ps, id%in%data_new_beforesev$Subject)
ps_beforesev<-transform_sample_counts(ps_beforesev,function (x) x/sum(x))
dis_beforesev<-phyloseq::distance(ps_beforesev,method="bray")

meta_beforesev<-get_variable(ps_beforesev)
meta_beforesev$CF<-as.factor(meta_beforesev$CF)
meta_beforesev$ab_atswab<-as.factor(meta_beforesev$ab_atswab)
meta_beforesev$bf_atswab<-as.factor(meta_beforesev$bf_atswab)
meta_beforesev$anf_sibshome<-as.factor(meta_beforesev$anf_sibshome)
meta_beforesev$season<-as.factor(meta_beforesev$season)
meta_beforesev$ag_sectio_n<-as.factor(meta_beforesev$ag_sectio_n)

table(meta_beforesev$CF) #435 CF, 357 Healthy
meta_beforesev_tab<-meta_beforesev %>% select(CF,uid) %>% distinct()
table(meta_beforesev_tab$CF) #48 CF; 30 Healthy

###PERMANOVA
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_beforesev,uid)
adonis2(dis_beforesev~CF+bf_atswab+season+anf_sibshome+ag_sectio_n+ns_week,by="margin", permutations = perm, data=meta_beforesev)#gives adonis2 same output?

# Subset to CF after AB and before LRTI -----------------------------------
data_new_beforesev_healthy<-data_new_beforesev %>% 
  filter(CF=="Healthy")##beforesev or never ab

data_CF_abbeforesev<-data_CF_ab %>%
  filter(era_ab=="After") %>% 
  filter(Subject%in%data_new_beforesev$Subject) # CF beforeSev,afterAB

ps_abnosev<-subset_samples(ps,id%in%data_CF_abbeforesev$Subject|id%in%data_new_beforesev_healthy$Subject)
ps_abnosev<-transform_sample_counts(ps_abnosev,function (x) x/sum(x))
dis_abnosev<-phyloseq::distance(ps_abnosev,method="bray")

meta_abnosev<-get_variable(ps_abnosev)
meta_abnosev$CF<-as.factor(meta_abnosev$CF)
meta_abnosev$ab_atswab<-as.factor(meta_abnosev$ab_atswab)
meta_abnosev$bf_atswab<-as.factor(meta_abnosev$bf_atswab)
meta_abnosev$anf_sibshome<-as.factor(meta_abnosev$anf_sibshome)
meta_abnosev$season<-as.factor(meta_abnosev$season)
meta_abnosev$ag_sectio_n<-as.factor(meta_abnosev$ag_sectio_n)

table(meta_abnosev$CF)
# CF Healthy 
# 75     357
meta_abnosev_tab<-meta_abnosev %>% select(uid,CF) %>% distinct()
table(meta_abnosev_tab$CF) #13 CF vs 30 Healthy

###PERMANOVA
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_abnosev,uid)
adonis2(dis_abnosev~CF+bf_atswab+season+anf_sibshome+ag_sectio_n+ns_week,by="margin", permutations = perm, data=meta_abnosev)

# Subset after RTI and before AB -----------------------------------
ps_noabsev<-subset_samples(ps, id%in%data_new_aftersevnaive$Subject)
ps_noabsev<-transform_sample_counts(ps_noabsev,function (x) x/sum(x))
dis_noabsev<-phyloseq::distance(ps_noabsev,method="bray")

meta_noabsev<-get_variable(ps_noabsev)
meta_noabsev$CF<-as.factor(meta_noabsev$CF)
meta_noabsev$ab_atswab<-as.factor(meta_noabsev$ab_atswab)
meta_noabsev$bf_atswab<-as.factor(meta_noabsev$bf_atswab)
meta_noabsev$anf_sibshome<-as.factor(meta_noabsev$anf_sibshome)
meta_noabsev$season<-as.factor(meta_noabsev$season)
meta_noabsev$ag_sectio_n<-as.factor(meta_noabsev$ag_sectio_n)

table(meta_noabsev$CF,meta_noabsev$ns_week)
table(meta_noabsev$CF)
# CF Healthy 
# 135     162 
meta_noabsev_tab<-meta_noabsev %>% select(uid,CF) %>% distinct()
table(meta_noabsev_tab$CF) #20 CF vs 19 Healthy

###PERMANOVA
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_noabsev,uid)
adonis2(dis_noabsev~CF+bf_atswab+season+anf_sibshome+ag_sectio_n+ns_week,by="margin", permutations = perm, data=meta_noabsev)

# Look before symptoms vs after symptoms ----------------------------------
#I) Before
ps_before<-subset_samples(ps_sympt,era_symptom=="Before") #441 in subset
ps_before<-transform_sample_counts(ps_before, function (x) x/sum(x))
dis_before<-phyloseq::distance(ps_before,method = "bray") #get Bray-distance

meta_before<-get_variable(ps_before) 
meta_before$era_symptom <-as.factor(meta_before$era_symptom)
meta_before$ab_atswab<-as.factor(meta_before$ab_atswab)
meta_before$bf_atswab<-as.factor(meta_before$bf_atswab)
meta_before$anf_sibshome<-as.factor(meta_before$anf_sibshome)
meta_before$season<-as.factor(meta_before$season)

#Calculate adonis
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_before,uid)
adonis2_before<-adonis2(dis_before~CF+bf_atswab+season+anf_sibshome,by="margin", permutations = perm, data=meta_before)
adonis2_before

###II) After
ps_after<-subset_samples(ps_2,era_symptom=="After") #
dis_after<-phyloseq::distance(ps_after,method = "bray") #get Bray-distance

meta_after<-get_variable(ps_after) 
meta_after$era_symptom <-as.factor(meta_after$era_symptom)
meta_after$ab_atswab<-as.factor(meta_after$ab_atswab)
meta_after$bf_atswab<-as.factor(meta_after$bf_atswab)
meta_after$anf_sibshome<-as.factor(meta_after$anf_sibshome)
meta_after$season<-as.factor(meta_after$season)

#Calculate adonis
adonis2_after<-adonis2(dis_after~CF+bf_atswab+season+anf_sibshome,by="margin", permutations = perm, data=meta_after)
adonis2_after

# CF and RTI frequency--------------------------------------------------------
ps_CF<-subset_samples(ps, CF=="CF")
ps_CF #933 samples
ps_CF<-transform_sample_counts(ps_CF,function (x) x/sum(x))
dis_CF<-phyloseq::distance(ps_CF,method = "bray") #get Bray-distance

meta_CF<-get_variable(ps_CF)
meta_CF$ab_atswab<-as.factor(meta_CF$ab_atswab)
meta_CF$bf_atswab<-as.factor(meta_CF$bf_atswab)
meta_CF$anf_sibshome<-as.factor(meta_CF$anf_sibshome)
meta_CF$season<-as.factor(meta_CF$season)
meta_CF$ag_sectio_n<-as.factor(meta_CF$ag_sectio_n)
meta_CF$sev_freq<-as.factor(meta_CF$sev_freq)
meta_CF$ab_freq<-as.factor(meta_CF$ab_freq)

#Adonis2
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_CF,uid)
adonis2_res_CF<-adonis2(dis_CF~sev_freq+ab_atswab+bf_atswab+season+anf_sibshome+ag_sectio_n+ns_week,by="margin", permutations = perm, data=meta_CF)#gives adonis2 same output?
adonis2_res_CF

# Antibiotics in CF ------------------------------------------
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_CF,uid)
adonis2_abfreq<-adonis2(dis_CF~ab_freq+bf_atswab+season+anf_sibshome+ag_sectio_n+ns_week,by="margin", permutations = perm, data=meta_CF)#gives adonis2 same output?
adonis2_abfreq

ps_ab<-subset_samples(ps,id%in%data_CF_ab$Subject)
ps_ab<-transform_sample_counts(ps_ab,function (x) x/sum(x))
data_CF_ab_ps<-data_CF_ab%>% 
  column_to_rownames("Subject") %>% 
  sample_data() 
ps_ab<-merge_phyloseq(ps_ab,data_CF_ab_ps)
meta_CF_ab<-get_variable(ps_ab) #635

dis_CF_ab<-phyloseq::distance(ps_ab, method = "bray") #get Bray-distance
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_CF_ab,uid)
adonis2_abera<-adonis2(dis_CF_ab~era_ab+bf_atswab+season+anf_sibshome+ag_sectio_n+ns_week,by="margin", permutations = perm, data=meta_CF_ab)#gives adonis2 same output?
adonis2_abera

# Now only CF with severe symptoms ---------------------------------------------
ps_CF_sevsympt<-subset_samples(ps_sevsympt,CF=="CF")
ps_CF_sevsympt<-transform_sample_counts(ps_CF_sevsympt,function (x) x/sum(x))
meta_CF_sevsympt<-get_variable(ps_CF_sevsympt)
dis_CF_sevsympt<-phyloseq::distance(ps_CF_sevsympt,method = "bray") #get Bray-distance

###Adonis2 (marginal)
perm<-how(nperm=999)
setBlocks(perm)<-with(meta_CF_sevsympt,uid)
adonis2_sev_CF<-adonis2(dis_CF_sevsympt ~ era_sevsymptom+ab_atswab+sev_freq+bf_atswab+season+anf_sibshome+ag_sectio_n,by="margin", permutations = perm, data=meta_CF_sevsympt)#gives adonis2 same output?
adonis2_sev_CF
