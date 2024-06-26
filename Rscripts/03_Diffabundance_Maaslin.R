###Microbiome infants### 
##Maaslin and plotting for Top families##
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
library(RColorBrewer)
library(ggrepel)

# Settings ----------------------------------------------------------------
theme_set(theme_bw())

# Show rel abundance of most ab families ----------------------------------
data_long<-data_new %>% 
  pivot_longer(Moraxellaceae:Others,names_to = "Family",values_to = "abundance") %>% 
  mutate(ns_month=ns_week/4.35)
data_long$ns_month<-round(data_long$ns_month,digits = 0)
data_long$ns_month

q <-data_long %>% 
  group_by(CF, Family, ns_week) %>% 
  dplyr::summarize(abundance=mean(abundance), .groups = "drop") %>% 
  ggplot(aes(y=abundance,x=ns_week, fill=Family))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~CF)+theme_bw()+
  labs(y="Relative abundance",x="Week after birth")
q

r<-data_long %>% 
  group_by(CF, Family, ns_month) %>% 
  dplyr::summarize(abundance=mean(abundance), .groups = "drop") %>%
  filter(ns_month<13) %>% 
  ggplot(aes(y=abundance,x=ns_month, fill=Family))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~CF)+theme_bw()+
  labs(y="Relative abundance",x="Month after birth")+
  scale_x_continuous(breaks = c(1:12), limits = c(0,13))+
  scale_colour_brewer(palette = "Spectral",aesthetics = "fill")+
theme(legend.position="bottom",legend.title=element_blank())
r

# Show rel abundance of most ab families by symptom era ----------------------------------
data_long_sympt<-data_new_symptoms %>% 
  pivot_longer(Moraxellaceae:Others,names_to = "Family",values_to = "abundance") %>% 
  mutate(ns_month=ns_week/4.35)
data_long_sympt$ns_month<-round(data_long_sympt$ns_month,digits = 0)
data_long_sympt$ns_month
data_long_sympt$era_symptom<-factor(data_long_sympt$era_symptom,labels = c("Before LRTI","At LRTI","After LRTI"))

###some plotting###

q <-data_long_sympt %>% 
  group_by(CF, Family, era_symptom) %>% 
  dplyr::summarize(abundance=mean(abundance), .groups = "drop") %>% 
  ggplot(aes(y=abundance,x=era_symptom, fill=Family))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~CF)+theme_bw()+
  labs(y="Relative abundance",x="First respiratory symptoms")+
  scale_color_brewer(palette = "Spectral",aesthetics = "fill")
q

s<-data_long_sympt %>% 
  group_by(CF, Family, ns_month,era_symptom) %>% 
  dplyr::summarize(abundance=mean(abundance), .groups = "drop") %>%
  filter(ns_month<13) %>% 
  filter(!era_symptom=="At LRTI") %>% 
  ggplot(aes(y=abundance,x=ns_month, fill=Family))+
  geom_bar(position="stack", stat="identity")+
  facet_grid(. ~ era_symptom+CF)+theme_bw()+
  labs(y="Relative abundance",x="Month after birth")+
  scale_x_continuous(limits = c(0,13))+
  scale_colour_brewer(palette = "Spectral",aesthetics = "fill")+
  theme(legend.position = "bottom",legend.title = element_blank())
s
t<-cowplot::plot_grid(r,s,nrow = 2,labels = "AUTO")
t


# -------------------------------------------------------------------------
# Maaslin analysis of most abundant families ------------------------------
Top10<-read_csv("R/Rfiles/20221019top10_families.csv")
df_input_data <-Top10 %>% 
  column_to_rownames("Subject")

df_input_metadata <- data_new[,-c(2:11)] %>% 
  column_to_rownames("Subject")

# Microbiome CF vs Healthy overtime -------------------------------------------------------------
fit_data = Maaslin2(
  input_data = df_input_data, #OTU Table
  input_metadata = df_input_metadata, # Metadata Table
  output = "CFvsHealthy", # Output Folder
  fixed_effects = c("CF", "ns_week", "anf_sibshome", "season", "bf_atswab","ag_sectio_n"),
  random_effects = c("uid"), # patient id
  reference = c("CF,Healthy","season,spring")) #reference if need a certain (e.g. before Trikafta..)

fit_data$results %>% 
  filter(metadata=="CF") %>% 
  ggplot(aes(coef, -log10(qval))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw()+
  labs(col="",x="Coefficient")+
  geom_text_repel(aes(label = feature), data = . %>% filter(qval < 0.05))

sigres<-read_tsv("CFvsHealthy/significant_results.tsv") %>% 
  filter(metadata=="CF")
sigres_df<-as.data.frame(sigres)
sigres_df %>% flextable()%>%
  flextable::save_as_docx(path = "CFvsHealthy/significant_results.docx")


fitdata<-readRDS("CFvsHealthy/fitted.rds")
fit<-data.frame(t(fitdata),check.names = F) %>% 
  rownames_to_column("Subject")
data_fit<-left_join(fit,Meta,by="Subject")


################################################################################
# Maaslin to symptom era --------------------------------------------------
################################################################################

# Compare CF vs Healthy MaAslin before first symptoms ---------------------
df_input_era<-data_new_symptoms %>% 
column_to_rownames("Subject") %>% 
  filter(era_symptom=="Before") ##only 427 samples

df_input_era_sev<-data_new_sevsymptoms %>% 
  column_to_rownames("Subject")%>% 
  filter(era_sevsymptom=="Before") ##only 477 samples

# Before any symptoms -----------------------------------------------------
fit_data = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_era,
  output = "Before_CFvsHealthy",
  fixed_effects = c("CF","ns_week","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring"))
results<-read_tsv("Before_CFvsHealthy/all_results.tsv")

fit_data$results %>% 
  filter(metadata=="CF") %>% 
  ggplot(aes(coef, -log10(qval))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw()+
  labs(col="",x="Coefficient")+
  geom_text_repel(aes(label = feature), data = . %>% filter(qval < 0.05))

#All results
allres_before<-read_tsv("Before_CFvsHealthy/all_results.tsv") %>%
  filter(metadata=="CF")

#Significant results
sigres_before<-read_tsv("Before_CFvsHealthy/significant_results.tsv") %>%
  filter(metadata=="CF")
sigres_df<-as.data.frame(sigres_before)
sigres_df %>% flextable()%>%
  flextable::save_as_docx(path = "Before_CFvsHealthy/significant_results.docx")

# Before sev symptoms occur -----------------------------------------------
fit_data = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_era_sev,
  output = "Beforesev_CFvsHealthy",
  fixed_effects = c("CF","ns_week","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring"))

fit_data$results %>% 
  filter(metadata=="CF") %>% 
  ggplot(aes(coef, -log10(qval))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw()+
  labs(col="",x="Coefficient (Healthy)")+
  geom_text_repel(aes(label = feature), data = . %>% filter(qval < 0.05))

#All results
allres_beforesev<-read_tsv("Beforesev_CFvsHealthy/all_results.tsv") %>%
  filter(metadata=="CF")

Fig2_C<-allres_beforesev %>% 
  filter(metadata=="CF") %>% 
  ggplot(aes(coef, -log10(qval))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw()+
  labs(col="",x="Coefficient (Healthy)")+
  geom_text_repel(aes(label = feature), data = . %>% filter(qval < 1))+
  annotate("segment", x = 0.1, xend = 3, y = 2, yend = 2, arrow = arrow(length = unit(.15, "inches"))) + 
  annotate("label", x = 1.5, y = 2, label = "Higher abundance in Healthy") +
  annotate("segment", x = -0.1, xend = -3, y = 2, yend = 2, arrow = arrow(length = unit(.15, "inches"))) + 
  annotate("label", x = -1.5, y = 2, label = "Higher abundance in CF") +
  annotate("label", x = 3, y = -log10(0.05), label = "q = 0.05")+
  # Some graphical tweaks
  scale_color_nejm(name = "feature") + #library(ggsci)
  labs(x="Effect size",y="P-value, FDR adjusted, -log10",tag = "C")+
  theme_bw() +
  theme(legend.position = "bottom",plot.tag = element_text(face = "bold")) +
  guides(size = "none",
         color = guide_legend(override.aes = list(size = 4)))
Fig2_C

#Significant results
sigres_beforesev<-read_tsv("Beforesev_CFvsHealthy/significant_results.tsv") %>%
  filter(metadata=="CF")
sigres_df<-as.data.frame(sigres_beforesev)
sigres_df %>% flextable()%>%
  flextable::save_as_docx(path = "Beforesev_CFvsHealthy/significant_results.docx")


# Compare CF vs Healthy MaAslin antibiotic naive---------------------
df_input_abnaive<-data_new_abnaive %>% 
  column_to_rownames("Subject")
df_input_abnaive ###1098 samples

# Antibiotic naive CF vs Healthy -----------------------------------------------------
fit_data = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_abnaive,
  output = "ABNaive_CFvsHealthy",
  fixed_effects = c("CF","ns_week","bf_atswab","season","anf_sibshome","ag_sectio_n"),
  random_effects = c("uid"),
  reference = c("season,spring"))

#visualize
fit_data$results %>% 
  filter(metadata=="CF") %>% 
  ggplot(aes(coef, -log10(qval))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw()+
  labs(col="",x="Coefficient")+
  geom_text_repel(aes(label = feature), data = . %>% filter(qval < 0.05))

#All results
allres_naive<-read_tsv("ABNaive_CFvsHealthy/all_results.tsv") %>%
  filter(metadata=="CF")

#Significant results
sigres_naive<-read_tsv("ABNaive_CFvsHealthy/significant_results.tsv") %>%
  filter(metadata=="CF")
sigres_df<-as.data.frame(sigres_naive)
sigres_df %>% flextable()%>%
  flextable::save_as_docx(path = "ABNaive_CFvsHealthy/significant_results.docx")

# After Sev Symptoms ------------------------------------------------------
df_input_after_sev<-data_new_sevsymptoms %>% 
  column_to_rownames("Subject")%>% 
  filter(era_sevsymptom=="After") ##only 594 samples

fit_data = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_after_sev,
  output = "Aftersev_CFvsHealthy",
  fixed_effects = c("CF","ns_week","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring"))

sigres_aftersev<-read_tsv("Aftersev_CFvsHealthy/significant_results.tsv") %>%
  filter(metadata=="CF")
sigres_df<-as.data.frame(sigres_aftersev)
sigres_df %>% flextable()%>%
  flextable::save_as_docx(path = "Aftersev_CFvsHealthy/significant_results.docx")

# After Sev symptoms, before AB -------------------------------------------
df_input_after_sev_naive<-data_new_sevsymptoms %>% 
  filter(Subject%in%data_new_abnaive$Subject) %>% 
  column_to_rownames("Subject")%>% 
  filter(era_sevsymptom=="After") ##only 297 samples

fit_data = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_after_sev_naive,
  output = "Aftersevnaive_CFvsHealthy",
  fixed_effects = c("CF","ns_week","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring"))

sigres_aftersevnaive<-read_tsv("Aftersevnaive_CFvsHealthy/significant_results.tsv") %>%
  filter(metadata=="CF")
sigres_df<-as.data.frame(sigres_aftersevnaive)
sigres_df %>% flextable()%>%
  flextable::save_as_docx(path = "Aftersevnaive_CFvsHealthy/significant_results.docx")

# CF vs healthy before RTI after AB --------------------------------------
df_input_before_sev_afterab<-data_new %>% 
  filter(Subject%in%data_new_beforesev$Subject|sev_frequency=="never")%>% 
  filter(!Subject%in%data_new_abnaive$Subject|CF=="Healthy")%>% 
  filter(ns_week>12) %>% 
  column_to_rownames("Subject") ##only 329 samples

fit_data = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_before_sev_afterab,
  output = "beforesevafterab_CFvsHealthy",
  fixed_effects = c("CF","ns_week","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring"))

sigres_beforesevafterab<-read_tsv("beforesevafterab_CFvsHealthy/significant_results.tsv") %>%
  filter(metadata=="CF")
sigres_df<-as.data.frame(sigres_beforesevafterab)
sigres_df %>% flextable()%>%
  flextable::save_as_docx(path = "beforesevafterab_CFvsHealthy/significant_results.docx")

# Microbiome changes depending on symptom frequency -----------------------
# CF / Symptom frequency --------------------------------------------------
fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_CF, 
  output = "SymptomsCF", 
  fixed_effects = c("symptom_frequency", "ns_week", "ab_atswab","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("symptom_frequency,never","season,spring"))

#Significant results
sigres_CF<-read_tsv("SymptomsCF/significant_results.tsv") %>% 
  filter(metadata=="symptom_frequency")

### In CF: significant changes in microbiome regarding symptom frequency:
#less Carnobacteriaceae

# CF/ Severe Symptom freq --------------------------------------------------------
fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_CF, 
  output = "SevSymptomFreqCF", 
  fixed_effects = c("sev_freq", "ns_week", "ab_atswab","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring"))
#Significant results
sigres_CF<-read_tsv("SevSymptomFreqCF/significant_results.tsv")
###Less neisseriaceae in infants with often RTI
allres_CF<-read_tsv("SevSymptomFreqCF/all_results.tsv") %>% 
  filter(metadata=="sev_freq")

# CF sev_freq in ab naive only -----------------------------------------------
df_input_CF_naive<-data_new_abnaive %>% 
  filter(CF=="CF") %>% 
  column_to_rownames("Subject")

fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_CF_naive, 
  output = "SevSymptomFreqCF_naive", 
  fixed_effects = c("sev_freq", "ns_week","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring"))
#Significant results
sigres_CF<-read_tsv("SevSymptomFreqCF_naive/significant_results.tsv")
###Less neisseriaceae in infants with often LRTI
allres_CF<-read_tsv("SevSymptomFreqCF_naive/all_results.tsv") %>% 
  filter(metadata=="sev_freq")
##less Neiss, less Propioni

# CF Sev Freq before RTI -------------------------------------------------
df_input_CF_beforesev<-data_new_beforesev %>% 
  filter(CF=="CF") %>% 
  column_to_rownames("Subject")

fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_CF_beforesev, 
  output = "SevSymptomFreqCF_beforeLRTI", 
  fixed_effects = c("sev_freq", "ns_week","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring"))
#Significant results
sigres_CF<-read_tsv("SevSymptomFreqCF_beforeLRTI/significant_results.tsv")
###Less Neisseriaceae in infants with often RTI

allres_CF<-read_tsv("SevSymptomFreqCF_beforeLRTI/all_results.tsv") %>% 
  filter(metadata=="sev_freq")

# Maaslin CF no AB and AB separated --------------------------------------------

fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_CF_noab, 
  output = "Results/SymptomEraCFnoAB", 
  fixed_effects = c("era_symptom", "time_firstsymptom", "bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("era_symptom,Before","season,spring"))

#All results
allres_EraCFnoAB<-read_tsv("Results/SymptomEraCFnoAB/all_results.tsv")

fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_CF_ab, 
  output = "Results/SymptomEraCFAB", 
  fixed_effects = c("era_symptom", "time_firstsymptom", "bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("era_symptom,Before","season,spring"))
sigres_EraCFnoAB<-read_tsv("Results/SymptomEraCFnoAB/significant_results.tsv")


#All results
allres_EraCFAB<-read_tsv("Results/SymptomEraCFAB/all_results.tsv") %>% 
  filter(metadata=="")

#Significant results
sigres_EraCFAB<-read_tsv("Results/SymptomEraCFAB/significant_results.tsv")

# Antibiotic era ----------------------------------------------------------
df_input_CF_ab <- data_CF_ab[,-c(1:10)]

fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_CF_ab, 
  output = "ABEraCF", 
  fixed_effects = c("era_ab","ns_week", "bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("era_ab,Before","season,spring"))

#All results
allres_ABera<-read_tsv("ABEraCF/all_results.tsv") %>% 
  filter(metadata=="era_ab")

#Significant results
sigres_ABera<-read_tsv("ABEraCF/significant_results.tsv") %>% 
  filter(metadata=="era_ab")

# CF/ AB freq --------------------------------------------------------
fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_CF, 
  output = "ABFreqCF", 
  fixed_effects = c("ab_freq", "ns_week","bf_atswab","season","anf_sibshome"),
  random_effects = c("uid"),
  reference = c("season,spring","ab_freq,never"))

#Significant results
sigres_CF<-read_tsv("ABFreqCF/significant_results.tsv")%>% 
  filter(metadata=="ab_freq")

allres_CF<-read_tsv("ABFreqCF/all_results.tsv") %>% 
  filter(metadata=="ab_freq")
