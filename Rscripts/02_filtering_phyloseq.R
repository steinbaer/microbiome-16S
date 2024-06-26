# After Using decontam, this script was used for filtering --------

library(BiocManager)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(readr)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(metagMisc)


sam_data_new<-read_csv("metadata.csv")
ps_decontam<-readRDS("R/R_Decontam/R_Decontam/phyloseq_decontam.rds")

#remove controls
ps_decontam_noneg = subset_samples(ps_decontam, Sample.name != "NTC-P1" & Sample.name != "NTC-P12" & Sample.name != "NTC-P13"& Sample.name != "NTC-P14"& 
                                     Sample.name != "NTC-P15"& Sample.name != "NTC-P19"& Sample.name != "NTC-P2"& Sample.name != "NTC-P3"& Sample.name != "NTC-P4"
                                   & Sample.name != "NTC-P6"& Sample.name != "NTC-P7"& Sample.name != "NTC-P8"& Sample.name != "NTC-P9")

##slice sam_data_new
sam_data_new<-slice(sam_data_new, -(1558:1570))

#row.names
row.names(sam_data_new)<-sam_data_new$Sample.name

#as phyloseq
sam_data_new <- sample_data(sam_data_new)

#merge phyloseq
phyloseq<-merge_phyloseq(ps_decontam_noneg,sam_data_new)

#Remove badTaxa (many reads in ntc but not picked up by decontam, == Mitochondria)
badTaxa = c("ASV14")
allTaxa = taxa_names(phyloseq)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps_ex1 = prune_taxa(allTaxa, phyloseq)
ps_ex1_nosingle <- prune_taxa(taxa_sums(ps_ex1) > 1, ps_ex1)

#In addition exclude low read numbers
ps_ex1_lowreads <- prune_samples(sample_sums(ps_ex1_nosingle)>=3000, ps_ex1_nosingle)

# Remove Class chloroplast -------------------------------------------------
ps<-ps_ex1 %>% subset_taxa(Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class)) 
ps_ex1<- metagMisc::phyloseq_filter_prevalence(ps, prev.trh=0.001, abund.trh=5, threshold_condition="OR") 