#Decontam#
#Vignette: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(lubridate)
library(readxl)
library(phyloseq)
library(decontam)

# read DADA2 output -------------------------------------------------------
ps<-readRDS("DADA2.rds")
# Prevalence of ASVs in neg controls --------------------------------------

sample_data(ps)$is.neg <- sample_data(ps)$sample_control == "Control"
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")


table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_control == "Sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
ps.noncontam
