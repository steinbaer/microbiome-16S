#!/usr/bin/env Rscript

##############
#Load librarys

library(optparse)
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(gridExtra)
library(ShortRead)

sessionInfo()

################
#Specify options
option_list = list(
  make_option(c("-e", "--extention"), type="character", default="fastq",
              help="input file name extention [default='fastq']"),
  make_option(c("-n", "--name"), type="integer", default=3,
              help="In which part of the name is the sample ID (separated by '_') [default=3]"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input folder with 16S amplicon reads."),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Path to output folder."),
  make_option(c("-d", "--db"), type="character", default=NULL,
              help="Path and name of the silva data base, as XX/XX/silva_nrXX_vXX_train_set.fa.gz"),
  make_option(c("-l", "--truncfor"), type="integer", default=NULL,
              help="To which length the forward reads should be truncated."),
  make_option(c("-r", "--truncrev"), type="integer", default=NULL,
              help="To which length the reverse reads should be truncated.")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# check inputs
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("An input folder must be provided.", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("An output folder must be provided.", call.=FALSE)
}
if (is.null(opt$db)){
  print_help(opt_parser)
  stop("The path to silva data base must be provided.", call.=FALSE)
}


####################
#Start DADA pipeline

# Path to the data
#setwd("~/16S_Amplicons")
path.in <- opt$input
dir.create(opt$output)
path.out <- paste(opt$output,"/",sep="")
#list.files(path)

###INSPECT READ QUALITY PROFILES
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fr <- paste("1.", opt$extention,sep="")
rr <- paste("2.", opt$extention,sep="")
fnFs <- sort(list.files(paste(path.in, sep=""), pattern=fr, full.names = TRUE))
fnRs <- sort(list.files(paste(path.in, sep=""), pattern=rr, full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, opt$name)
print(paste("ID example: ", sample.names[1], " out of ", length(sample.names)," samples.", sep=""))

#plot quality profiles of all samples in one pdf
plot_qp <- function(fn){
  plots.list <- list()
  i = 1
  while (i <= length(fn))
  {
    if (i < length(fn))
    {
      p.1 <- plotQualityProfile(fn[i:(i+1)])
   } else {
     p.1 <- plotQualityProfile(fn[i])
   }
    plots.list[[ceiling(i/2)]] <- p.1
    i <- i+2
  }
  return(plots.list)
}
plots.list <- plot_qp(fnFs)
ggsave(paste(path.out, "qulity_profiles_fr.pdf", sep=""), 
       marrangeGrob(grobs = plots.list, nrow=2, ncol=1))
plots.list <- plot_qp(fnRs)
ggsave(paste(path.out, "qulity_profiles_rr.pdf", sep=""), 
       marrangeGrob(grobs = plots.list, nrow=2, ncol=1))


###FILTER AND TRIM
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.out, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.out, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

if (is.null(opt$truncfor) || is.null(opt$truncrev)){
print("No truncation of the reads.")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=50,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
} else {
print(paste("Reads are truncated for:", opt$truncfor, " rev:", opt$truncrev,sep=""))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=50, 
                     truncLen=c(opt$truncfor,opt$truncrev),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
}

head(out)
write.table(out, file=paste(path.out,"trim_statistics.txt", sep=""), sep="\t", quote = F)

plots.list <- plot_qp(filtFs)
ggsave(paste(path.out, "qulity_profiles_filt_fr.pdf", sep=""), 
       marrangeGrob(grobs = plots.list, nrow=2, ncol=1))
plots.list <- plot_qp(filtRs)
ggsave(paste(path.out, "qulity_profiles_filt_rr.pdf", sep=""), 
       marrangeGrob(grobs = plots.list, nrow=2, ncol=1))

###LEARN THE ERROR RATES
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plot.error <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste(path.out, "errors_distribution_f.pdf", sep=""), plot.error, device="pdf")
plot.error <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste(path.out, "errors_distribution_r.pdf", sep=""), plot.error, device="pdf")


###SAMPLE INFERENCE
# apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object
dadaFs[[1]]

###MERGING PAIRED READS
#By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases,
#and are identical to each other in the overlap region (but these conditions can be changed via function arguments).
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

###CONSTRUCT SEQUENCE TABLE
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
seq.len <- table(nchar(getSequences(seqtab)))
write.table(seq.len, file=paste(path.out, "sequence_statistics.txt", sep=""), sep="\t", quote = F)
#Sequences that are much longer or shorter than expected may be the result of non-specific priming. 
#You can remove non-target-length sequences from your sequence table 
#(eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). 
#seqtab.lnorm <- seqtab[,nchar(colnames(seqtab)) %in% 457:459]
seqtab.lnorm <- seqtab

#write(paste("\nPercentage of remaining sequences after too short/long sequences were removed: ", sep=""), 
#            file=paste(path.out, "sequence_statistics.txt", sep=""), sep="/t", append = T)
#write(round(sum(seqtab.lnorm)/sum(seqtab)*100, digits=2), file=paste(path.out, "sequence_statistics.txt", sep=""), 
#      sep="/t", append = T)

###REMOVE CHIMERAS
seqtab.nochim <- removeBimeraDenovo(seqtab.lnorm, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

write(paste("\nPercentage of remaining sequences after chimeras were removed: ", sep=""), 
      file=paste(path.out, "sequence_statistics.txt", sep=""), sep="/t", append = T)
write(round(sum(seqtab.nochim)/sum(seqtab.lnorm)*100, digits=2), file=paste(path.out, "sequence_statistics.txt", sep=""), 
      sep="/t", append = T)

#save library
require(mgcv)
saveRDS(seqtab.nochim, paste(path.out,"seqtab_nochim.rds", sep=""))

seqtab.print <- seqtab.nochim # Removing sequence colnames for display only
colnames(seqtab.print) <- NULL
head(seqtab.print)[,1:5]

write.table(seqtab.print, file=paste(path.out, "ASV.txt", sep=""), sep="\t", quote = F)


###TRACK READS THROUGH THE PIPELINE, how many of them made it through?
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, file=paste(path.out, "reads_track.txt", sep=""), sep="\t", quote = F)

###ASSIGN TAXONOMY
#To follow along, download the silva_nr_v132_train_set.fa.gz file, and place it in the directory with the fastq files.
taxa <- assignTaxonomy(seqtab.nochim, opt$db, multithread=TRUE)
t2 <- strsplit(opt$db, "/")
t2 <- t2[[1]]
t3 <- paste(t2[1:(length(t2)-1)], collapse = "/")
t4 <- t2[length(t2)]
t4 <- sapply(strsplit(t4, "_"), `[`, 3)
t4 <- paste("silva_species_assignment_", t4, ".fa.gz", sep = "")
silva_species <- paste(t3,t4,sep = "/")

taxa <- addSpecies(taxa, silva_species)

#inspect the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, paste(path.out,"taxa.rds", sep=""))
write.table(taxa.print, file=paste(path.out, "taxa.txt", sep=""), sep="\t", quote = F)

##PRINT OUT ABBUNDANCE TABLES
#transpose seqtab_mouse_trial.rds
seqtab.trans <- as.data.frame(t(seqtab.nochim))
taxa.info <- as.data.frame(taxa)

#bind taxa info and number of reads per ASV 
taxa.seqtab <- merge(taxa.info, seqtab.trans, by="row.names", all=TRUE, sort=FALSE)
seqtab.bact <- subset(taxa.seqtab, Kingdom == "Bacteria")

#save sequences 
sequences <- as.data.frame(seqtab.bact$Row.names)
colnames(sequences)[1]<- "sequence"
seqtab.bact <- seqtab.bact[,-1]

#rename rows
sv.list <- c()
for (i in 1:nrow(seqtab.bact)) {
  sv.list[i] <- paste("ASV",i,sep = "")
}
head(sv.list)
rownames(seqtab.bact) <- sv.list

seqtab.final <- cbind(sequences, seqtab.bact)
saveRDS(seqtab.final, paste(path.out, "seqtab_final.rds", sep=""))

write.table(seqtab.bact, file=paste(path.out, "taxa_num_reads.txt", sep=""), sep="\t", quote = F)


#get relative abundance
seqtab.numb <- seqtab.bact[,-c(1:7)]
seqtab.ra <- scale(seqtab.numb, center=F, scale=colSums(seqtab.numb))
# check if each column has sum of 1; than it is correct
colSums(seqtab.ra)
seqtab.fra <- cbind(seqtab.bact[,c(1:7)], seqtab.ra)
saveRDS(seqtab.fra, paste(path.out, "seqtab_frac.rds", sep=""))
write.table(seqtab.fra, file=paste(path.out, "taxa_frac_reads.txt", sep=""), sep="\t", quote = F)


###HANDOFF TO PHYLOSEQ

theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)
subject <- samples.out
samdf <- data.frame(Subject=subject)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, paste(path.out, "phyloseq.rds", sep=""))

