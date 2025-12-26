# Danshen-analysis
#Reproducible code for: Streptomyces aurantiacus YS-4 / metacycloprodigiosin / Danshen root rot resistance
#DADA2 Analysis#
#2022/11/30_LuyangSong#
library(dada2)
#library(dplyr)
#library(ggplot2)
#library(DECIPHER)
#library(jsonlite)
packageVersion("dada2")
loger=function(msg){
  t=strftime(Sys.time(),"%Y-%m-%d %H-%M-%S")
cat("[",t,"]",msg,"\n")
}
#file parsing
loger("file parsing")
pathF<-"/public/home/qqwang/SLY/data2_any/danshen_raw_data/Bac_F"
pathR<-"/public/home/qqwang/SLY/data2_any/danshen_raw_data/Bac_R"
filtpathF<-file.path(pathF,"filtered")
filtpathR<-file.path(pathR,"filtered")
fastqFs<-sort(list.files(pathF,pattern = "fq.gz"))
fastqRs<-sort(list.files(pathR,pattern = "fq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
#filter sequence
out<-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(220,220), maxEE=c(2,2), truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
#INFER SEQUENCE VARIANTS
filtpathF<-"/public/home/qqwang/SLY/data2_any/danshen_raw_data/Bac_F/filtered"
filtpathR<-"/public/home/qqwang/SLY/data2_any/danshen_raw_data/Bac_R/filtered"
filtFs<-list.files(filtpathF,pattern = ".fq.gz",full.names = TRUE)
filtRs<-list.files(filtpathR,pattern = ".fq.gz",full.names = TRUE)

sample_names<-sapply(strsplit(basename(filtFs),"_"),'[',1)

sample_namesR<-sapply(strsplit(basename(filtRs),"_"),'[',1)

if(!identical(sample_names,sample_namesR)) stop("Forward and reverse files do not match.")


names(filtFs)<-sample_names
names(filtRs) <- sample_names
set.seed(100)
# Learn forward error rates
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=2e6, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=2e6, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample_names))
names(mergers) <- sample_names

for(sam in sample_names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  cat("Processing:", sam, "denoisedF",sum(getUniques(ddF)), "\n")
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  cat("Processing:", sam, "denoisedR",sum(getUniques(ddR)), "\n")
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
  cat("Processing:", sam, "merged",sum(getUniques(merger)), "\n")
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
dim(seqtab)#sample and ASVs number
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#plot sequence distribution
#plot(table(nchar(getSequences(seqtab))))

saveRDS(seqtab, "/public/home/qqwang/SLY/data2_any/danshen_raw_data/seqtab.rds") # CHANGE ME to where you want sequence table saved

# Remove chimeras
seqtab2 <- removeBimeraDenovo(seqtab,method="consensus", multithread=TRUE)
#cacluate removed chimeras ratio
sum(seqtab2)/sum(seqtab)

# Assign taxonomy
tax <- assignTaxonomy(seqtab2,"/public/home/qqwang/SLY/data2_any/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)#????ע?͵?????ˮƽ
taxa<-addSpecies(tax,"/public/home/qqwang/SLY/data2_any/silva_species_assignment_v138.1.fa.gz")
# Write to disk
saveRDS(out,"/public/home/qqwang/SLY/data2_any/danshen_raw_data/out.rds")
saveRDS(seqtab2,"/public/home/qqwang/SLY/data2_any/danshen_raw_data/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax,"/public/home/qqwang/SLY/data2_any/danshen_raw_data/tax_final.rds") # CHANGE ME ...
saveRDS(taxa, "/public/home/qqwang/SLY/data2_any/danshen_raw_data/taxa_final.rds")
# track reads through the pipeline; looking at the number of 
#reads that progressed to each step
#```{r}
#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, rowSums(seqtab2))
#saveRDS(taxa, "/public/home/wyzhang/Bac_R/track.rds")
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#colnames(track) <- c("input","nonchim")
#rownames(track) <- sample_names
#write.csv(track,"/public/home/wyzhang/Bac_R/track.csv")
save.image("/public/home/qqwang/SLY/data2_any/danshen_raw_data/DADA2_ASVs.RData")
#head(track)
#```




