library(dplyr)
library(vcfR)
library(data.table)

## Interested in Family I73T: Samples DV02-DV06
setwd("L:/medex/askol/Hamvas/NUSeqData/Hamvas01")

u <- fread(file = "GATK/Hamvas01.hg19_multianno.txt", header=T, fill=T, quote="")

## Based on the shell script for joining the individual vcf files, we are 
## assigning IDs to u based on the ording in the shell script.
## DV01, DV02, DV03, DV04, DV05, DV06, DV07, DV08, "DV09, DV10, DV11, DV12
ids <- c("DV01", "DV02", "DV03", "DV04", "DV05", "DV06", "DV07", "DV08", "DV09", "DV10", "DV11", "DV12")

names(u)[ncol(u) + 1 - length(ids):1] <- ids

## probably want to use files in Mol"Dx_Alamut_Annotated/
files <- dir("MolDx_Alamut_Annotated/", full.names = TRUE)

genesInterest <- read.table("../ListofRelevantProteinsandGenes-RR2-4-22.txt", as.is=TRUE)
genesInterest <- unlist(genesInterest)
allData <- c()

for (file in files){
  
  cat("Working on file ", file, "\n\n")
  
  data <- fread(file = file, header=T, fill=T, quote="")
  
  data <- data %>% select(chrom, inputPos, inputRef, inputAlt, gene, varType, codingEffect, varLocation, omimId, gnomadAltFreq_all)
  
  data <- data %>% filter(gene %in% genesInterest)
  
  allData <- rbind(allData, data)   
  
}

allData <- unique(allData)
allData <- allData %>% mutate(chrom = paste0("chr",chrom)) %>% filter(varLocation =="exon") %>% filter(codingEffect != "synonymous")


uFiltered <- left_join(allData, u, by = c("chrom"="Chr", "inputPos"="Start", "inputRef"="Ref", "inputAlt"="Alt"))

uFiltered <- uFiltered %>% select(chrom, inputPos, inputRef, inputAlt, gene, varType, codingEffect, varLocation, gnomadAltFreq_all, 
                                  CADD_raw, CADD_phred, DV02, DV03, DV04, DV05, DV06)

for (id in ids[2:6]){
  
  newColName = paste0(id,"_geno")
  uFiltered <- uFiltered %>% mutate(!!as.name(newColName) = getGene(!!as.name(id)))
}


getGene <- function(x){
  
  geno <- gsub(":.+", "", x)
}