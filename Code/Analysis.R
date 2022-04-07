library(dplyr)
library(vcfR)
library(data.table)
library(tidyverse)
library(ggtext)

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

genesInterest <- read.table("C:/Users/askol/Dropbox/GitHub/spc-variant-project/ListofRelevantProteinsandGenes-RR3-17-22.txt", as.is=TRUE)
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
  uFiltered <- uFiltered %>% mutate(!!as.name(newColName) := getGene(!!as.name(id)))
}

input <-  uFiltered %>% select(gnomadAltFreq_all, DV02_geno:DV06_geno)
input <- input %>% mutate(gnomadAltFreq_all =
                                ifelse(is.na(gnomadAltFreq_all), 0, gnomadAltFreq_all))
ind <- which(rowSums(is.na(input[,-1])) == (ncol(input)-1) )
input <- input[-ind, ]

weights <- as.numeric(uFiltered$CADD_phred[-ind])
weights <- weights / sum(weights, na.rm=T)

Scores <- calcGenoScore(input, freqCutoffs = c(.01, .05, .1, .2))

## Show Genotypes
genotypes_data <- uFiltered %>% select(inputRef, inputAlt, gene, gnomadAltFreq_all, DV02_geno:DV06_geno) %>% mutate(varName = paste(gene, inputRef, inputAlt, sep = "")) %>% select(varName, gnomadAltFreq_all, DV02_geno:DV06_geno)
genotypes_data <- genotypes_data %>% mutate(gnomadAltFreq_all = ifelse(is.na(gnomadAltFreq_all), 0, gnomadAltFreq_all))
na_geno_indices <- which(rowSums(is.na(genotypes_data[,-1])) == (ncol(genotypes_data)-2) )
genotypes_data = genotypes_data[-na_geno_indices, ]

genotypes_data <- genotypes_data %>%
  pivot_longer(cols = starts_with("DV"), names_to = "DV_genotypes", values_to = "genotype") %>%
  mutate(DV_genotypes = as.character(DV_genotypes), genotype = as.character(genotype))


showGenos(genotypes_data)

getGene <- function(x){

    tranTable <- rbind( c("0/0", 0), c("0/1", 1), c("1/1", 2))
    geno <- gsub(":.+", "", x)
    for (i in 1:nrow(tranTable)){

        ind <- which(geno == tranTable[i,1])
        geno[ind] <- as.numeric(tranTable[i,2])
    }
    return(geno)
}


calcGenoScore <- function(x, freqCutoffs, weights=c()){

    scores <- c()
 
    for (freq in freqCutoffs){
        
        z <- x %>% filter(gnomadAltFreq_all <= freq)
        if (length(weights) == nrow(x) ){
            w <- weights[x$gnomadAltFreq_all <= freq]
            w <- kronecker(w, t(rep(1, ncol(z)-1)))
        }
        
        if (nrow(z) > 0){
            sc <- z %>% select(-gnomadAltFreq_all)
            sc <- apply(sc, 2, as.numeric)
            if (length(weights) > 0){
                sc <- 1 * (sc > 0) * w
            }
            sc <- colSums(sc, na.rm=T)
            scores <- rbind(scores, c(freq, sc))
        }else{
            scores <- rbind(scores, c(freq, rep(0, ncol(x) - 1)))
        }
    }
        
    return(scores)
}
    

showGenos <- function(x){
    
  p <- ggplot(data = a, aes(x = varName, y = DV_genotypes, size = genotype, color = gnomadAltFreq_all)) + geom_point() + theme(axis.text.x = element_markdown(angle=90))
  
  
  return(p)
}
