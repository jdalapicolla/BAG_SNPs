#Script by Jeronymo Dalapicolla, 2021

###PRE-ANALYSIS ----

### Libraries ----
library(tidyverse)
library(r2vcftools)
library(vcfR)
library(adegenet)
library(ggrepel)
library(ggpubr)
library(PopGenReport)
library(strataG)


### Load auxilary functions ----
source("Utilities.R")

###1. LOAD FILES ----
#A. Project name
project_name = "pilocarpus"

#B. Load neutral .vcf file with genetic clusters ID:
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_LEA.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5266 SNPs.
names(snps_neutral@meta)

#C. Number of populations and method:
snps_neutral@meta
optimal_K = 4
method = "SNMF"

#D. If you need split samples by population. Position of samples by population by sNMF approach. Choose one method and change it on script:
for (i in 1:length(unique(snps_neutral@meta$PopID_snmf))){
  pop = which(snps_neutral@meta$PopID_snmf == i)
  pop_vcf = Subset(snps_neutral, samples=pop)
  assign(paste0("pop_SNMF_", i), pop)
  assign(paste0("pop", i), pop_vcf)
}

#E. Create a Genind object based on neutral vcfR
vcf_neutral = read.vcfR(paste0("vcf/", project_name,"_filtered_neutral_LEA.vcf"), verbose = FALSE)
#convert
genind = vcfR2genind(vcf_neutral)

#F. If you need split samples by population. Set populations
pop1_gi = genind[pop_SNMF_1,]
pop2_gi = genind[pop_SNMF_2,]
pop3_gi = genind[pop_SNMF_3,]
pop4_gi = genind[pop_SNMF_4,]

#G. Define maximum number of indivuidual by population
table(snps_neutral@meta$PopID_snmf)
# 1   2   3   4 
#47  65  47 118


### 2. ALLELE FREQUENCY USING GENIND OBJECTS-----

#2.1. Population 1 ----
#A. Define inputs for simulations:
genind = pop1_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications
copies = 5 # number of SNPs copies to be consider as rare allele 

#B. Run simulations:
Allele_freq_pop1 = boot.alleles (genind, range, nboot, copies, replacement =  FALSE)

#C. Verify and save as csv
Allele_freq_pop1
lapply(Allele_freq_pop1, function(x) write.table(data.frame(x), "Allele_Freq_POP1.csv", append = T, sep = ","))

#D. Calculate for Real Data:
total_freq_obs = as.data.frame(colSums(genind@tab, na.rm =T))
total_freq_obs = total_freq_obs[grep(pattern = '.1$', x = rownames(total_freq_obs), value = T),] #selecting reference allele
null_allele_pop1 = length(total_freq_obs[total_freq_obs  == 0]) # alleles in species not present in population
rare_allele_pop1 = length(total_freq_obs[total_freq_obs < copies]) # number of rare alleles

#E. Verify
null_allele_pop1
rare_allele_pop1

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = Allele_freq_pop1[[1]]#change pop
df$pvalue = round(Allele_freq_pop1[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=null_allele_pop1), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Null Alleles") +
  xlab("Sample Size")+
  ggtitle("Population 1") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

#define df for graphics for RARE ALLELES
df = Allele_freq_pop1[[3]]#change pop
df$pvalue = round(Allele_freq_pop1[[4]]$`p-value`,3)

plotB = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=rare_allele_pop1), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Rare Alleles") +
  xlab("Sample Size")+
  ggtitle("Population 1") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotB

# Save as pdf
pdf("Graphics_SeedBanks_AlleleFreq_POP1.pdf") #change the name
plotA
plotB
dev.off()




#2.2. Population 2 ----
#A. Define inputs for simulations:
genind = pop2_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45,50,55,60) #change the last number for sample size of wild population/species
nboot = 100 #replications
copies = 5 # number of SNPs copies to be consider as rare allele 

#B. Run simulations:
Allele_freq_pop2 = boot.alleles (genind, range, nboot, copies, replacement =  FALSE)

#C. Verify
Allele_freq_pop2
lapply(Allele_freq_pop2, function(x) write.table(data.frame(x), "Allele_Freq_POP2.csv", append = T, sep = ","))

#D. Calculate for Real Data:
total_freq_obs = as.data.frame(colSums(genind@tab, na.rm =T))
total_freq_obs = total_freq_obs[grep(pattern = '.1$', x = rownames(total_freq_obs), value = T),] #selecting reference allele
null_allele_pop2 = length(total_freq_obs[total_freq_obs  == 0]) # alleles in species not present in population
rare_allele_pop2 = length(total_freq_obs[total_freq_obs < copies]) # number of rare alleles

#E. Verify
null_allele_pop2
rare_allele_pop2

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = Allele_freq_pop2[[1]]#change pop
df$pvalue = round(Allele_freq_pop2[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=null_allele_pop2), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Null Alleles") +
  xlab("Sample Size")+
  ggtitle("Population 2") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

#define df for graphics for RARE ALLELES
df = Allele_freq_pop2[[3]]#change pop
df$pvalue = round(Allele_freq_pop2[[4]]$`p-value`,3)

plotB = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=rare_allele_pop2), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Rare Alleles") +
  xlab("Sample Size")+
  ggtitle("Population 2") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotB

# Save as pdf
pdf("Graphics_SeedBanks_AlleleFreq_POP2.pdf") #change the name
plotA
plotB
dev.off()



#2.3. Population 3 ----
#A. Define inputs for simulations:
genind = pop3_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications
copies = 5 # number of SNPs copies to be consider as rare allele 

#B. Run simulations:
Allele_freq_pop3 = boot.alleles (genind, range, nboot, copies, replacement =  FALSE)

#C. Verify
Allele_freq_pop3
lapply(Allele_freq_pop3, function(x) write.table(data.frame(x), "Allele_Freq_POP3.csv", append = T, sep = ","))

#D. Calculate for Real Data:
total_freq_obs = as.data.frame(colSums(genind@tab, na.rm =T))
total_freq_obs = total_freq_obs[grep(pattern = '.1$', x = rownames(total_freq_obs), value = T),] #selecting reference allele
null_allele_pop3 = length(total_freq_obs[total_freq_obs  == 0]) # alleles in species not present in population
rare_allele_pop3 = length(total_freq_obs[total_freq_obs < copies]) # number of rare alleles

#E. Verify
null_allele_pop3
rare_allele_pop3

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = Allele_freq_pop3[[1]]#change pop
df$pvalue = round(Allele_freq_pop3[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=null_allele_pop3), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Null Alleles") +
  xlab("Sample Size")+
  ggtitle("Population 3") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

#define df for graphics for RARE ALLELES
df = Allele_freq_pop3[[3]]#change pop
df$pvalue = round(Allele_freq_pop3[[4]]$`p-value`,3)

plotB = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=rare_allele_pop3), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Rare Alleles") +
  xlab("Sample Size")+
  ggtitle("Population 3") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotB

# Save as pdf
pdf("Graphics_SeedBanks_AlleleFreq_POP3.pdf") #change the name
plotA
plotB
dev.off()



#2.4. Population 4 ----
#A. Define inputs for simulations:
genind = pop4_gi #genind by population/species
range = c(10,20,30,40,50,60,70,80,90,100,110) #change the last number for sample size of wild population/species
nboot = 100 #replications
copies = 5 # number of SNPs copies to be consider as rare allele 

#B. Run simulations:
Allele_freq_pop4 = boot.alleles (genind, range, nboot, copies, replacement =  FALSE)

#C. Verify
Allele_freq_pop4
lapply(Allele_freq_pop4, function(x) write.table(data.frame(x), "Allele_Freq_POP4.csv", append = T, sep = ","))

#D. Calculate for Real Data:
total_freq_obs = as.data.frame(colSums(genind@tab, na.rm =T))
total_freq_obs = total_freq_obs[grep(pattern = '.1$', x = rownames(total_freq_obs), value = T),] #selecting reference allele
null_allele_pop4 = length(total_freq_obs[total_freq_obs  == 0]) # alleles in species not present in population
rare_allele_pop4 = length(total_freq_obs[total_freq_obs < copies]) # number of rare alleles

#E. Verify
null_allele_pop4
rare_allele_pop4

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = Allele_freq_pop4[[1]]#change pop
df$pvalue = round(Allele_freq_pop4[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=null_allele_pop4), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Null Alleles") +
  xlab("Sample Size")+
  ggtitle("Population 4") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

#define df for graphics for RARE ALLELES
df = Allele_freq_pop4[[3]]#change pop
df$pvalue = round(Allele_freq_pop4[[4]]$`p-value`,3)

plotB = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=rare_allele_pop4), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Rare Alleles") +
  xlab("Sample Size")+
  ggtitle("Population 4") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotB

# Save as pdf
pdf("Graphics_SeedBanks_AlleleFreq_POP4.pdf") #change the name
plotA
plotB
dev.off()











### 3. GENETIC DIVERSITY (HE) USING GENIND OBJECTS-----

#3.1. Population 1 ----
#A. Define inputs for simulations:
genind = pop1_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications


#B. Run simulations:
he_pop1 = boot.he (genind, range, nboot, replacement =  FALSE)

#C. Verify
he_pop1
lapply(he_pop1, function(x) write.table(data.frame(x), "Genetic_Diversity_POP1.csv", append = T, sep = ","))

#D. Calculate for Real Data:
he_obs_pop1 = Hs(genind)

#E. Verify
he_obs_pop1

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = he_pop1[[1]]#change pop
df$pvalue = round(he_pop1[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=he_obs_pop1), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("He") +
  xlab("Sample Size")+
  ggtitle("Population 1") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_GeneticDiversity_POP1.pdf") #change the name
plotA
dev.off()



#3.2. Population 2 ----
#A. Define inputs for simulations:
genind = pop2_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45,50,55,60) #change the last number for sample size of wild population/species
nboot = 100 #replications

#B. Run simulations:
he_pop2 = boot.he (genind, range, nboot, replacement =  FALSE)

#C. Verify
he_pop2
lapply(he_pop2, function(x) write.table(data.frame(x), "Genetic_Diversity_POP2.csv", append = T, sep = ","))

#D. Calculate for Real Data:
he_obs_pop2 = Hs(genind)

#E. Verify
he_obs_pop2

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = he_pop2[[1]]#change pop
df$pvalue = round(he_pop2[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=he_obs_pop2), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("He") +
  xlab("Sample Size")+
  ggtitle("Population 2") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_GeneticDiversity_POP2.pdf") #change the name
plotA
dev.off()





#3.3. Population 3 ----
#A. Define inputs for simulations:
genind = pop3_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications

#B. Run simulations:
he_pop3 = boot.he (genind, range, nboot, replacement =  FALSE)

#C. Verify
he_pop3
lapply(he_pop3, function(x) write.table(data.frame(x), "Genetic_Diversity_POP3.csv", append = T, sep = ","))

#D. Calculate for Real Data:
he_obs_pop3 = Hs(genind)

#E. Verify
he_obs_pop3

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = he_pop3[[1]]#change pop
df$pvalue = round(he_pop3[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=he_obs_pop3), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("He") +
  xlab("Sample Size")+
  ggtitle("Population 3") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_GeneticDiversity_POP3.pdf") #change the name
plotA
dev.off()



#3.4. Population 4 ----
#A. Define inputs for simulations:
genind = pop4_gi #genind by population/species
range = c(10,20,30,40,50,60,70,80,90,100,110) #change the last number for sample size of wild population/species
nboot = 100 #replications

#B. Run simulations:
he_pop4 = boot.he (genind, range, nboot, replacement =  FALSE)

#C. Verify
he_pop4
lapply(he_pop4, function(x) write.table(data.frame(x), "Genetic_Diversity_POP4.csv", append = T, sep = ","))

#D. Calculate for Real Data:
he_obs_pop4 = Hs(genind)

#E. Verify
he_obs_pop4

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = he_pop4[[1]]#change pop
df$pvalue = round(he_pop4[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=he_obs_pop4), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("He") +
  xlab("Sample Size")+
  ggtitle("Population 4") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_GeneticDiversity_POP4.pdf") #change the name
plotA
dev.off()









### 4. ALLELE RICHNESS (AR) USING GENIND OBJECTS-----

#4.1. Population 1 ----
#A. Define inputs for simulations:
genind = pop1_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications
pop = rep(1,47) #population assessement

#B. Run simulations:
AR_pop1 = boot.AR (genind, range, nboot, pop, replacement =  FALSE)
warnings()

#C. Verify
AR_pop1
lapply(AR_pop1, function(x) write.table(data.frame(x), "Allele_Richness_POP1.csv", append = T, sep = ","))

#D. Calculate for Real Data:
genind@pop = as.factor(pop)
alle_rich_obs = allel.rich(genind, min.alleles = NULL)
AR_obs_pop1 = alle_rich_obs$mean.richness

#E. Verify
AR_obs_pop1

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = AR_pop1[[1]]#change pop
df$pvalue = round(AR_pop1[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=AR_obs_pop1), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Allele Richness") +
  xlab("Sample Size")+
  ggtitle("Population 1") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
 # geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_AlleleRichness_POP1.pdf") #change the name
plotA
dev.off()









#4.2. Population 2 ----
#A. Define inputs for simulations:
genind = pop2_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45,50,55,60) #change the last number for sample size of wild population/species
nboot = 100 #replications
pop = rep(1,65) #population assessement

#B. Run simulations:
AR_pop2 = boot.AR (genind, range, nboot, pop, replacement =  FALSE)


#C. Verify
AR_pop2
lapply(AR_pop2, function(x) write.table(data.frame(x), "Allele_Richness_POP2.csv", append = T, sep = ","))

#D. Calculate for Real Data:
genind@pop = as.factor(pop)
alle_rich_obs = allel.rich(genind, min.alleles = NULL)
AR_obs_pop2 = alle_rich_obs$mean.richness

#E. Verify
AR_obs_pop2

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = AR_pop2[[1]]#change pop
df$pvalue = round(AR_pop2[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=AR_obs_pop2), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Allele Richness") +
  xlab("Sample Size")+
  ggtitle("Population 2") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_AlleleRichness_POP2.pdf") #change the name
plotA
dev.off()






#4.3. Population 3 ----
#A. Define inputs for simulations:
genind = pop3_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications
pop = rep(1,47) #population assessement

#B. Run simulations:
AR_pop3 = boot.AR (genind, range, nboot, pop, replacement =  FALSE)

#C. Verify
AR_pop3
lapply(AR_pop3, function(x) write.table(data.frame(x), "Allele_Richness_POP3.csv", append = T, sep = ","))

#D. Calculate for Real Data:
genind@pop = as.factor(pop)
alle_rich_obs = allel.rich(genind, min.alleles = NULL)
AR_obs_pop3 = alle_rich_obs$mean.richness

#E. Verify
AR_obs_pop3

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = AR_pop3[[1]]#change pop
df$pvalue = round(AR_pop3[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=AR_obs_pop3), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Allele Richness") +
  xlab("Sample Size")+
  ggtitle("Population 3") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_AlleleRichness_POP3.pdf") #change the name
plotA
dev.off()





#4.4. Population 4 ----
#A. Define inputs for simulations:
genind = pop4_gi #genind by population/species
range =  c(10,20,30,40,50,60,70,80,90,100,110) #change the last number for sample size of wild population/species
nboot = 100 #replications
pop = rep(1,118) #population assessement

#B. Run simulations:
AR_pop4 = boot.AR (genind, range, nboot, pop, replacement =  FALSE)

#C. Verify
AR_pop4
lapply(AR_pop4, function(x) write.table(data.frame(x), "Allele_Richness_POP4.csv", append = T, sep = ","))

#D. Calculate for Real Data:
genind@pop = as.factor(pop)
alle_rich_obs = allel.rich(genind, min.alleles = NULL)
AR_obs_pop4 = alle_rich_obs$mean.richness

#E. Verify
AR_obs_pop4

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = AR_pop4[[1]]#change pop
df$pvalue = round(AR_pop4[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=AR_obs_pop4), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Allele Richness") +
  xlab("Sample Size")+
  ggtitle("Population 4") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_AlleleRichness_POP4.pdf") #change the name
plotA
dev.off()







### 5. NUCLEOTIDE DIVERSITY (PI) USING VCFLINK OBJECTS-----

#5.1. Population 1 ----
#A. Define inputs for simulations:
vcf = pop1 #VCF by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications

#B. Run simulations:
PI_pop1 = boot.pi (vcf, range, nboot, replacement =  FALSE)


#C. Verify
PI_pop1
lapply(PI_pop1, function(x) write.table(data.frame(x), "Nucleodite_Diversity_POP1.csv", append = T, sep = ","))

#D. Calculate for Real Data:
PI_obs = Query(vcf, type="site-pi")
PI_obs_pop1 = mean(PI_obs$PI, na.rm = T)

#E. Verify
PI_obs_pop1

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = PI_pop1[[1]]#change pop
df$pvalue = round(PI_pop1[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=PI_obs_pop1), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Nucleotide Diversity") +
  xlab("Sample Size")+
  ggtitle("Population 1") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
#  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_NucleotideDiversity_POP1.pdf") #change the name
plotA
dev.off()








#5.2. Population 2 ----
#A. Define inputs for simulations:
vcf = pop2 #VCF by population/species
range = c(5,10,15,20,25,30,35,40,45,50,55,60) #change the last number for sample size of wild population/species
nboot = 100 #replications

#B. Run simulations:
PI_pop2 = boot.pi (vcf, range, nboot, replacement =  FALSE)

#C. Verify
PI_pop2
lapply(PI_pop2, function(x) write.table(data.frame(x), "Nucleodite_Diversity_POP2.csv", append = T, sep = ","))

#D. Calculate for Real Data:
PI_obs = Query(vcf, type="site-pi")
PI_obs_pop2 = mean(PI_obs$PI, na.rm = T)

#E. Verify
PI_obs_pop2

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = PI_pop2[[1]]#change pop
df$pvalue = round(PI_pop2[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=PI_obs_pop2), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Nucleotide Diversity") +
  xlab("Sample Size")+
  ggtitle("Population 2") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_NucleotideDiversity_POP2.pdf") #change the name
plotA
dev.off()


#5.3. Population 3 ----
#A. Define inputs for simulations:
vcf = pop3 #VCF by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications

#B. Run simulations:
PI_pop3 = boot.pi (vcf, range, nboot, replacement =  FALSE)


#C. Verify
PI_pop3
lapply(PI_pop3, function(x) write.table(data.frame(x), "Nucleodite_Diversity_POP3.csv", append = T, sep = ","))

#D. Calculate for Real Data:
PI_obs = Query(vcf, type="site-pi")
PI_obs_pop3 = mean(PI_obs$PI, na.rm = T)

#E. Verify
PI_obs_pop3

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = PI_pop3[[1]]#change pop
df$pvalue = round(PI_pop3[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=PI_obs_pop3), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Nucleotide Diversity") +
  xlab("Sample Size")+
  ggtitle("Population 3") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_NucleotideDiversity_POP3.pdf") #change the name
plotA
dev.off()



#5.4. Population 4 ----
#A. Define inputs for simulations:
vcf = pop4 #VCF by population/species
range = c(10,20,30,40,50,60,70,80,90,100,110) #change the last number for sample size of wild population/species
nboot = 100 #replications

#B. Run simulations:
PI_pop4 = boot.pi (vcf, range, nboot, replacement =  FALSE)

#C. Verify
PI_pop4
lapply(PI_pop4, function(x) write.table(data.frame(x), "Nucleodite_Diversity_POP4.csv", append = T, sep = ","))

#D. Calculate for Real Data:
PI_obs = Query(vcf, type="site-pi")
PI_obs_pop4 = mean(PI_obs$PI, na.rm = T)

#E. Verify
PI_obs_pop4

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = PI_pop4[[1]]#change pop
df$pvalue = round(PI_pop4[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=PI_obs_pop4), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Nucleotide Diversity") +
  xlab("Sample Size")+
  ggtitle("Population 4") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_NucleotideDiversity_POP4.pdf") #change the name
plotA
dev.off()




### 6. POPULATION SIZE (Ne) USING GENIND OBJECTS-----

#6.1 Population 1 ----
#A. Define inputs for simulations:
genind = pop1_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications


#B. Run simulations:
Ne_pop1 = boot.Ne (genind, range, nboot, replacement =  FALSE)


#C. Verify
Ne_pop1
lapply(Ne_pop1, function(x) write.table(data.frame(x), "Population_Sizes_POP1.csv", append = T, sep = ","))

#D. Calculate for Real Data:
gtypes = genind2gtypes(genind)
Ne_obs = ldNe(gtypes, maf.threshold = 0, by.strata = FALSE, ci = 0.95, drop.missing = TRUE, num.cores = 3)
Ne_obs_pop1 = Ne_obs[, c(6)] 
Ne_obs_pop1

#E. Verify
Ne_obs_pop1

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = Ne_pop1[[1]]#change pop
df$pvalue = round(Ne_pop1[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=Ne_obs_pop1), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Effective Population Size (Ne)") +
  xlab("Sample Size")+
  ggtitle("Population 1") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  #geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_PopulationSize_POP1.pdf") #change the name
plotA
dev.off()




















#6.2 Population 2 ----
#A. Define inputs for simulations:
genind = pop2_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45,50,55,60) #change the last number for sample size of wild population/species
nboot = 100 #replications


#B. Run simulations:
Ne_pop2 = boot.Ne (genind, range, nboot, replacement =  FALSE)


#C. Verify
Ne_pop2
lapply(Ne_pop2, function(x) write.table(data.frame(x), "Population_Sizes_POP2.csv", append = T, sep = ","))

#D. Calculate for Real Data:
gtypes = genind2gtypes(genind)
Ne_obs = ldNe(gtypes, maf.threshold = 0, by.strata = FALSE, ci = 0.95, drop.missing = TRUE, num.cores = 3)
Ne_obs_pop2 = Ne_obs[, c(6)] 
Ne_obs_pop2

#E. Verify
Ne_obs_pop2

#D. Graphics for Results
#define df for graphics
df = Ne_pop2[[1]]#change pop
df$pvalue = round(Ne_pop2[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=Ne_obs_pop2), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Effective Population Size (Ne)") +
  xlab("Sample Size")+
  ggtitle("Population 2") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_PopulationSize_POP2.pdf") #change the name
plotA
dev.off()





#6.3 Population 3 ----
#A. Define inputs for simulations:
genind = pop3_gi #genind by population/species
range = c(5,10,15,20,25,30,35,40,45) #change the last number for sample size of wild population/species
nboot = 100 #replications


#B. Run simulations:
Ne_pop3 = boot.Ne (genind, range, nboot, replacement =  FALSE)

#C. Verify
Ne_pop3
lapply(Ne_pop3, function(x) write.table(data.frame(x), "Population_Sizes_POP3.csv", append = T, sep = ","))

#D. Calculate for Real Data:
gtypes = genind2gtypes(genind)
Ne_obs = ldNe(gtypes, maf.threshold = 0, by.strata = FALSE, ci = 0.95, drop.missing = TRUE, num.cores = 3)
Ne_obs_pop3 = Ne_obs[, c(6)] 
Ne_obs_pop3

#E. Verify
Ne_obs_pop3

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = Ne_pop3[[1]]#change pop
df$pvalue = round(Ne_pop3[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=Ne_obs_pop3), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Effective Population Size (Ne)") +
  xlab("Sample Size")+
  ggtitle("Population 3") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_PopulationSize_POP3.pdf") #change the name
plotA
dev.off()








#6.4 Population 4 ----
#A. Define inputs for simulations:
genind = pop4_gi #genind by population/species
range = c(10,20,30,40,50,60,70,80,90,100,110) #change the last number for sample size of wild population/species
nboot = 100 #replications

#B. Run simulations:
Ne_pop4 = boot.Ne (genind, range, nboot, replacement =  FALSE)

#C. Verify
Ne_pop4
lapply(Ne_pop4, function(x) write.table(data.frame(x), "Population_Sizes_POP4.csv", append = T, sep = ","))

#D. Calculate for Real Data:
gtypes = genind2gtypes(genind)
Ne_obs = ldNe(gtypes, maf.threshold = 0, by.strata = FALSE, ci = 0.95, drop.missing = TRUE, num.cores = 3)
Ne_obs_pop4 = Ne_obs[, c(6)] 
Ne_obs_pop4

#E. Verify
Ne_obs_pop4

#D. Graphics for Results
#define df for graphics for NULL ALLELES
df = Ne_pop4[[1]]#change pop
df$pvalue = round(Ne_pop4[[2]]$`p-value`,3)

plotA = ggplot(df, aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=2) +
  geom_point(size=3.5) +
  geom_hline(mapping = aes(yintercept=Ne_obs_pop4), colour="red", linetype="dashed") + #change pop number
  theme_bw() +
  ylab("Effective Population Size (Ne)") +
  xlab("Sample Size")+
  ggtitle("Population 4") + #change pop number
  scale_x_continuous(breaks = range, labels = as.character(range)) +
  geom_label_repel(data= df, mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold") + # in case to show p-values 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16),
        axis.title = element_text(hjust = 0.5, face="bold", size = 12))

#Verify
plotA

# Save as pdf
pdf("Graphics_SeedBanks_PopulationSize_POP4.pdf") #change the name
plotA
dev.off()

##END
