##Script by Jeronymo Dalapicolla, 2020

##-------------------------------Libraries
library(tidyverse)
library(r2vcftools)
library(vcfR)
library(adegenet)
library(ggrepel)


##---------------------------------Inputs
#Load functions
source("Utilities.R")

# Project name:
project_name = "ipomoea"

# DOWNLOAD A VCF FILE AS EXAMPLE "Icavalcantei.vcf" FROM FIGSHARE: https://doi.org/10.6084/m9.figshare.6100004.v1
#A. CREATE A FOLDER NAMED "vcf" IN YOUR WORKING DIRECTORY AND SAVE THE .vcf FILE THERE and cahnge de name if you prefer.

# Load neutral .vcf file with genetic clusters ID:
snps_neutral = vcfLink(paste0("vcf/", project_name, ".vcf"), overwriteID=T)
VCFsummary(snps_neutral) #115 individuals and 13167 SNPs.

# Number of populations and method:
snps_neutral@meta
optimal_K = 1
method = "SNMF"

#If you need split samples by population. Position of samples by population by sNMF approach. Choose one method and change it on script:
for (i in 1:length(unique(snps_neutral@meta$PopID_snmf))){
  pop = which(snps_neutral@meta$PopID_snmf == i)
  pop_vcf = Subset(snps_neutral, samples=pop)
  assign(paste0("pop_SNMF_", i), pop)
  assign(paste0("pop", i), pop_vcf)
}


#Create a Genind object based on neutral vcfR
vcf_neutral = read.vcfR(paste0("vcf/", project_name,".vcf"), verbose = FALSE)
#convert
genind = vcfR2genind(vcf_neutral)

#If you need split samples by population. Set populations
pop1_gi = genind[pop_SNMF_1,]
pop2_gi = genind[pop_SNMF_2,]
pop3_gi = genind[pop_SNMF_3,]
pop4_gi = genind[pop_SNMF_4,]


##################################################################################
####----------- 1. Allelic Frequency in Wild Species/Population
##################################################################################
      
#define datasets and their names by population or species, in the case, 1 population
y = snps_neutral
freq.categories = c("NA", "0.05", "0.1", "0.25", "0.5") #1st category is always NA

#calculate frequencies in original populations
freq_ipomoea = allele.freq(y, freq.categories)
#verify
freq_ipomoea

#create a final tibble
freq_obser = freq_ipomoea %>% ##change dataset name here
  add_column(n_sample = c(length(snps_neutral@sample_id))) %>% #change dataset name here
  rename(obser = n) #change dataset name here
#verify
freq_obser


##################################################################################
####----------- 2. Allelic Frequency in Resampled Datasets
##################################################################################

## Analyses by Sample Size, one population by turn (estimated). Maximum range is the original sample size, in my case 115:
#define inputs:
y = snps_neutral #vcfLink by population
freq.categories = c("NA", "0.05", "0.1", "0.25", "0.5") #1st category is always NA
range= c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110) #sample size to be tested
nboot = 100 #replications

ipomoea_resampled = pmt.allele.freq(y, range, nboot, freq.categories)


##--------------------------------------Graphics for results

pdf("Graphics_Ipomoea_Freq_100Permutation.pdf")

ggplot(ipomoea_resampled, aes(n_sample, smean)) + #change data name
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1.5) + #change for min and max values
  geom_point(size=1) +
  #geom_point(data=freq_obser, aes(n_sample, y=obser), colour="red", size=3) +
  geom_hline(data=freq_obser, aes(yintercept=obser), colour="red", linetype="dashed") +
  facet_wrap(vars(class), scales = "free")+
  theme_bw() +
  ylab("Number of SNPs per Frequency Category") +
  xlab("Sample Size")+
  ggtitle("Ipomoea - 100 Permutations w/o replacement") + #change pop
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

dev.off()


##################################################################################
#### ----------- 3. Genetic Diversity in Resampled Datasets
##################################################################################

#Define range of sample size to be tested.
#define inputs:
genind = genind #genind by population/species
range = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 115) #change the last number for sample size of wild population/species
nboot = 100 #replications

#This function has "replace = TRUE", for permutation (w/o replacement) edit Utilities.R file to "replace = FALSE"
he_boots_ipomoea = boot.he (genind, range, nboot)

#verify
he_boots_ipomoea

##-------------------------------------- Graphics for results
#define df for graphics
df = he_boots_ipomoea[[1]]#change pop
df$pvalue = round(he_boots_ipomoea[[2]]$`p-value`,3)

pdf("Graphics_Ipomoea_He_100Bootstrap.pdf") #change the name

ggplot(df[-length(df[,1]),], aes(n_sample, mean)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=1.5) +
  geom_point(size=3) +
  geom_hline(df[length(df[,1]),], mapping = aes(yintercept=min), colour="red", linetype="dashed") +
  geom_hline(df[length(df[,1]),], mapping = aes(yintercept=max), colour="red", linetype="dashed") +
  #geom_hline(df[length(df[,1]),], mapping = aes(yintercept=mean), colour="black", linetype="solid") +
  theme_bw() +
  ylab("He Exp (Mean of 100 Bootstrap)") +
  xlab("Sample Size")+
  ggtitle("POP1 - 100 Bootstrap w/ replacement") + #change the name
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  geom_label_repel(data= df[-length(df[,1]),], mapping = aes(x=n_sample, y=mean, label=pvalue), fontface="bold")
  #geom_label_repel(data = subset(df[-length(df[,1]),], df[-length(df[,1]),][,5] > 0.05), mapping = aes(x=n_sample, y=mean, label=df$pvalue[df$pvalue > 0.05][-5]), fontface="bold", color ="blue") +
  #geom_label_repel(data = subset(df[-length(df[,1]),], df[-length(df[,1]),][,5] <= 0.05), mapping = aes(x=n_sample, y=mean, label=df$pvalue[df$pvalue <= 0.05]), fontface="bold", color ="red")

dev.off()

#ENDED
