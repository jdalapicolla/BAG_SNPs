##--------------------------------Functions
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

allele.freq = function(y, freq.categories){
  freq = Query(y, type="freq2")
  
  freq_result = freq %>%
    rownames_to_column(var = "snp_index") %>%
    mutate(class = case_when(.$N_CHR == 0 ~ freq.categories[1],
                             .$X.FREQ. < as.numeric(freq.categories[2]) ~ paste0("< ",freq.categories[2]),
                             .$X.FREQ. >= as.numeric(freq.categories[2]) & .$X.FREQ. < as.numeric(freq.categories[3]) ~ paste0("< ",freq.categories[3]),
                             .$X.FREQ. >= as.numeric(freq.categories[3]) & .$X.FREQ. < as.numeric(freq.categories[4]) ~ paste0("< ",freq.categories[4]),
                             .$X.FREQ. >= as.numeric(freq.categories[4]) & .$X.FREQ. < as.numeric(freq.categories[5]) ~ paste0("< ",freq.categories[5]),
                             .$X.FREQ. >= as.numeric(freq.categories[5]) ~ paste0("> ",freq.categories[5])
    )) %>%
    group_by(class, .drop=F) %>%
    count(.drop=F)
  
  return(freq_result)
}

pmt.allele.freq = function(y, range, nboot, freq.categories) {
  
  
  #lists to save results
  resB = list() #save results by permutation
  resS = list() #save results by sample size
  
  #name the categories:
  categories = c(freq.categories[1],
                 paste0("< ",freq.categories[2]),
                 paste0("< ",freq.categories[3]),
                 paste0("< ",freq.categories[4]),
                 paste0("< ",freq.categories[5]),
                 paste0("> ",freq.categories[5])
  )
  
  #counting loops:
  loop = 0
  
  for(i in range) { #vary sample size
    loop=loop+1 #counting number of different sample size tested
    res_sample = tibble(as.data.frame(matrix(categories,6, 1))) %>% # save results by sample size
      rename(class = V1)
    
    #loop for bootstrap
    for(k in 1:nboot){ #replications
      positions = sample(c(1:length(y@sample_id)), i, replace = F) #sample random indiv
      UNIND = y@sample_id[positions] #select indiv
      pop_subset = Subset(y, samples=UNIND) # subset in a vcf
      freq = Query(pop_subset, type="freq2") #estimating allelic frequencies
      
      res_freq = freq %>%
        rownames_to_column(var = "snp_index") %>%
        mutate(class = case_when(.$N_CHR == 0 ~ freq.categories[1],
                                 .$X.FREQ. < as.numeric(freq.categories[2]) ~ paste0("< ",freq.categories[2]),
                                 .$X.FREQ. >= as.numeric(freq.categories[2]) & .$X.FREQ. < as.numeric(freq.categories[3]) ~ paste0("< ",freq.categories[3]),
                                 .$X.FREQ. >= as.numeric(freq.categories[3]) & .$X.FREQ. < as.numeric(freq.categories[4]) ~ paste0("< ",freq.categories[4]),
                                 .$X.FREQ. >= as.numeric(freq.categories[4]) & .$X.FREQ. < as.numeric(freq.categories[5]) ~ paste0("< ",freq.categories[5]),
                                 .$X.FREQ. >= as.numeric(freq.categories[5]) ~ paste0("> ",freq.categories[5])
        )) %>%
        group_by(class, .drop=F) %>%
        count(.drop=F)
      
      res_sample = res_sample %>%
        full_join (res_freq, by = "class")
      
    }
    
    #permutations by a sample size in a single table
    tab = as.data.frame(dplyr::bind_rows(res_sample, .id = NULL)) 
    
    #Transpose the results
    tab = t(tab)
    colnames(tab) = tab[1,]
    rownames(tab) = NULL
    tab = tab[-1,]
    tab = tibble(as.data.frame(tab))
    tab = tab %>%
      mutate_if(is.factor, ~ as.numeric(as.character(.)))
    
    #summarizing the results
    result = bind_rows(summarise_all(tab, mean), summarise_all(tab, min), summarise_all(tab, max) , summarise_all(tab, sd), summarise_all(tab, length))
    result=t(result)
    colnames(result) = c("smean", "min2", "max2", "ssd", "count")
    result = tibble(as.data.frame(result))
    
    #Calculate CI and select columns to final results
    resS[[loop]] = result %>%
      mutate(class = categories) %>%
      mutate_if(is.factor, ~ as.numeric(as.character(.))) %>%
      mutate(se = ssd / sqrt(count),
             lower_ci = lower_ci(smean, se, count),
             upper_ci = upper_ci(smean, se, count),
             n_sample = i) %>%
      select(class, smean, min2, max2, lower_ci, upper_ci, count, n_sample)
  }
  
  #Save all results in a final table
  final = as.data.frame(bind_rows(resS, .id = NULL))
  return(final)
}

boot.pi = function(vcf, range, nboot, replacement =  FALSE) {
  resB = list() #save results by bootstrap
  resS = list() #save results by sample size
  
  sig_test = matrix(NA, length(range), 2)
  colnames(sig_test) = c("n_sample", "p-value")
  final = matrix(NA, length(range), 6)
  colnames(final) = c("mean", "min", "max", "n_sample", "lower_ci", "upper_ci")
  
  PI_obs = Query(vcf, type="site-pi")
  PI_obs = mean(PI_obs$PI, na.rm = T) ## Mean nucleotide divergency per-site
  loop = 0
  
  for(i in range) { #vary sample size
    cat (paste0("Sample size = ", i,"\n"))
    loop=loop+1
    
    for(k in 1:nboot){ #bootstraps
      cat (paste0("Bootstraping... ", k,"\n"))
      positions = sample(c(1:length(vcf@sample_id)), i, replace = replacement) #sample random indiv
      UNIND = vcf@sample_id[positions] #select indiv
      pop_subset = Subset(vcf, samples=UNIND) # subset in a vcf
      PI = Query(pop_subset, type="site-pi") #estimating allelic frequencies
      resB[[k]] = as.data.frame(mean(PI$PI, na.rm = T))
    }
    
    resS[[i]] <- dplyr::bind_rows(resB, .id = NULL) #summarize Bootstraps by a sample size
    #randon test
   # w =  wilcox.test(resS[[i]]$`mean(PI$PI, na.rm = T)`, mu = PI_obs, alternative = "two.sided", exact = FALSE)
    #sig_test[loop,] = c(i, w$p.value)
    w = as.randtest(sim = resS[[i]]$`mean(PI$PI, na.rm = T)`, obs = PI_obs, alter = "two-sided")
    sig_test[loop,] = c(i, w$pvalue)
    
    #save results
    Lower_ci = mean(resS[[i]]$`mean(PI$PI, na.rm = T)`) - qt(1 - ((0.05) / 2), length(resS[[i]]$`mean(PI$PI, na.rm = T)`) - 1) * (sd(resS[[i]]$`mean(PI$PI, na.rm = T)`) / sqrt(length(resS[[i]]$`mean(PI$PI, na.rm = T)`)))
    Upper_ci = mean(resS[[i]]$`mean(PI$PI, na.rm = T)`) + qt(1 - ((0.05) / 2), length(resS[[i]]$`mean(PI$PI, na.rm = T)`) - 1) * (sd(resS[[i]]$`mean(PI$PI, na.rm = T)`) / sqrt(length(resS[[i]]$`mean(PI$PI, na.rm = T)`)))
    
    final[loop,] =cbind(mean(resS[[i]]$`mean(PI$PI, na.rm = T)`), min(resS[[i]]$`mean(PI$PI, na.rm = T)`), max(resS[[i]]$`mean(PI$PI, na.rm = T)`), i[1], Lower_ci, Upper_ci)
    
    
  }
  
  result_final = list(as.data.frame(final), as.data.frame(sig_test))
  names(result_final) = c("Simulations Mean", "Significance Tests")
  cat ("DONE! \n")
  return(result_final)
}

boot.Ne = function(genind, range, nboot, replacement =  FALSE) {
  resB = list() #save results by bootstrap
  resS = list() #save results by sample size
  
  sig_test = matrix(NA, length(range), 2)
  colnames(sig_test) = c("n_sample", "p-value")
  final = matrix(NA, length(range), 6)
  colnames(final) = c("mean", "min", "max", "n_sample", "lower_ci", "upper_ci")
  
  gtypes = genind2gtypes(genind)
  Ne_obs = ldNe(gtypes, maf.threshold = 0, by.strata = TRUE, ci = 0.95, drop.missing = TRUE, num.cores = 3)
  Ne_obs = Ne_obs[, c(6)] # results, using the lower value
  loop = 0
  
  for(i in range) { #vary sample size
    cat (paste0("Sample size = ", i,"\n"))
    loop=loop+1
    
    for(k in 1:nboot){ #bootstraps
      cat (paste0("Bootstraping... ", k,"\n"))
      positions = sample(c(1:length(row.names(genind@tab))), i, replace = replacement) #sample random indiv
      subset = genind[positions,]  #select indiv
      gtypes = genind2gtypes(subset)
      Ne = ldNe(gtypes, maf.threshold = 0, by.strata = TRUE, ci = 0.95, drop.missing = TRUE, num.cores = 3)
      resB[[k]] = as.data.frame(Ne[, c(6)])
    }
    
    resS[[i]] <- dplyr::bind_rows(resB, .id = NULL) #summarize Bootstraps by a sample size
   
     #randon test
    #w =  wilcox.test(resS[[i]]$`Ne[, c(6)]`, mu = Ne_obs, alternative = "two.sided", exact = FALSE)
  #  sig_test[loop,] = c(i, w$p.value)
    w = as.randtest(sim = resS[[i]]$`Ne[, c(6)]`, obs = Ne_obs, alter = "two-sided")
    sig_test[loop,] = c(i, w$pvalue)
    
    #save results
    Lower_ci = mean(resS[[i]]$`Ne[, c(6)]`) - qt(1 - ((0.05) / 2), length(resS[[i]]$`Ne[, c(6)]`) - 1) * (sd(resS[[i]]$`Ne[, c(6)]`) / sqrt(length(resS[[i]]$`Ne[, c(6)]`)))
    Upper_ci = mean(resS[[i]]$`Ne[, c(6)]`) + qt(1 - ((0.05) / 2), length(resS[[i]]$`Ne[, c(6)]`) - 1) * (sd(resS[[i]]$`Ne[, c(6)]`) / sqrt(length(resS[[i]]$`Ne[, c(6)]`)))
    
    final[loop,] =cbind(mean(resS[[i]]$`Ne[, c(6)]`), min(resS[[i]]$`Ne[, c(6)]`), max(resS[[i]]$`Ne[, c(6)]`), i[1], Lower_ci, Upper_ci)
    
    
  }
  
  result_final = list(as.data.frame(final), as.data.frame(sig_test))
  names(result_final) = c("Simulations Mean", "Significance Tests")
  cat ("DONE! \n")
  return(result_final)
}

boot.AR = function(genind, range, nboot, pop, replacement =  FALSE) {
  resB = list() #save results by bootstrap
  resS = list() #save results by sample size
  
  sig_test = matrix(NA, length(range), 2)
  colnames(sig_test) = c("n_sample", "p-value")
  final = matrix(NA, length(range), 6)
  colnames(final) = c("mean", "min", "max", "n_sample", "lower_ci", "upper_ci")
  
  genind@pop = as.factor(pop)
  alle_rich_obs = allel.rich(genind, min.alleles = NULL)
  AR_obs = alle_rich_obs$mean.richness
  loop = 0
  
  for(i in range) { #vary sample size
    cat (paste0("Sample size = ", i,"\n"))
    loop=loop+1
    
    for(k in 1:nboot){ #bootstraps
      cat (paste0("Bootstraping... ", k,"\n"))
      positions = sample(c(1:length(row.names(genind@tab))), i, replace = replacement) #sample random indiv
      subset = genind[positions,]  #select indiv
      alle_rich = allel.rich(subset, min.alleles = NULL)
      resB[[k]] = alle_rich$mean.richness
    }
    
    resS[[i]] <- dplyr::bind_rows(resB, .id = NULL) #summarize Bootstraps by a sample size
    
    
    #randon test
  # w =  wilcox.test(resS[[i]]$`1`, mu = AR_obs, alternative = "two.sided", exact = FALSE)
   # sig_test[loop,] = c(i, w$p.value)
    w = as.randtest(sim = resS[[i]]$`1`, obs = AR_obs, alter = "two-sided")
    sig_test[loop,] = c(i, w$pvalue)
    
    #save results
    Lower_ci = mean(resS[[i]]$`1`) - qt(1 - ((0.05) / 2), length(resS[[i]]$`1`) - 1) * (sd(resS[[i]]$`1`) / sqrt(length(resS[[i]]$`1`)))
    Upper_ci = mean(resS[[i]]$`1`) + qt(1 - ((0.05) / 2), length(resS[[i]]$`1`) - 1) * (sd(resS[[i]]$`1`) / sqrt(length(resS[[i]]$`1`)))
    
    final[loop,] =cbind(mean(resS[[i]]$`1`), min(resS[[i]]$`1`), max(resS[[i]]$`1`), i[1], Lower_ci, Upper_ci)
    
    
  }
  
  result_final = list(as.data.frame(final), as.data.frame(sig_test))
  names(result_final) = c("Simulations Mean", "Significance Tests")
  cat ("DONE! \n")
  return(result_final)
}

boot.alleles = function(genind, range, nboot, copies = 50, replacement = FALSE) {
  resS_N = list() #save results by sample size for null alleles
  resS_R = list() #save results by sample size for rare alleles
  resN = list() #save results by bootstrap for null alleles
  resR = list() #save results by bootstrap for rare alleles
  
  sig_test = matrix(NA, length(range), 2)
  colnames(sig_test) = c("n_sample", "p-value")
  final = matrix(NA, length(range), 6)
  colnames(final) = c("mean", "min", "max", "n_sample", "lower_ci", "upper_ci")
  
  sig_test2 = matrix(NA, length(range), 2)
  colnames(sig_test2) = c("n_sample", "p-value")
  final2 = matrix(NA, length(range), 6)
  colnames(final2) = c("mean", "min", "max", "n_sample", "lower_ci", "upper_ci")
  
  
  total_freq_obs = as.data.frame(colSums(genind@tab, na.rm =T))
  total_freq_obs = total_freq_obs[grep(pattern = '.1$', x = rownames(total_freq_obs), value = T),] #selecting reference allele
  null_allele_obs = length(total_freq_obs[total_freq_obs  == 0]) # alleles not present in population
  rare_allele_obs = length(total_freq_obs[total_freq_obs  < copies]) # rare alleles
  loop = 0
  
  for(i in range) { #vary sample size
    cat (paste0("Sample size = ", i,"\n"))
    loop=loop+1
    
    for(k in 1:nboot){ #bootstraps
      cat (paste0("Bootstraping... ", k,"\n"))
      positions = sample(c(1:length(row.names(genind@tab))), i, replace = replacement) #sample random indiv
      subset = genind[positions,]  #select indiv
      total_freq = as.data.frame(colSums(subset@tab, na.rm =T))
      total_freq = total_freq[grep(pattern = '.1$', x = rownames(total_freq), value = T),]
      resN[[k]] = as.data.frame(length(total_freq[total_freq == 0]))
      resR[[k]] = as.data.frame(length(total_freq[total_freq < copies]))
    }
    
    #NULL ALLELES
    resS_N[[i]] <- dplyr::bind_rows(resN, .id = NULL) #summarize Bootstraps by a sample size
    
    #randon test
    w = as.randtest(sim = resS_N[[i]]$`length(total_freq[total_freq == 0])`, obs = null_allele_obs, alter = "two-sided")
    sig_test[loop,] = c(i, w$pvalue)
  #  w = wilcox.test(resS_N[[i]]$`length(total_freq[total_freq == 0])`, mu = null_allele_obs, alternative = "two.sided", exact = FALSE)
   # sig_test[loop,] = c(i, w$p.value)
    
    #save results
    Lower_ci = mean(resS_N[[i]]$`length(total_freq[total_freq == 0])`) - qt(1 - ((0.05) / 2), length(resS_N[[i]]$`length(total_freq[total_freq == 0])`) - 1) * (sd(resS_N[[i]]$`length(total_freq[total_freq == 0])`) / sqrt(length(resS_N[[i]]$`length(total_freq[total_freq == 0])`)))
    Upper_ci = mean(resS_N[[i]]$`length(total_freq[total_freq == 0])`) + qt(1 - ((0.05) / 2), length(resS_N[[i]]$`length(total_freq[total_freq == 0])`) - 1) * (sd(resS_N[[i]]$`length(total_freq[total_freq == 0])`) / sqrt(length(resS_N[[i]]$`length(total_freq[total_freq == 0])`)))
    
    final[loop,] =cbind(mean(resS_N[[i]]$`length(total_freq[total_freq == 0])`), min(resS_N[[i]]$`length(total_freq[total_freq == 0])`), max(resS_N[[i]]$`length(total_freq[total_freq == 0])`), i[1], Lower_ci, Upper_ci)
    
    
    
    #RARE ALLELES
    resS_R[[i]] <- dplyr::bind_rows(resR, .id = NULL) #summarize Bootstraps by a sample size
    
    #randon test
    w2 = as.randtest(sim = resS_R[[i]]$`length(total_freq[total_freq < copies])`, obs = rare_allele_obs, alter = "two-sided")
    sig_test2[loop,] = c(i, w$pvalue)
     #w2 = wilcox.test(resS_R[[i]]$`length(total_freq[total_freq < copies])`, mu = rare_allele_obs, alternative = "two.sided", exact = FALSE)
    #sig_test2[loop,] = c(i, w2$p.value)
    
    #save results
    Lower_ci2 = mean(resS_R[[i]]$`length(total_freq[total_freq < copies])`) - qt(1 - ((0.05) / 2), length(resS_R[[i]]$`length(total_freq[total_freq < copies])`) - 1) * (sd(resS_R[[i]]$`length(total_freq[total_freq < copies])`) / sqrt(length(resS_R[[i]]$`length(total_freq[total_freq < copies])`)))
    Upper_ci2 = mean(resS_R[[i]]$`length(total_freq[total_freq < copies])`) + qt(1 - ((0.05) / 2), length(resS_R[[i]]$`length(total_freq[total_freq < copies])`) - 1) * (sd(resS_R[[i]]$`length(total_freq[total_freq < copies])`) / sqrt(length(resS_R[[i]]$`length(total_freq[total_freq < copies])`)))
    
    final2[loop,] =cbind(mean(resS_R[[i]]$`length(total_freq[total_freq < copies])`), min(resS_R[[i]]$`length(total_freq[total_freq < copies])`), max(resS_R[[i]]$`length(total_freq[total_freq < copies])`), i[1], Lower_ci2, Upper_ci2)
    
    
  }
  
  result_final = list(as.data.frame(final), as.data.frame(sig_test), as.data.frame(final2), as.data.frame(sig_test2))
  names(result_final) = c("Simulations Mean Null Alleles", "Significance Tests Null Alleles", "Simulations Mean Rare Alleles", "Significance Tests Rare Alleles")
  cat ("DONE! \n")
   return(result_final)
}

boot.he = function(genind, range, nboot, replacement =  FALSE) {
  resB = list() #save results by bootstrap
  resS = list() #save results by sample size
  sig_test = matrix(NA, length(range), 2)
  colnames(sig_test) = c("n_sample", "p-value")
  final = matrix(NA, length(range), 6)
  colnames(final) = c("mean", "min", "max", "n_sample", "lower_ci", "upper_ci")
  he_obs = Hs(genind)
  loop = 0
  
  for(i in range) { #vary sample size
    cat (paste0("Sample size = ", i,"\n"))
    loop=loop+1
    
    for(k in 1:nboot){ #bootstraps
      cat (paste0("Bootstraping... ", k,"\n"))
      positions = sample(c(1:length(row.names(genind@tab))), i, replace = replacement) #sample random indiv
      subset = genind[positions,]  #select indiv
      resB[[k]] = Hs(subset, pop = NULL)
    }
    
    resS[[i]] <- dplyr::bind_rows(resB, .id = NULL) #summarize Bootstraps by a sample size
    #randon test
    w = as.randtest(sim = resS[[i]]$`1`, obs = he_obs, alter = "two-sided")
    sig_test[loop,] = c(i, w$pvalue)
   # w =  wilcox.test(resS[[i]]$`1`, mu = he_obs, alternative = "two.sided", exact = FALSE)
    #sig_test[loop,] = c(i, w$p.value)
    
    #save results
    Lower_ci = mean(resS[[i]]$`1`) - qt(1 - ((0.05) / 2), length(resS[[i]]$`1`) - 1) * (sd(resS[[i]]$`1`) / sqrt(length(resS[[i]]$`1`)))
    Upper_ci = mean(resS[[i]]$`1`) + qt(1 - ((0.05) / 2), length(resS[[i]]$`1`) - 1) * (sd(resS[[i]]$`1`) / sqrt(length(resS[[i]]$`1`)))
    
    final[loop,] =cbind(mean(resS[[i]]$`1`), min(resS[[i]]$`1`), max(resS[[i]]$`1`), i[1], Lower_ci, Upper_ci)
    
    
  }
  
  result_final = list(as.data.frame(final), as.data.frame(sig_test))
  names(result_final) = c("Simulations Mean", "Significance Tests")
  cat ("DONE! \n")
  return(result_final)
}
