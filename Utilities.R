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
      positions = sample(c(1:length(y@sample_id)), i, replace = FALSE) #sample random indiv
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

boot.he = function(genind, range, nboot) {
  resB = list() #save results by bootstrap
  resS = list() #save results by sample size
  sig_test = matrix(NA, length(range), 2)
  colnames(sig_test) = c("n_sample", "p-value")
  final = matrix(NA, length(range), 6)
  colnames(final) = c("mean", "min", "max", "n_sample", "lower_ci", "upper_ci")
  he_obs = Hs(genind)
  loop = 0
  
  for(i in range) { #vary sample size
    loop=loop+1
    
    for(k in 1:nboot){ #bootstraps
      positions = sample(c(1:length(row.names(genind@tab))), i, replace = TRUE) #sample random indiv
      subset = genind[positions,]  #select indiv
      resB[[k]] = Hs(subset, pop = NULL)
    }
    
    resS[[i]] <- dplyr::bind_rows(resB, .id = NULL) #summarize Bootstraps by a sample size
    #randon test
    w = as.randtest(sim = resS[[i]]$`1`, obs = he_obs, alter = "two-sided")
    sig_test[loop,] = c(i, w$pvalue)
    
    #save results
    Lower_ci = mean(resS[[i]]$`1`) - qt(1 - ((0.05) / 2), length(resS[[i]]$`1`) - 1) * (sd(resS[[i]]$`1`) / sqrt(length(resS[[i]]$`1`)))
    Upper_ci = mean(resS[[i]]$`1`) + qt(1 - ((0.05) / 2), length(resS[[i]]$`1`) - 1) * (sd(resS[[i]]$`1`) / sqrt(length(resS[[i]]$`1`)))
    
    final[loop,] =cbind(mean(resS[[i]]$`1`), min(resS[[i]]$`1`), max(resS[[i]]$`1`), i[1], Lower_ci, Upper_ci)
    
    
  }
  
  result_final = list(as.data.frame(final), as.data.frame(sig_test))
  names(result_final) = c("Simulations Mean", "Significance Tests")
  return(result_final)
}




