# =============================
# performs permutations - doPermutations() calls create.random.distribution()- creates random distribution for given population - which calls sample.data() which performs sampling

do.Permutations <- function(data=observed_data, sampleSize=sample.size,nPerms=2, biallelic=FALSE){
  
  if(!require(dplyr)) {
    install.packages("dplyr")
    library(dplyr)
  }
  
  data$count=0
  # results are generated here
  for(i in 1:nPerms){
    expected.nVar.by.gene <- create.random.distribution(sampleSize,
                                                        biallelic=biallelic)
    
    data <- data %>% 
      tbl_df %>%
      left_join(expected.nVar.by.gene, by=c("Feature"="Gene")) %>% 
      removeNAs %>%
      mutate(test=Expected.nVar>=nSampleWithVar) %>%
      transmute(Feature,nSampleWithVar, count=count+test)
  }
  
  data %>% 
    transmute(Feature,
              observed=nSampleWithVar,
              pVal=count/nPerms) %>%
    arrange(pVal,desc(observed))
  
}

# =============================
# performs sampling - sample.data() is called by create.random.distribution()
# input is dataframe with gene names (df$ensembl_transcript_id) and probability or gene lenght - (df$cds_length)
# number of samples to take is defined by size, i.e. it corresponds to the number of variants from each patient

sample.data <- function(size, df = df_scores) {
  data.frame(
    Gene = sample(
      x = df$ensembl_transcript_id,
      size = size,
      replace = T,
      prob = df$cds_length),
    stringsAsFactors = F) %>%
    group_by(Gene) %>%
    summarize(Freq = n()) %>%
    return
}

# =============================
# creates random distribution for given population - create.random.distribution() calls sample.data() which performs sampling

create.random.distribution <- function(sample.size, biallelic = FALSE) {
  # sample.size is a vector, for each sample contains n of variants in that sample
  number.of.patients <- length(sample.size)
  # create a list into which loop results will go
  list.all <- list()
  # repeat the sampling for each patient
  # sample size for each iteration corresponds to the number of variants for a given patient
  # the list.all output will contain one random distribution for all patients
  for(i in 1:number.of.patients) {
    sampling <- sample.data(size = sample.size[i])
    list.all[[i]] <- sampling
  }
  
  # convert the list to df
  list2dataframe <- do.call(rbind.data.frame, list.all)
  # if biallelic = TRUE, only count genes with 2 hits, ie recessive model
  if (biallelic) {
    list2dataframe <- filter(list2dataframe, Freq >= 2)
  }
  # count how many times a given gene appears in the random population
  list2dataframe %>%
    group_by(Gene) %>%
    summarise(Expected.nVar = n()) %>%
    return
}

# =============================
# replaces NAs with 0 - removeNAs() is called by doPermutations(), which also calls create.random.distribution() - creates random distribution for given population - which calls sample.data() which performs sampling

removeNAs <- function(x){
  x[is.na(x)] <- 0;return(x)
}

