# =============================
# returns number of genes common to any subset of X exomes; argument list.of.samples needs to be specified before executing the function

count.gene <- function(df, list.of.samples) {
  list.of.samples <- list.of.samples
  
  nVarPerGene_perSample <- data.frame(Feature = c("xxx", "xxx"))
  
  for (i in 1:length(list.of.samples)) {
    # remove duplicates & restrict to canonical transcripts
    nVar_PerGene =
      df %>%
      filter(!duplicated(IDENTIFIER1)) %>%
      filter(CANONICAL == "YES")
    
    # count variants per gene in individual samples
    PerSampVar =
      nVar_PerGene %>%
      filter(samples.id == list.of.samples[i])
    
    # group variants per gene and count for each sample
    PerSampGene =
      PerSampVar %>%
      group_by(Feature) %>%
      summarise(nVarPerGene = n())
    
    colnames(PerSampGene)[2] <- paste(list.of.samples[i])
    
    # populate the 'nVarPerGene_perSample' table
    nVarPerGene_perSample <-
      full_join(nVarPerGene_perSample,
                PerSampGene,
                by = "Feature",
                copy = TRUE)
  }
  
  nVarPerGene_perSample = nVarPerGene_perSample[-c(1:2),] # removes 2 first rows
  
  # count NA to work out how many exomes have at least one novel variant in a given gene
  logical_nVarPerGene_perSample = as.data.frame(is.na(nVarPerGene_perSample[,2:ncol(nVarPerGene_perSample)]))
  
  logical_nVarPerGene_perSample$Count.NA <-
    rowSums(logical_nVarPerGene_perSample , na.rm = TRUE)
  
  # Count nExome with at least one variant (or more) in a given gene
  nVarPerGene_perSample_Count = mutate(logical_nVarPerGene_perSample,
                                       nSampleWithVar = (length(list.of.samples) - Count.NA))
  
  nVarPerGene_perSample_Count %>%
    group_by(nSampleWithVar) %>%
    summarise(nGene = n()) %>%
    return()
}

# =============================
# defines function that plots number of genes common to any subset of X affected exomes

plot.nGene <- function(x, main) {
  nGene_nPeople = x %>% select(Feature, nSampleWithVar)
  barplot(
    table(nGene_nPeople$nSampleWithVar), 
    xlab = "Any nPeople",
    ylab = "nGene",
    col = "red",
    main = main)
  return()
}
