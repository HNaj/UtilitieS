# =============================
# returns matrix with QC measures; argument row.name needs to be specified before executing the function

do.QC <- function(df, row.name = row.name) {
  
  if(!require(dplyr)) {
    install.packages("dplyr")
    library(dplyr)
  }
  
  dimnames = list(paste(row.name),c("Var",
                                    "Het",
                                    "Hom",
                                    "Het/Hom",
                                    "SNVs",
                                    "Insertions",
                                    "Deletions",
                                    "Indels",
                                    "Transition",
                                    "Transvertion",
                                    "Ti/Tv",
                                    "FPR",
                                    "NS",
                                    "SS",
                                    "NS/SS"))
  
  QC_after_subsetting <- matrix(NA, 1 , 15, dimnames = dimnames)
  
  # remove duplicates
  variants = df[!duplicated(df$Variant), ]
  
  # total number of coding variants
  QC_after_subsetting[1, 1] <- nrow(variants)
  
  # number of coding heterozygous
  het.var = variants %>% filter(Genotype %in% c("0/1", "0/2", "0/3", "0/4", "0/5", "0/6", "1/2", "1/3", "1/4", "1/5", "1/6", "2/3", "2/4", "2/5", "3/4", "3/6", "4/5", "5/6"))
  QC_after_subsetting[1, 2] <- nrow(het.var)
  
  # number of coding homozygous
  homo.var = variants %>%  filter(Genotype %in% c("1/1", "2/2", "3/3", "4/4", "5/5"))
  QC_after_subsetting[1, 3] <- nrow(homo.var)
  
  # heterozygous/homozygous ratio
  ratio <- nrow(het.var) / nrow(homo.var)
  QC_after_subsetting[1, 4] <- round(ratio, digits = 2)
  
  # number of coding SNVs
  SNV <-
    filter(variants,
           ALT %in% c("A", "C", "G", "T"),
           REF %in% c("A", "C", "G", "T"))
  QC_after_subsetting[1, 5] <- nrow(SNV)
  
  # number of insertions
  Insertions = filter(variants,
                      ALT != "A",
                      ALT != "C",
                      ALT != "G",
                      ALT != "T",
                      REF %in% c("A", "C", "G", "T"))
  QC_after_subsetting[1, 6] <- nrow(Insertions)
  
  # number of deletions
  Deletions = filter(variants,
                     ALT %in% c("A", "C", "G", "T"),
                     REF != "A",
                     REF != "C",
                     REF != "G",
                     REF != "T")
  QC_after_subsetting[1, 7] <- nrow(Deletions)
  
  # number of indels
  QC_after_subsetting[1, 8] <- nrow(Insertions) + nrow(Deletions)
  
  ## Number of transitions
  # Two-ring purines (A <-> G)
  var.Transition.AG = filter(SNV, ALT %in% c("A", "G"), REF %in% c("A", "G"))
  nrow(var.Transition.AG)
  # One-ring pyrimidines (C <-> T)
  var.Transition.CT = filter(SNV, ALT %in% c("C", "T"), REF %in% c("C", "T"))
  nrow(var.Transition.CT)
  # Total number of transitions per sample
  var.Transition <-
    ((nrow(var.Transition.AG)) + (nrow(var.Transition.CT)))
  QC_after_subsetting[1, 9] <- var.Transition
  
  ## Number of transversions
  var.Transver1 = filter(SNV, REF %in% c("C", "T"), ALT %in% c("A", "G"))
  nrow(var.Transver1)
  var.Transver2 = filter(SNV, REF %in% c("A", "G"), ALT %in% c("C", "T"))
  nrow(var.Transver2)
  var.Transver <- ((nrow(var.Transver1)) + (nrow(var.Transver2)))
  QC_after_subsetting[1, 10] <- var.Transver
  
  ## Ratio Tsi/Tsver
  QC_after_subsetting[1, 11] <-
    round(var.Transition / var.Transver, digits = 2)
  
  ## An estimate of the false positive rate
  FPR = (1 - ((var.Transition / var.Transver) - 0.5) / (3.3 - 0.5))
  QC_after_subsetting[1, 12] <- round(FPR, digits = 2)
  
  ## Non-synonymous SNVs
  nonsynon.SNV = filter(SNV, Consequence %in% c("stop_gained", "missense_variant"))
  QC_after_subsetting[1, 13] <- nrow(nonsynon.SNV)
  
  ## Synonymous SNVs
  synon.SNV = filter(SNV, Consequence == "synonymous_variant")
  QC_after_subsetting[1, 14] <- nrow(synon.SNV)
  
  ## Ratio Ns/S SNVs
  nS.S_SNVs_ratio  = (nrow(nonsynon.SNV) / nrow(synon.SNV))
  QC_after_subsetting[1, 15] <- round(nS.S_SNVs_ratio, digits = 2)
  
  return(QC_after_subsetting)
}
