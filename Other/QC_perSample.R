# QC per Sample

library(dplyr)

protein_coding = df %>% filter(
  FILTER == "PASS",
  Consequence %in% c("frameshift_variant",
                     "stop_gained",
                     "splice_donor_variant",
                     "missense_variant",
                     "splice_acceptor_variant",
                     "protein_altering_variant",
                     "start_lost",
                     "stop_lost",
                     "inframe_deletion",
                     "inframe_insertion",
                     "transcript_ablation",
                     "splice_region_variant",
                     "stop_retained_variant",
                     "synonymous_variant",
                     "coding_sequence_variant"),
  BIOTYPE == "protein_coding")

# Remove duplicats
Protein_Coding_Var = protein_coding[!duplicated(protein_coding$Variant), ]

# Restrict to novel protein-coding variants
Novel_Coding_Var <-
  Protein_Coding_Var[is.na(Protein_Coding_Var$ExAC_MAF),]

# Define list of samples
list.of.samples <- c("sample1", "sample2", "sample3", "sample4")

dimnames = list(paste(list.of.samples),c("var",
                                         "var_novel",
                                         "%novel",
                                         "Het",
                                         "Hom",
                                         "Het/Hom",
                                         "Het_novel",
                                         "Hom_novel",
                                         "Het/Hom_novel",
                                         "SNV",
                                         "SNV_novel",
                                         "%SNV_novel",
                                         "Insertion",
                                         "Deletion",
                                         "Indel",
                                         "Inser_novel",
                                         "Del_novel",
                                         "Indel novel",
                                         "%Indel_novel",
                                         "Transition",
                                         "Transition_novel",
                                         "%novel",
                                         "Transvertion",
                                         "Transvertion_novel",
                                         "%novel",
                                         "Ratio_Tsi/Tsver",
                                         "FPR",
                                         "Ratio_novel_Tsi/Tsver",
                                         "FPR_novel",
                                         "Ns/S_SNVs_ratio",
                                         "Ns/S_novelSNVs_ratio",
                                         "Ns/S_ratio",
                                         "Ns/S_novelVar_ratio"))

result <-
  matrix(NA, length(list.of.samples), 33, dimnames = dimnames)

for (i in 1:length(list.of.samples)) {
  ## Total number of coding variants per sample
  var.per.sample = Protein_Coding_Var[Protein_Coding_Var[paste(list.of.samples[i])] != "./." &
                                        Protein_Coding_Var[paste(list.of.samples[i])] != "0/0",]
  result[i, 1] <- nrow(var.per.sample)
  
  ## Number of novel coding variants per sample
  novel.per.sample <-
    Novel_Coding_Var[Novel_Coding_Var[paste(list.of.samples[i])] != "./." &
                       Novel_Coding_Var[paste(list.of.samples[i])] != "0/0",]
  result[i, 2] <- nrow(novel.per.sample)
  
  ## % of novel coding variants per sample
  procent.novel <-
    (100 * (nrow(novel.per.sample)) / (nrow(var.per.sample)))
  result[i, 3] <- round(procent.novel, digits = 2)
  
  ## Number of coding heterozygous per sample
  het.var.per.sample <-
    var.per.sample[var.per.sample[paste(list.of.samples[i])] == "0/1" |
                     var.per.sample[paste(list.of.samples[i])] == "0/2" |
                     var.per.sample[paste(list.of.samples[i])] == "0/3" |
                     var.per.sample[paste(list.of.samples[i])] == "0/4" |
                     var.per.sample[paste(list.of.samples[i])] == "1/2" |
                     var.per.sample[paste(list.of.samples[i])] == "1/3" |
                     var.per.sample[paste(list.of.samples[i])] == "1/4" |
                     var.per.sample[paste(list.of.samples[i])] == "2/3" |
                     var.per.sample[paste(list.of.samples[i])] == "2/4" |
                     var.per.sample[paste(list.of.samples[i])] == "3/4",]
  
  result[i, 4] <- nrow(het.var.per.sample)
  
  ## Number of coding homozygous calls per sample
  homo.var.per.sample <-
    var.per.sample[var.per.sample[paste(list.of.samples[i])] == "1/1" |
                     var.per.sample[paste(list.of.samples[i])] == "2/2" |
                     var.per.sample[paste(list.of.samples[i])] == "3/3" |
                     var.per.sample[paste(list.of.samples[i])] == "4/4",]
  
  result[i, 5] <- nrow(homo.var.per.sample)
  
  ## Heterozygous/homozygous ratio per sample
  ratio <- nrow(het.var.per.sample) / nrow(homo.var.per.sample)
  result[i, 6] <- round(ratio, digits = 2)
  
  ## Number of 'novel' coding heterozygous per sample
  het.novel.per.sample <-
    novel.per.sample[novel.per.sample[paste(list.of.samples[i])] == "0/1" |
                       novel.per.sample[paste(list.of.samples[i])] == "0/2" |
                       novel.per.sample[paste(list.of.samples[i])] == "0/3" |
                       novel.per.sample[paste(list.of.samples[i])] == "0/4" |
                       novel.per.sample[paste(list.of.samples[i])] == "1/2" |
                       novel.per.sample[paste(list.of.samples[i])] == "1/3" |
                       novel.per.sample[paste(list.of.samples[i])] == "1/4" |
                       novel.per.sample[paste(list.of.samples[i])] == "2/3" |
                       novel.per.sample[paste(list.of.samples[i])] == "2/4" |
                       novel.per.sample[paste(list.of.samples[i])] == "3/4",]
  
  result[i, 7] <-  nrow(het.novel.per.sample)
  
  ## Number of 'novel' coding homozygous per sample
  homo.novel.per.sample = novel.per.sample[novel.per.sample[paste(list.of.samples[i])] == "1/1" |
                                             novel.per.sample[paste(list.of.samples[i])] == "2/2" |
                                             novel.per.sample[paste(list.of.samples[i])] == "3/3" |
                                             novel.per.sample[paste(list.of.samples[i])] == "4/4",]
  
  result[i, 8] <- nrow(homo.novel.per.sample)
  
  ## 'Novel' heterozygous/homozygous ratio per sample
  ratio.novel <-
    (nrow(het.novel.per.sample) / nrow(homo.novel.per.sample))
  
  result[i, 9] <- round(ratio.novel, digits = 2)
  
  ## Number of coding SNVs per sample
  SNV.per.sample <-
    filter(var.per.sample,
           ALT %in% c("A", "C", "G", "T"),
           REF %in% c("A", "C", "G", "T"))
  
  result[i, 10] <- nrow(SNV.per.sample)
  
  ## Number of novel coding SNVs per sample
  novel.SNV.per.sample = filter(novel.per.sample,
                                ALT %in% c("A", "C", "G", "T"),
                                REF %in% c("A", "C", "G", "T"))
  
  result[i, 11] <- nrow(novel.SNV.per.sample)
  
  ## % Novel SNVs per sample
  proc.SNV.novel = (100 * (nrow(novel.SNV.per.sample)) / nrow(SNV.per.sample))
  
  result[i, 12] <- round(proc.SNV.novel, digits = 2)
  
  ## Number of insertions per sample
  Inser.per.sample = filter(var.per.sample,
                            ALT != "A",
                            ALT != "C",
                            ALT != "G",
                            ALT != "T",
                            REF %in% c("A", "C", "G", "T"))
  
  result[i, 13] <- nrow(Inser.per.sample)
  
  ## Number of deletions per sample
  Del.per.sample = filter(var.per.sample,
                          ALT %in% c("A", "C", "G", "T"),
                          REF != "A",
                          REF != "C",
                          REF != "G",
                          REF != "T")
  
  result[i, 14] <- nrow(Del.per.sample)
  
  ## Number of indels per sample
  result[i, 15] <- nrow(Inser.per.sample) + nrow(Del.per.sample)
  
  ## Number of novel insertions per sample
  novel.Inser.per.sample = filter(
    novel.per.sample,
    ALT != "A",
    ALT != "C",
    ALT != "G",
    ALT != "T",
    REF %in% c("A", "C", "G", "T"))
  
  result[i, 16] <- nrow(novel.Inser.per.sample)
  
  ## Number of novel deletions per sample
  novel.Del.per.sample = filter(
    novel.per.sample,
    ALT %in% c("A", "C", "G", "T"),
    REF != "A",
    REF != "C",
    REF != "G",
    REF != "T")
  
  result[i, 17] <- nrow(novel.Del.per.sample)
  
  ## Number of novel indels per sample
  result[i, 18] <-
    nrow(novel.Inser.per.sample) + nrow(novel.Del.per.sample)
  
  ## % Novel Indels per sample
  proc.indels = (100 * (nrow(novel.Inser.per.sample) + nrow(novel.Del.per.sample)) /
                   (nrow(Inser.per.sample) + nrow(Del.per.sample)))
  
  result[i, 19] <- round(proc.indels, digits = 2)
  
  ## Number of transitions per sample
  # Two-ring purines (A <-> G)
  var.Transition.AG = filter(SNV.per.sample, ALT %in% c("A", "G"), REF %in% c("A", "G"))
  nrow(var.Transition.AG)
  # One-ring pyrimidines (C <-> T)
  var.Transition.CT = filter(SNV.per.sample, ALT %in% c("C", "T"), REF %in% c("C", "T"))
  nrow(var.Transition.CT)
  # Total number of transitions per sample
  var.Transition <-
    ((nrow(var.Transition.AG)) + (nrow(var.Transition.CT)))
  
  result[i, 20] <- var.Transition
  
  ## Number of novel transitions per sample
  # Two-ring purines (A <-> G)
  novel.Transition.AG = filter(novel.SNV.per.sample, ALT %in% c("A", "G"), REF %in% c("A", "G"))
  nrow(novel.Transition.AG)
  # One-ring pyrimidines (C <-> T)
  novel.Transition.CT = filter(novel.SNV.per.sample, ALT %in% c("C", "T"), REF %in% c("C", "T"))
  nrow(novel.Transition.CT)
  # Total number of novel transitions per sample
  novel.Transition <-
    ((nrow(novel.Transition.AG)) + (nrow(novel.Transition.CT)))
  
  result[i, 21] <- novel.Transition
  
  ## % Novel transitions per sample
  proc.transitions = (100 * novel.Transition / var.Transition)
  
  result[i, 22] <- round(proc.transitions, digits = 2)
  
  ## Number of transversions per sample
  var.Transver1 = filter(SNV.per.sample, REF %in% c("C", "T"), ALT %in% c("A", "G"))
  nrow(var.Transver1)
  var.Transver2 = filter(SNV.per.sample, REF %in% c("A", "G"), ALT %in% c("C", "T"))
  nrow(var.Transver2)
  var.Transver <- ((nrow(var.Transver1)) + (nrow(var.Transver2)))
  
  result[i, 23] <- var.Transver
  
  #### Number of novel transversions per sample
  novel.Transver1 = filter(novel.SNV.per.sample, REF %in% c("C", "T"), ALT %in% c("A", "G"))
  nrow(novel.Transver1)
  novel.Transver2 = filter(novel.SNV.per.sample, REF %in% c("A", "G"), ALT %in% c("C", "T"))
  nrow(novel.Transver2)
  novel.Transver <-
    ((nrow(novel.Transver1)) + (nrow(novel.Transver2)))
  
  result[i, 24] <- novel.Transver
  
  ## % Novel transversions per sample
  proc.transversions = (100 * novel.Transver / var.Transver)
  
  result[i, 25] <- round(proc.transversions, digits = 2)
  
  ## Ratio Tsi/Tsver
  result[i, 26] <- round(var.Transition / var.Transver, digits = 2)
  
  ## An estimate of the false positive rate
  FPR = (1 - ((var.Transition / var.Transver) - 0.5) / (3.3 - 0.5))
  result[i, 27] <- round(FPR, digits = 2)
  
  ## Ratio novel Tsi/Tsver
  result[i, 28] <- round(novel.Transition / novel.Transver, digits = 2)
  
  ## An estimate of the false positive rate
  FPR.novel = (1 - ((novel.Transition / novel.Transver) - 0.5) / (3.3 -
                                                                    0.5))
  result[i, 29] <- round(FPR.novel, digits = 2)
  
  ## Synonymous SNVs
  synon.SNV.per.sample = filter(SNV.per.sample, Consequence == "synonymous_variant")
  nrow(synon.SNV.per.sample)
  
  ## Non-synonymous SNVs
  nonsynon.SNV.per.sample = filter(SNV.per.sample,
                                   Consequence %in% c("stop_gained", "missense_variant"))
  nrow(nonsynon.SNV.per.sample)
  
  ## Ratio Ns/S SNVs
  nS.S_SNVs_ratio  = (nrow(nonsynon.SNV.per.sample) / nrow(synon.SNV.per.sample))
  
  result[i, 30] <- round(nS.S_SNVs_ratio, digits = 2)
  
  ## Synonymous novel SNVs
  synon.SNVnovel.per.sample = filter(novel.SNV.per.sample, Consequence == "synonymous_variant")
  nrow(synon.SNVnovel.per.sample)
  
  ## Non-synonymous novel SNVs
  nonsyn.SNVnovel.per.sample = filter(novel.SNV.per.sample,
                                      Consequence %in% c("stop_gained", "missense_variant"))
  nrow(nonsyn.SNVnovel.per.sample)
  
  ## Ratio Ns/S novel SNVs
  nS.S_novelSNVs_ratio  = (nrow(nonsyn.SNVnovel.per.sample) / nrow(synon.SNVnovel.per.sample))
  
  result[i, 31] <- round(nS.S_novelSNVs_ratio, digits = 2)
  
  ## Synonymous variants
  synonymous.per.sample = filter(var.per.sample, Consequence == "synonymous_variant")
  nrow(synonymous.per.sample)
  
  ## Non-synonymous variants
  nonsynonymous.per.sample = filter(var.per.sample,
                                    Consequence %in% c("stop_gained", "missense_variant"))
  nrow(nonsynonymous.per.sample)
  
  ## Ratio Ns/S variants
  nS.S_var_ratio  = (nrow(nonsynonymous.per.sample) / nrow(synonymous.per.sample))
  
  result[i, 32] <- round(nS.S_var_ratio, digits = 2)
  
  ## Synonymous novel variants
  synon.novel.per.sample = filter(novel.per.sample, Consequence == "synonymous_variant")
  nrow(synon.novel.per.sample)
  
  ## Non-synonymous novel variants
  nonsyn.novel.per.sample = filter(novel.per.sample,
                                   Consequence %in% c("stop_gained", "missense_variant"))
  nrow(nonsyn.novel.per.sample)
  
  ## Ratio Ns/S novel variants
  nS.S_novel.var_ratio  = (nrow(nonsyn.novel.per.sample) / nrow(synon.novel.per.sample))
  
  result[i, 33] <- round(nS.S_novel.var_ratio, digits = 2)
}

result
