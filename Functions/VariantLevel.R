
if(!require(dplyr)) {
    install.packages("dplyr")
    library(dplyr)
}

# =============================
# returns subset of data frame that contains PROTEIN ALTERING variants (missense, nonsense, indel, frameshift, splice)

select.protein.altering <- function(df) {
  df %>%
    filter(Consequence %in% c(
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "missense_variant",
        "protein_altering_variant",
        "start_lost",
        "stop_lost")) %>%
    return()
}

# =============================
# returns count of PROTEIN ALTERING variants in distinct 'Consequence' categories (missense, nonsense, indel, frameshift, splice)

count.protein.altering <- function(df) {
  protein_altering_count = df[!duplicated(df$Variant), ]
  protein_altering_count  %>%
    filter(Consequence %in% c(
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "missense_variant",
        "protein_altering_variant",
        "start_lost",
        "stop_lost")) %>%
    group_by(Consequence) %>%
    summarise(total = n()) %>%
    return()
}

# =============================
# returns subset of data frame that contains LOSS OF FUNCTION variants (nonsense, splice, frameshift)

select.LoF <- function(df) {
  df  %>%
    filter(Consequence == "frameshift_variant" |
        Consequence == "stop_gained" |
        Consequence == "transcript_ablation" |
        Consequence == "splice_acceptor_variant" |
        Consequence == "splice_donor_variant" |
        Consequence == "start_lost" | Consequence == "stop_lost") %>%
    return()
}

# =============================
# returns count of LoF variants in distinct 'Consequence' categories

count.LoF <- function(df) {
  LoF_count = df[!duplicated(df$Variant), ]
  LoF_count  %>%
    filter(Consequence == "frameshift_variant" |
        Consequence == "stop_gained" |
        Consequence == "transcript_ablation" |
        Consequence == "splice_acceptor_variant" |
        Consequence == "splice_donor_variant" |
        Consequence == "start_lost" | Consequence == "stop_lost") %>%
    group_by(Consequence) %>%
    summarise(total = n()) %>%
    return()
}

# =============================
# returns non-cumulative number of variants (bottom row) in given number of exomes (top row)

count.variant = function(df) {
  df_unique = df[!duplicated(df$Variant), ]
  return(table(df_unique$SampleCount))
}

