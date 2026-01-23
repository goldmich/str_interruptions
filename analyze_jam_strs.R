# Beginning to analyze the Jam et al. 1kgp & h3a data

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Documents/Quinlan/Interruptions/")

jam_genotypes_enhancers <- read.delim("./interruptions_full_annotations.txt", header = F)
colnames(jam_genotypes_enhancers) <-
  c("chrom", "start", "end", "ru", "len", "max_homo_ru", "levenshtein_distance", "maf", 
    "in_vista", "gh_element", "elite", "min_elite_loeuf", "genehancer_gene_is_ad", "overlap_coding_genes", 
    "overlap_coding_transcripts", "coding_gene_loeuf", "coding_gene_is_ad", "overlap_tf_chip",
    "cpgi", "fire", "ipsc_atac", "k562_gm12878_atac", "problematic_region", "vista",
    "chromhmm_numbers", "chromhmm_state")

tf_list <-
  c("BRCA1", "CTCF", "E2F1", "E2F4", "ELK1", "MAX", "MYC", "NFE2L2", "PRDM1", "SMARCA4")

# write.table(
#   jam_genotypes_enhancers %>%
#     group_by(., chrom, start, end, in_genehancer, elite, min_elite_loeuf) %>%
#     filter(., maf == max(maf) & (levenshtein_distance / len) < 0.85 & min_elite_loeuf != -1),
#   "./constrained_strs_for_hope.txt", sep = "\t", row.names = F
# )

jam_genotypes_enhancers$"fraction_gc" <- 
  sapply(jam_genotypes_enhancers$ru,
         function(x)
           nchar(gsub("[AT]+", "", x)) / 
           nchar(x))
jam_genotypes_enhancers$"in_genehancer" <- jam_genotypes_enhancers$gh_element != "."
jam_genotypes_enhancers$"max_len_homology" <-
  jam_genotypes_enhancers$max_homo_ru * nchar(jam_genotypes_enhancers$ru)
jam_genotypes_enhancers$"purity" <-
  1 - (jam_genotypes_enhancers$levenshtein_distance / jam_genotypes_enhancers$len)
jam_genotypes_enhancers$elite <-
  sapply(
    jam_genotypes_enhancers$elite,
    function(x)
      if(x == "0,1") "1" else x
  )
jam_genotypes_enhancers$"overlap_coding" <- jam_genotypes_enhancers$overlap_coding_transcripts != "."
jam_genotypes_enhancers$overlap_coding_transcripts <- NULL
jam_genotypes_enhancers$coding_gene_loeuf <- 
  sapply(
    jam_genotypes_enhancers$coding_gene_loeuf,
    function(x)
      if(x == -1) NA else x
  )
jam_genotypes_enhancers$"contains_cpg" <-
  sapply(
    jam_genotypes_enhancers$ru,
    function(x)
      grepl("CG", paste0(x, x))
  )

jam_genotypes_enhancers$"n_unique_nucleotides" <-
  sapply(
    jam_genotypes_enhancers$ru,
    function(x)
      length(unique(strsplit(x, "")[[1]]))
  )
jam_genotypes_enhancers <-
  subset(jam_genotypes_enhancers,
         jam_genotypes_enhancers$n_unique_nucleotides > 1)
jam_genotypes_enhancers$n_unique_nucleotides <- NULL

jam_genotypes_enhancers <- subset(
  jam_genotypes_enhancers,
  jam_genotypes_enhancers$len >= 10
)

# detect_cpg <- function(ru) {
#   grepl("CG", paste0(c(ru, ru)))
# }


# for loci where the major allele is perfect, are loci with higher MAF interrupted alleles more likely to fall within enhancers?

# relative to the major allele, are minor alleles also more perfect at constrained sites?

# summary(
#   glm(
#     in_genehancer ~ levenshtein_distance + nchar(ru) + len + fraction_gc,
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & (levenshtein_distance / len) < 0.85 ),
#     family = binomial(link = "logit")
#   )
# )

# What do these noncoding loci look like?
library(ggridges)
ggsave(
  "./figures/purity_by_ru_length_all_noncoding_jam_20241206.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & !overlap_coding),
    aes(x = (levenshtein_distance / len))
  ) +
    geom_density(adjust = 5, fill = "gray") +
    theme_minimal() +
    facet_grid(rows = vars(nchar(ru)))
)

jam_genotypes_enhancers %>%
  group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
  filter(., maf == max(maf) & !overlap_coding) %>%
  mutate(., interrupted = levenshtein_distance > 0) %>%
  group_by(., interrupted) %>%
  count(.)

######## using chromHMM states over genehancer
jam_genotypes_enhancers %>%
  group_by(., chrom, start, end, ru, fraction_gc) %>%
  filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
  group_by(., chromhmm_state) %>%
  count(.)

summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state_factor + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "het")) %>%
      mutate(., chromhmm_state_factor = as.factor(chromhmm_state)),
    family = poisson(link = "log")
  )
)

# all enhancers vs heterochromatin
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk", "het")) %>%
      mutate(., chromhmm_state = ifelse(chromhmm_state == "het", "het", "enh")),
    family = poisson(link = "log")
  )
)
# plotting this
ggsave(
  "./figures/noncoding_chromhmm_enh_vs_het_20250501.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 &
               !overlap_coding & !is.na(chromhmm_state) &
               chromhmm_state %in% c("enh_str", "enh_wk", "het")) %>%
      mutate(., is_enhancer = chromhmm_state != "het") %>%
      group_by(., is_enhancer) %>%
      slice_sample(., n = 25000),
    aes(x = levenshtein_distance / len)
  ) +
    geom_density(aes(color = is_enhancer), adjust = 4) +
    theme_minimal()
)
# het STRs are less pure
summary(
  glm(
    is_pure ~ nchar(ru) + fraction_gc + chromhmm_state + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk", "het")) %>%
      mutate(., chromhmm_state = ifelse(chromhmm_state == "het", "het", "enh"),
             is_pure = levenshtein_distance == 0),
    family = binomial()
  )
)

# longest pure stretch?
summary(
  glm(
    longest_pure_stretch ~ nchar(ru) + fraction_gc + chromhmm_state + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk", "het")) %>%
      mutate(., chromhmm_state = ifelse(chromhmm_state == "het", "het", "enh"),
             longest_pure_stretch = nchar(ru) * max_homo_ru),
    family = poisson(link = "log")
  )
)

# Regressing out purity
summary(
  glm(
    longest_pure_stretch ~ nchar(ru) + fraction_gc + chromhmm_state + levenshtein_distance + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk", "het")) %>%
      mutate(., chromhmm_state = ifelse(chromhmm_state == "het", "het", "enh"),
             longest_pure_stretch = nchar(ru) * max_homo_ru),
    family = poisson(link = "log")
  )
)

# Breaking down enhancers by their strength
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state_factor + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk")) %>%
      mutate(., chromhmm_state_factor = as.factor(chromhmm_state)),
    family = poisson(link = "log")
  )
)

# Plotting this one
ggsave(
  "./figures/noncoding_chromhmm_enh_str_vs_enh_wk_20250501.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 &
               !overlap_coding & !is.na(chromhmm_state) &
               chromhmm_state %in% c("enh_str", "enh_wk")) %>%
      group_by(., chromhmm_state) %>%
      slice_sample(., n = 2500),
    aes(x = levenshtein_distance / len)
  ) +
    geom_density(aes(color = chromhmm_state), adjust = 4) +
    theme_minimal()
)

ggsave(
  "./figures/noncoding_chromhmm_states_20250411.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & !is.na(chromhmm_state)) %>%
      mutate(., chromhmm_state_factor = as.factor(chromhmm_state)),
    aes(x = chromhmm_state_factor, y = levenshtein_distance / len)
  ) +
    geom_boxplot()
)

# downsampling to the smallest n (strong enhancers)
ggsave(
  "./figures/noncoding_chromhmm_enh_str_vs_het_20250417.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 &
               !overlap_coding & !is.na(chromhmm_state) &
               chromhmm_state %in% c("enh_str", "enh_wk", "het")) %>%
      mutate(., chromhmm_state_factor = as.factor(chromhmm_state)) %>%
      group_by(., chromhmm_state_factor) %>%
      slice_sample(., n = 2500),
    aes(x = levenshtein_distance / len)
  ) +
    geom_density(aes(color = chromhmm_state), adjust = 2)
)

# ipsc with chromhmm
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_ipsc_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_wk", "enh_str")) %>%
      mutate(., in_ipsc_peak = ipsc_atac > 0),
    family = poisson(link = "log")
  )
)

ggsave(
  "./figures/noncoding_ipsc_peaks_in_chromhmm_enh_20250501.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk")) %>%
      mutate(., in_ipsc_peak = ipsc_atac > 0) %>%
      group_by(., in_ipsc_peak) %>%
      slice_sample(., n = 6250),
    aes(x = levenshtein_distance / len)
  ) +
    geom_density(aes(color = in_ipsc_peak), adjust = 3) +
    theme_minimal()
)

# promoters are purer
ggsave(
  "./figures/noncoding_enh_vs_prom_chromhmm_20250424.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk", "prom")) %>%
      mutate(., is_promoter = chromhmm_state == "prom") %>%
      group_by(., is_promoter) %>%
      slice_sample(., n = 2748),
    aes(x = levenshtein_distance / len)
  ) +
    geom_density(aes(color = is_promoter), adjust = 3)
)

summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + is_promoter + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk", "prom")) %>%
      mutate(., is_promoter = chromhmm_state == "prom"),
    family = poisson(link = "log")
  )
)

# in fire peak
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_fire_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_wk", "enh_str")) %>%
      mutate(., in_fire_peak = fire > 0),
    family = poisson(link = "log")
  )
)

# in atac peak
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_atac_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "enh_wk")) %>%
      mutate(., in_atac_peak = k562_gm12878_atac > 0),
    family = poisson(link = "log")
  )
)

# in vista
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_vista + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_wk", "enh_str")) %>%
      mutate(., in_vista = in_vista > 0),
    family = poisson(link = "log")
  )
) # no effect

# testing the specific DNA binding proteins
sapply(
  tf_list,
  function(x)
    print(summary(
      glm(
        levenshtein_distance ~ ru + offset(log(len)),
        data = jam_genotypes_enhancers %>%
          group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
          filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
                   !overlap_coding & ru %in% c("AG", "AC") & grepl(x, overlap_tf_chip) & chromhmm_state %in% c("enh_str", "enh_wk")) %>%
          mutate(., ru = as.factor(ru)),
        family = poisson(link = "log")
      )
    ))
)
# Interesting - only ELK1 and MAX show significant effects and it is in the opposite direction

library(tidyr)
ggsave(
  "./figures/tf_chip_overlap_purity_ag_ac_20250501.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
               !overlap_coding & ru %in% c("AG", "AC") & chromhmm_state %in% c("enh_str", "enh_wk")) %>%
      mutate(., ru = as.factor(ru),
             BRCA1 = grepl("BRCA1", overlap_tf_chip),
             CTCF = grepl("CTCF", overlap_tf_chip),
             E2F1 = grepl("E2F1", overlap_tf_chip),
             E2F4 = grepl("E2F4", overlap_tf_chip),
             ELK1 = grepl("ELK1", overlap_tf_chip),
             MAX = grepl("MAX", overlap_tf_chip),
             MYC = grepl("MYC", overlap_tf_chip),
             NFE2L2 = grepl("NFE2L2", overlap_tf_chip),
             PRDM1 = grepl("PRDM1", overlap_tf_chip),
             SMARCA4 = grepl("SMARCA4", overlap_tf_chip)
      ) %>%
      pivot_longer(., cols = tf_list, names_to = "tf", values_to = "near_dbd") %>%
      filter(., near_dbd),
    aes(y = levenshtein_distance / len, x = as.factor(tf))
  ) +
    geom_boxplot(aes(fill = ru)) +
    theme_minimal()
)

# what are the properties of AG vs AC?
sapply(
  tf_list,
  function(x)
    print(
      summary(
        glm(
          levenshtein_distance ~ ru * in_tf_binding_site + offset(log(len)),
          data = jam_genotypes_enhancers %>%
            group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
            filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
                     !overlap_coding & ru %in% c("AG", "AC")) %>%
            mutate(., ru = as.factor(ru), 
                   in_tf_binding_site = grepl(x, overlap_tf_chip) & chromhmm_state %in% c("enh_str", "enh_wk")) %>%
            filter(., in_tf_binding_site | !chromhmm_state %in% c("enh_str", "enh_wk")),
          family = poisson(link = "log")
        )
      )
    )
)

# Specifically correcting for AG mutability
summary(
  glm(
    levenshtein_distance ~ ru + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
               !overlap_coding & ru %in% c("AG", "AC") & chromhmm_state %in% c("enh_str", "enh_wk")) %>%
      mutate(., ru = as.factor(ru)),
    family = poisson(link = "log")
  )
)

sapply(
  tf_list,
  function(x)
    print(
      summary(
        glm(
          levenshtein_distance ~ ru * in_tf_binding_site + offset(log(len)),
          data = jam_genotypes_enhancers %>%
            group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
            filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
                     !overlap_coding & ru %in% c("AG", "AC") & chromhmm_state %in% c("enh_str", "enh_wk")) %>%
            mutate(., ru = as.factor(ru), 
                   in_tf_binding_site = grepl(x, overlap_tf_chip)),
          family = poisson(link = "log")
        )
      )
    )
)

# Revcomp/cycle motifs as necessary
# AG are favored, AC are disfavored
sapply(
  tf_list,
  function(x)
    print(summary(
      glm(
        levenshtein_distance ~ ru + offset(log(len)),
        data = jam_genotypes_enhancers %>%
          group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
          filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
                   !overlap_coding & ru %in% c("AG", "AC") & grepl(x, overlap_tf_chip)) %>%
          mutate(., ru = as.factor(ru)),
        family = poisson(link = "log")
      )
    ))
)


######## Replicating now in genehancer

# nrow(
#   jam_genotypes_enhancers %>%
#     group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#     filter(., maf == max(maf) & !overlap_coding)
#   )

summary(
  glm(
    levenshtein_distance ~ in_genehancer + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
               !overlap_coding),
    family = poisson(link = "log")
  )
)
ggsave(
  "./figures/noncoding_purity_by_genehancer_20250502.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & !overlap_coding) %>%
      group_by(., in_genehancer) %>%
      slice_sample(., n = 125000),
    aes(x = (levenshtein_distance / len))
  ) +
    geom_density(aes(col = in_genehancer), adjust = 6) +
    theme_minimal()
)

ggplot(
  data = jam_genotypes_enhancers %>%
    group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
    filter(., maf == max(maf) & !overlap_coding),
  aes(x = (levenshtein_distance / len))
) +
  geom_density(aes(col = in_genehancer), adjust = 5) +
  theme_minimal() +
  facet_grid(rows = vars(nchar(ru)))


summary(
  glm(
    is_pure ~ in_genehancer + nchar(ru) + fraction_gc + len,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., is_pure = (levenshtein_distance == 0)),
    family = "binomial"
  )
)

summary(
  glm(
    is_pure ~ elite + nchar(ru) + fraction_gc + len,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & in_genehancer) %>%
      mutate(., is_pure = (levenshtein_distance == 0)),
    family = "binomial"
  )
)

summary(
  glm(
    levenshtein_distance ~ elite + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & in_genehancer),
    family = poisson(link = "log")
  )
)

ggsave(
  "./figures/elite_vs_nonelite_genehancer_interruption_density_20250317.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) &  !overlap_coding & (levenshtein_distance / len) < 0.25 & in_genehancer),
    aes(x = (levenshtein_distance / len))
  ) +
    geom_density(aes(col = elite), adjust = 3) +
    theme_minimal() +
    facet_grid(rows = vars(nchar(ru)))
)

# Stratifying by elite score - very weak effect
ggsave(
  "./figures/purity_by_loeuf_elite_genehancers_jam_20241120.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & min_elite_loeuf != -1),
    aes(x = min_elite_loeuf, y = (levenshtein_distance / len))
  ) +
    geom_point(alpha = 0.1) +
    theme_minimal() +
    stat_smooth() +
    facet_grid(rows = vars(nchar(ru)), cols = vars(fraction_gc > 0))
)

## what about just stratifying by whether the STR is in an elite enhancer vs not in an enhancer
summary(
  glm(
    levenshtein_distance ~ elite_enhancer + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., elite_enhancer = ifelse(in_genehancer, ifelse(elite == 1, TRUE, NA), FALSE)),
    family = poisson(link = "log")
  )
)

ggsave(
  "./figures/noncoding_elite_enhancers_vs_not_enhancers_by_ru_gc_112024.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., elite_enhancer = ifelse(in_genehancer, ifelse(elite == 1, TRUE, NA), FALSE),
             contains_gc = fraction_gc > 0) %>%
      filter(., !is.na(elite_enhancer)),
    aes(x = (levenshtein_distance / len))
  ) +
    geom_density(aes(col = elite_enhancer), adjust = 5) +
    theme_minimal() +
    facet_grid(rows = vars(nchar(ru)), cols = vars(contains_gc))
)

# Stratifying by loeuf score
summary(
  glm(
    levenshtein_distance ~ min_elite_loeuf + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & min_elite_loeuf != -1),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    is_pure ~ min_elite_loeuf + nchar(ru) + fraction_gc + len,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & min_elite_loeuf != -1) %>%
      mutate(., is_pure = (levenshtein_distance == 0)),
    family = "binomial"
  )
)

# Weird - the STRs in elite enhancers linked to autosomal dominant genes are not significantly more likely to be pure
summary(
  glm(
    is_pure ~ genehancer_gene_is_ad + nchar(ru) + fraction_gc + len,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & genehancer_gene_is_ad != -1 & elite == 1) %>%
      mutate(., is_pure = (levenshtein_distance == 0)),
    family = "binomial"
  )
)

summary(
  glm(
    levenshtein_distance ~ genehancer_gene_is_ad + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & genehancer_gene_is_ad != -1 & elite == 1) %>%
      mutate(., genehancer_gene_is_ad = as.logical(genehancer_gene_is_ad)),
    family = poisson(link = "log")
  )
)

ggsave(
  "./figures/noncoding_elite_enhancers_is_gene_ad_20241120.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & genehancer_gene_is_ad != -1 & elite == 1) %>%
      mutate(., genehancer_gene_is_ad = as.logical(genehancer_gene_is_ad)),
    aes(x = (levenshtein_distance / len))
  ) +
    geom_density(aes(color = genehancer_gene_is_ad), adjust = 3) +
    theme_minimal()
)

# Any difference between enhancers and promoter/enhancers?
summary(
  glm(
    levenshtein_distance ~ element_type + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & elite == 1 &
               gh_element %in% c("Enhancer", "Promoter")) %>%
      mutate(., element_type = as.factor(gh_element)),
    family = poisson(link = "log")
  )
)

# FIRE peaks
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_fire_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, fire) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & elite == "1") %>%
      mutate(., in_fire_peak = fire > 0),
    family = poisson(link = "log")
  )
)

jam_genotypes_enhancers %>%
  group_by(., chrom, start, end, ru, fraction_gc, fire) %>%
  filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & elite == "1") %>%
  mutate(., in_fire_peak = fire > 0) %>%
  group_by(., in_fire_peak) %>%
  count()

# does this replicate in atac?
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_atac_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & elite == "1") %>%
      mutate(., in_atac_peak = k562_gm12878_atac > 0),
    family = poisson(link = "log")
  )
)

# IPSC chip
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_ipsc_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & elite == "1") %>%
      mutate(., in_ipsc_peak = ipsc_atac > 0),
    family = poisson(link = "log")
  )
)

# Any association with problematic regions?
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + problematic + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & elite == "1") %>%
      mutate(., problematic = problematic_region > 0),
    family = poisson(link = "log")
  )
)

# What about VISTA enhancers?
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_vista + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., in_vista = in_vista > 0),
    family = poisson(link = "log")
  )
) # no effect

#### FIRE peaks - accessibility
ggsave(
  "./figures/fire_peaks_by_purity_20250225.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, fire) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., in_fire_peak = fire > 0),
    aes(x = levenshtein_distance / len)
  ) +
    geom_density(aes(color = in_fire_peak), adjust = 5)
)

summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_fire_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, fire) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., in_fire_peak = fire > 0),
    family = poisson(link = "log")
  )
)

# does this replicate in ATAC?
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_atac_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., in_atac_peak = k562_gm12878_atac > 0),
    family = poisson(link = "log")
  )
)

# what about in relevant cell type?
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + in_atac_peak + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., in_atac_peak = ipsc_atac > 0),
    family = poisson(link = "log")
  )
)

###########################################
# any association with specific non-B structures?
major_alleles_shapes <- read.delim("./major_alleles_shapes.txt", header = F)
colnames(major_alleles_shapes) <-
  c("chrom", "start", "DR_obs", "DR_delta_perfect", "MR_obs", "MR_delta_perfect", "IR_obs", "IR_delta_perfect",
    "APR_obs", "APR_delta_perfect", "GQ_obs", "GQ_delta_perfect", "Z_obs", "Z_delta_perfect")

major_alleles_shapes <-
  merge(
    major_alleles_shapes,
    jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf)) %>%
      slice(., 1L),
    by = c("chrom", "start"), all.x = T, all.y = F)

library(stringr)
library(purrr)
library(broom)
major_alleles_shapes <-
  major_alleles_shapes %>%
  mutate(., across(ends_with("_obs"), ~ ifelse(. > 0, 1, 0), .names = "{.col}_contains")) %>%
  rename_at(., vars(ends_with("_contains")), ~ str_remove(., "\\_obs"))

ggsave(
  "./figures/nonb_shape_frequency_20260106.pdf",
  ggplot(
    data = major_alleles_shapes %>%
      filter(., (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      select(., contains("_contains")) %>%
      pivot_longer(., cols = contains("_contains"), names_to = "nonb_shape", values_to = "contains_shape"),
    aes(x = nonb_shape, fill = as.factor(contains_shape))
  ) +
    geom_bar(position = "fill") +
    theme_minimal()
)

major_alleles_shapes %>%
  mutate(interruption_density = levenshtein_distance / len) %>%
  filter(., interruption_density < 0.25 & !overlap_coding) %>%
  select(., interruption_density, contains("_contains")) %>%
  pivot_longer(., cols = contains("_contains"), names_to = "nonb_shape", values_to = "contains_shape") %>%
  group_by(., nonb_shape) %>%
  nest() %>%
  mutate(
    test = map(data, ~ wilcox.test(interruption_density ~ contains_shape, data = .x)),
    tidy = map(test, tidy)
  ) %>%
  unnest(tidy)

major_alleles_shapes %>%
  mutate(interruption_density = levenshtein_distance / len) %>%
  filter(., interruption_density < 0.25 & !overlap_coding) %>%
  select(., levenshtein_distance, contains("_contains"), ru, fraction_gc, len) %>%
  pivot_longer(., cols = contains("_contains"), names_to = "nonb_shape", values_to = "contains_shape") %>%
  group_by(., nonb_shape) %>%
  nest() %>%
  mutate(
    test = map(data, ~ glm(levenshtein_distance ~ nchar(ru) + fraction_gc + contains_shape + offset(log(len)), data = .x)),
    tidy = map(test, tidy)
  ) %>%
  unnest(tidy) %>%
  print(n = 25)

major_alleles_shapes %>%
  mutate(interruption_density = levenshtein_distance / len) %>%
  filter(., interruption_density < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "het")) %>%
  mutate(., chromhmm_state_factor = as.factor(chromhmm_state)) %>%
  select(., levenshtein_distance, contains("_contains"), ru, fraction_gc, len, chromhmm_state_factor) %>%
  pivot_longer(., cols = contains("_contains"), names_to = "nonb_shape", values_to = "contains_shape") %>%
  group_by(., nonb_shape) %>%
  nest() %>%
  mutate(
    test = map(data, ~ glm(levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state_factor * contains_shape + offset(log(len)), 
                           data = .x, family = poisson(link = "log"))),
    tidy = map(test, tidy)
  ) %>%
  unnest(tidy) %>%
  print(n = 50)

major_alleles_shapes %>%
  mutate(interruption_density = levenshtein_distance / len) %>%
  filter(., interruption_density < 0.25 & !overlap_coding & chromhmm_state %in% c("prom", "het")) %>%
  mutate(., chromhmm_state_factor = as.factor(chromhmm_state)) %>%
  select(., levenshtein_distance, contains("_contains"), ru, fraction_gc, len, chromhmm_state_factor) %>%
  pivot_longer(., cols = contains("_contains"), names_to = "nonb_shape", values_to = "contains_shape") %>%
  group_by(., nonb_shape) %>%
  nest() %>%
  mutate(
    test = map(data, ~ glm(levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state_factor * contains_shape + offset(log(len)), 
                           data = .x, family = poisson(link = "log"))),
    tidy = map(test, tidy)
  ) %>%
  unnest(tidy) %>%
  print(n = 50)

# AIC values for the models without shape are 1706100 for prom v het and 1705403 for enh_str v het
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state_factor + offset(log(len)),
    data = 
      major_alleles_shapes %>%
      mutate(interruption_density = levenshtein_distance / len) %>%
      filter(., interruption_density < 0.25 & !overlap_coding & chromhmm_state %in% c("prom", "het")) %>%
      mutate(., chromhmm_state_factor = as.factor(chromhmm_state)),
    family = poisson(link = "log")
  )
)
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state_factor + offset(log(len)),
    data = 
      major_alleles_shapes %>%
      mutate(interruption_density = levenshtein_distance / len) %>%
      filter(., interruption_density < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "het")) %>%
      mutate(., chromhmm_state_factor = as.factor(chromhmm_state)),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state_factor * IR_contains + offset(log(len)),
    data = 
      major_alleles_shapes %>%
      mutate(interruption_density = levenshtein_distance / len) %>%
      filter(., interruption_density < 0.25 & !overlap_coding & chromhmm_state %in% c("prom", "het")) %>%
      mutate(., chromhmm_state_factor = as.factor(chromhmm_state)) %>%
      select(., levenshtein_distance, contains("_contains"), ru, fraction_gc, len, chromhmm_state_factor),
    family = poisson(link = "log")
  )
)
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + chromhmm_state_factor * IR_contains + offset(log(len)),
    data = 
      major_alleles_shapes %>%
      mutate(interruption_density = levenshtein_distance / len) %>%
      filter(., interruption_density < 0.25 & !overlap_coding & chromhmm_state %in% c("enh_str", "het")) %>%
      mutate(., chromhmm_state_factor = as.factor(chromhmm_state)) %>%
      select(., levenshtein_distance, contains("_contains"), ru, fraction_gc, len, chromhmm_state_factor),
    family = poisson(link = "log")
  )
)

ggsave(
  "./figures/nonb_shape_interruptions_20260106.pdf",
  ggplot(
    data = 
      major_alleles_shapes %>%
      select(., levenshtein_distance, len, contains("_contains")) %>%
      pivot_longer(., cols = contains("_contains"), names_to = "nonb_shape", values_to = "contains_shape"),
    aes(x = nonb_shape, y = levenshtein_distance / len)
  ) +
    geom_boxplot(aes(fill = as.factor(contains_shape)))
)