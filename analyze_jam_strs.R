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


european_eqtls <- read.csv("./jam_supp14_estrs_eur.csv")
# some strs are eqtls for multiple sites. i am interested if any of those associations
# are significant or are the best str. currently i'm averaging the slopes of the associations...
# maybe this is not the most precise bet but it'll work for now i guess because 
# slope of association hasn't really been much of a problem
jam_genotypes_enhancers <-
  merge(
    jam_genotypes_enhancers,
    european_eqtls %>% select(., chrom, start, slope, sig.genelevel, best.str) %>%
      group_by(., chrom, start) %>%
      mutate(., best.str = any(best.str), sig.genelevel = any(sig.genelevel), slope = mean(slope)) %>%
      distinct(., chrom, start, slope, sig.genelevel, best.str),
    by.x = c("chrom", "start"), by.y = c("chrom", "start"),
    all.x = T
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

##### QTL analysis
# genome wide significant eqtls
summary(
  glm(
    levenshtein_distance ~ is_sig_eqtl + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, sig.genelevel) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & !contains_cpg) %>%
      mutate(., is_sig_eqtl = !is.na(sig.genelevel) & sig.genelevel),
    family = poisson(link = "log")
  )
)

# best eqtl from gene
summary(
  glm(
    levenshtein_distance ~ best.str + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, sig.genelevel, best.str) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & !contains_cpg) %>%
      mutate(., best.str = !is.na(best.str) & best.str),
    family = poisson(link = "log")
  )
)

# conditional on eqtl being significant, what about the sign of the eqtl?
summary(
  glm(
    levenshtein_distance ~ slope_positive + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, sig.genelevel, slope) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & !contains_cpg & sig.genelevel) %>%
      mutate(., slope_positive = slope > 0),
    family = poisson(link = "log")
  )
)
# no, these fit poorly
summary(
  glm(
    levenshtein_distance ~ slope_positive + nchar(ru) + fraction_gc + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, best.str, slope) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & !contains_cpg & best.str) %>%
      mutate(., slope_positive = slope > 0),
    family = poisson(link = "log")
  )
)

# density plot of purity at significant eqtls
ggsave(
  "./figures/significant_eqtl_purity_density_20240404.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, sig.genelevel) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.85 & !overlap_coding & !contains_cpg) %>%
      mutate(., is_sig_eqtl = !is.na(sig.genelevel) & sig.genelevel),
    aes(x = (levenshtein_distance / len))
  ) +
    geom_density(aes(color = is_sig_eqtl), adjust = 5) +
    facet_grid(rows = vars(nchar(ru))) +
    theme_minimal()
)

# Revcomp/cycle motifs as necessary
# AG are favored, AC are disfavored
# Seems like regularization has already been done
# summary(
#   glm(
#   levenshtein_distance ~ ru + offset(log(len)),
#   data = jam_genotypes_enhancers %>%
#     group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
#     filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
#              !overlap_coding & ru %in% c("AG", "AC") & grepl("MAX", overlap_tf_chip)) %>%
#     mutate(., ru = as.factor(ru)),
#   family = poisson(link = "log")
# )
# )

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

sapply(
  tf_list,
  function(x)
    print(summary(
      glm(
        max_len_homology ~ levenshtein_distance + ru + offset(log(len)),
        data = jam_genotypes_enhancers %>%
          group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
          filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
                   !overlap_coding & ru %in% c("AG", "AC") & grepl(x, overlap_tf_chip)) %>%
          mutate(., max_len_homology = nchar(ru) * max_homo_ru,
                 ru = as.factor(ru)),
        family = poisson(link = "log")
      )
    ))
)

sapply(
  tf_list,
  function(x)
    summary(
      glm(
        max_len_homology ~ levenshtein_distance + ru + offset(log(len)),
        data = jam_genotypes_enhancers %>%
          group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
          filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
                   !overlap_coding & ru %in% c("AG", "AC") & grepl(x, overlap_tf_chip)) %>%
          mutate(., max_len_homology = nchar(ru) * max_homo_ru,
                 ru = as.factor(ru)),
        family = poisson(link = "log")
      )
    )$coefficients[2, c(1, 4)]
)


######## CpG islands
ggsave(
  "./figures/cpgi_purity_by_nuc_content_20241209.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., cpgi = cpgi > 0),
    aes(x = levenshtein_distance / len)
  ) +
    geom_density(aes(color = cpgi), adjust = 5) +
    facet_grid(rows = vars(contains_cpg))
)

summary(
  glm(
    levenshtein_distance ~ nchar(ru) + contains_cpg * cpgi + offset(log(len)),
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., cpgi = cpgi > 0),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    is_pure ~ nchar(ru) + contains_cpg * cpgi + len,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, overlap_tf_chip) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding) %>%
      mutate(., cpgi = cpgi > 0, is_pure = levenshtein_distance == 0),
    family = binomial()
  )
)

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

#### 

# jam_genotypes_enhancers %>%
#   group_by(., chrom, start, end, ru, fraction_gc, overlap_max_peak) %>%
#   filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
#            !overlap_coding & overlap_max_peak > 0 & ru %in% c("AG", "AC")) %>%
#   mutate(., ru = as.factor(ru)) %>%
#   group_by(., ru) %>%
#   count()

ggsave(
  "./figures/noncoding_max_overlap_strs_AG_vs_AC_20241011.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, overlap_max_peak) %>%
      filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & 
               !overlap_coding & overlap_max_peak > 0 & ru %in% c("AG", "AC")) %>%
      mutate(., ru = as.factor(ru)),
    aes(x = ru, y = levenshtein_distance / len)
  ) +
    geom_boxplot() 
)


# summary(
#   glm(
#     levenshtein_distance ~ in_genehancer + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & !overlap_coding & contains_cpg),
#     family = poisson(link = "log")
#   )
# )

# summary(
#   glm(
#     levenshtein_distance ~ min_elite_loeuf + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & min_elite_loeuf != -1),
#     family = poisson(link = "log")
#   )
# )

# summary(
#   glm(
#     levenshtein_distance ~ elite + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & in_genehancer),
#     family = poisson(link = "log")
#   )
# )

#### What happens when i exclude pure ones?
# summary(
#   glm(
#     levenshtein_distance ~ in_genehancer + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & levenshtein_distance > 0 & !overlap_coding & contains_cpg),
#     family = poisson(link = "log")
#   )
# )
# 
# summary(
#   glm(
#     levenshtein_distance ~ elite + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & levenshtein_distance > 0 & in_genehancer),
#     family = poisson(link = "log")
#   )
# )
# 


# summary(
#   glm(
#     levenshtein_distance ~ min_elite_loeuf + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25 & levenshtein_distance > 0 & min_elite_loeuf != -1),
#     family = poisson(link = "log")
#   )
# )


# ggplot(
#   data = jam_genotypes_enhancers %>%
#     group_by(., chrom, start, end, ru, fraction_gc, in_genehancer, elite) %>%
#     filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25) %>%
#     mutate(., enhancer_status = ifelse(in_genehancer, ifelse(elite == 1, "elite", "non-elite"), "non-enhancer"),
#            nucleotide_content = ifelse(fraction_gc == 0, "AT-only", "GC-containing")),
#   aes(x = (levenshtein_distance))
# ) +
#   geom_density(aes(color = enhancer_status), adjust = 5) +
#   # xlim(c(0, 0.2)) +
#   facet_grid(rows = vars(nucleotide_content))

# ggplot(
#   data = jam_genotypes_enhancers %>%
#     group_by(., chrom, start, end, ru, fraction_gc, in_genehancer, elite) %>%
#     filter(., maf == max(maf) & (levenshtein_distance / len) < 0.25) %>%
#     mutate(., enhancer_status = ifelse(in_genehancer, ifelse(elite == 1, "elite", "non-elite"), "non-enhancer"),
#            nucleotide_content = ifelse(fraction_gc == 0, "AT-only", "GC-containing")),
#   aes(x = (levenshtein_distance))
# ) +
#   geom_density(aes(color = enhancer_status), adjust = 5) +
#   # xlim(c(0, 0.2)) +
#   facet_grid(rows = vars(nucleotide_content))



# library(betareg)
# summary(
#   betareg(
#     purity ~ in_vista + nchar(ru) + len + fraction_gc,
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15) %>%
#       mutate(., purity = purity - 0.0001)
#   )
# )

# summary(
#   glm(
#     purity ~ in_genehancer + nchar(ru) + len + fraction_gc,
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15),
#     family = quasibinomial(link = "logit")
#   )
# )

# what if i tried a regression of the number of bases unaccounted for by a repeat offset by total length in bp of repeat?

# summary(
#   glm(
#     max_len_homology ~ elite + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.25 & in_genehancer) %>%
#       mutate(., max_len_homology = nchar(ru) * max_homo_ru),
#     family = poisson(link = "log")
#   )
# )

# summary(
#   glm(
#     max_len_homology ~ in_genehancer + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15) %>%
#       mutate(., max_len_homology = nchar(ru) * max_homo_ru),
#     family = poisson(link = "log")
#   )
# )

# summary(
#   glm(
#     max_len_homology ~ min_elite_loeuf + nchar(ru) + fraction_gc + offset(log(len)),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15 & min_elite_loeuf != -1) %>%
#       mutate(., max_len_homology = nchar(ru) * max_homo_ru),
#     family = poisson(link = "log")
#   )
# )

summary(
  glm(
    max_len_homology ~ in_genehancer + nchar(ru) + fraction_gc + offset(log(len)) + levenshtein_distance,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & purity > 0.25 & !overlap_coding) %>%
      mutate(., max_len_homology = nchar(ru) * max_homo_ru),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    max_len_homology ~ min_elite_loeuf + nchar(ru) + fraction_gc + offset(log(len)) + levenshtein_distance,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      filter(., maf == max(maf) & purity > 0.25 & min_elite_loeuf != -1 & !overlap_coding) %>%
      mutate(., max_len_homology = nchar(ru) * max_homo_ru),
    family = poisson(link = "log")
  )
)

# summary(
#   glm(
#     max_len_homology ~ elite + nchar(ru) + fraction_gc + offset(log(len)) + levenshtein_distance,
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15 & in_genehancer) %>%
#       mutate(., max_len_homology = nchar(ru) * max_homo_ru),
#     family = poisson(link = "log")
#   )
# )


# ggplot(
#   data = jam_genotypes_enhancers %>%
#     group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#     filter(., maf == max(maf) & purity > 0.15),
#   aes(x = len)
# ) +
#   geom_bar(aes(fill = in_genehancer), position = "fill")

# ggplot(
#   data = jam_genotypes_enhancers %>%
#     group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#     filter(., maf == max(maf) & purity > 0.15 & len > 10),
#   aes(x = len)
# ) +
#   geom_histogram()

# summary(
#   glm(
#     len ~ in_genehancer + fraction_gc + nchar(ru),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15),
#     family = poisson(link = "log")
#   )
# )

# summary(
#   glm(
#     len ~ elite + fraction_gc + nchar(ru),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15 & in_genehancer),
#     family = poisson(link = "log")
#   )
# )

# summary(
#   glm(
#     len ~ min_elite_loeuf + fraction_gc + nchar(ru),
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15 & min_elite_loeuf != -1),
#     family = poisson(link = "log")
#   )
# )

# summary(
#   glm(
#     max_len_homology ~ min_elite_loeuf + nchar(ru) + fraction_gc + offset(log(len)) + levenshtein_distance,
#     data = jam_genotypes_enhancers %>%
#       group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
#       filter(., maf == max(maf) & purity > 0.15 & min_elite_loeuf != -1) %>%
#       mutate(., max_len_homology = nchar(ru) * max_homo_ru),
#     family = poisson(link = "log")
#   )
# )



# Some analyses on how many nonsingleton alleles are observed ... these need work i think
summary(
  glm(
    n_alleles ~ in_genehancer + nchar(ru) + fraction_gc + len + levenshtein_distance,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      mutate(., n_alleles = n()) %>%
      filter(., maf == max(maf) & purity > 0.15) %>%
      mutate(., max_len_homology = nchar(ru) * max_homo_ru),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    n_alleles ~ in_genehancer + nchar(ru) + fraction_gc + len + levenshtein_distance,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      mutate(., n_alleles = n()) %>%
      filter(., maf == max(maf) & purity > 0.15) %>%
      mutate(., max_len_homology = nchar(ru) * max_homo_ru),
    family = poisson(link = "log")
  )
)

ggplot(
  data = jam_genotypes_enhancers,
  aes(x = purity)
) +
  geom_histogram()

ggplot(
  data = jam_genotypes_enhancers %>%
    mutate(., max_len_homology = nchar(ru) * max_homo_ru),
  aes(x = purity, y = log(max_len_homology))
) + geom_point()

summary(
  glm(
    n_alleles ~ in_genehancer + nchar(ru) + fraction_gc + len + levenshtein_distance,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      mutate(., n_alleles = n()) %>%
      filter(., maf == max(maf) & purity > 0.15) %>%
      mutate(., max_len_homology = nchar(ru) * max_homo_ru),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    n_alleles ~ in_genehancer + nchar(ru) + fraction_gc + len + levenshtein_distance,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, in_genehancer) %>%
      mutate(., n_alleles = n()) %>%
      filter(., maf == max(maf) & purity > 0.15) %>%
      mutate(., max_len_homology = nchar(ru) * max_homo_ru),
    family = poisson(link = "log")
  )
)

# summary(
#   glm(
#     levenshtein_distance ~ fraction_GC + offset(log(len)) + 
#   )
# )

table(
  subset(
    nchar(jam_genotypes_enhancers$ru),
    !is.na(jam_genotypes_enhancers$coding_transcript_loeuf)
  )
)

# do we want to use the same purity standards for the coding strs?
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + offset(log(len)) + coding_transcript_loeuf,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
      filter(., maf == max(maf) & purity > 0.75 & !is.na(coding_transcript_loeuf) & nchar(ru) %in% c(3, 6)),
    family = poisson(link = "log")
  )
)

###########################################
# any association with specific shapes?
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

nonb_shapes <- c("APR", "DR", "GQ", "IR", "MR", "Z")

##############################################################################
# coding analysis: how do synonymous variants covary with loeuf

# How interrupted are coding STRs in general
library(ggridges)
ggsave(
  "./figures/purity_by_ru_length_all_coding_jam_20240605.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc) %>%
      filter(., maf == max(maf) & overlap_coding),
    aes(x = (levenshtein_distance / len))
  ) +
    geom_density(adjust = 5, fill = "gray") +
    theme_minimal() +
    facet_grid(rows = vars(nchar(ru)))
)

#
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + offset(log(len)) + coding_gene_loeuf,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               (interruption_coding_annotation == "synonymous_variant" |
                  levenshtein_distance == 0)),
    family = poisson(link = "log")
  )
)

# AD genes?
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + offset(log(len)) + coding_gene_is_ad,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_is_ad != -1 & 
               nchar(ru) %in% c(3, 6) &
               (interruption_coding_annotation == "synonymous_variant" |
                  levenshtein_distance == 0)),
    family = poisson(link = "log")
  )
)
# plot this?
ggsave(
  "./figures/purity_coding_strs_ad_genes_20240731.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_is_ad != -1 & 
               nchar(ru) %in% c(3, 6) &
               (interruption_coding_annotation == "synonymous_variant" |
                  levenshtein_distance == 0)) %>%
      mutate(., distance_quantile = cut(levenshtein_distance / len, breaks = seq(0, 1, 0.05))) %>%
      mutate(., distance_quantile = factor(distance_quantile, levels = rev(levels(distance_quantile)))),
    aes(y = levenshtein_distance / len, x = as.factor(coding_gene_is_ad))
  ) +
    geom_bar(aes(fill = distance_quantile), position = "fill", stat = "identity") +
    facet_grid(rows = vars(as.factor(nchar(ru))))
)


# do we see the same thing with missense variants?
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + offset(log(len)) + coding_transcript_loeuf,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               !is.na(coding_transcript_loeuf) & 
               nchar(ru) %in% c(3, 6) &
               (interruption_coding_annotation == "missense_variant" |
                  levenshtein_distance == 0)),
    family = poisson(link = "log")
  )
)

# pfam domains
summary(
  glm(
    levenshtein_distance ~ nchar(ru) + fraction_gc + offset(log(len)) + in_domain * coding_transcript_loeuf,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               !is.na(coding_transcript_loeuf) & 
               nchar(ru) %in% c(3, 6) &
               (interruption_coding_annotation == "synonymous_variant" |
                  levenshtein_distance == 0)) %>%
      mutate(., in_domain = coding_pfam_domain != "."),
    family = poisson(link = "log")
  )
)


# change regression line to linear?
ggsave("./figures/coding_strs_leouf_distance_20240119.pdf",
       ggplot(
         data = jam_genotypes_enhancers %>%
           group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
           filter(., maf == max(maf) & purity > 0.75 & !is.na(coding_transcript_loeuf) & nchar(ru) %in% c(3, 6)),
         aes(x = coding_transcript_loeuf, y = levenshtein_distance / len)
       ) + 
         geom_jitter(height = 0.01, alpha = 0.1) + 
         theme_minimal() +
         stat_smooth())
ggsave("./figures/coding_strs_leouf_distance_synonymous_20240307.pdf",
       ggplot(
         data = jam_genotypes_enhancers %>%
           group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
           filter(., maf == max(maf) & 
                    purity > 0.75 & 
                    !is.na(coding_transcript_loeuf) & 
                    nchar(ru) %in% c(3, 6) &
                    (interruption_coding_annotation == "synonymous_variant" |
                       levenshtein_distance == 0)),
         aes(x = coding_transcript_loeuf, y = levenshtein_distance / len)
       ) + 
         geom_point(alpha = 0.1) + 
         theme_minimal() +
         stat_smooth()
)
ggsave("./figures/coding_strs_leouf_distance_synonymous_in_domain_20240418.pdf",
       ggplot(
         data = jam_genotypes_enhancers %>%
           group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
           filter(., maf == max(maf) & 
                    purity > 0.75 & 
                    !is.na(coding_transcript_loeuf) & 
                    nchar(ru) %in% c(3, 6) &
                    (interruption_coding_annotation == "synonymous_variant" |
                       levenshtein_distance == 0))  %>%
           mutate(., in_domain = coding_pfam_domain != "."),
         aes(x = coding_transcript_loeuf, y = levenshtein_distance / len)
       ) + 
         geom_point(aes(col = in_domain), alpha = 0.1) + 
         theme_minimal() +
         stat_smooth(formula = y ~ x, method = 'glm', aes(col = in_domain))
)

# 
jam_genotypes_enhancers %>%
    group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
    filter(., maf == max(maf) & 
             purity > 0.75 & 
             !is.na(coding_transcript_loeuf) & 
             nchar(ru) %in% c(3, 6) &
             (interruption_coding_annotation == "synonymous_variant" |
                levenshtein_distance == 0))  %>%
    mutate(., in_domain = coding_pfam_domain != ".") %>%
  group_by(., in_domain) %>%
  count()

median_coding_loeuf <- 
  median(jam_genotypes_enhancers %>%
           group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
           filter(., maf == max(maf) & 
                    purity > 0.75 & 
                    !is.na(coding_transcript_loeuf) & 
                    nchar(ru) %in% c(3, 6) &
                    (interruption_coding_annotation == "synonymous_variant" |
                       levenshtein_distance == 0)) %>%
           pull(., coding_transcript_loeuf)
  )
ggsave(
  "./figures/synonymous_interruption_purity_constraint_density_20240307.pdf",
  ggplot(
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               !is.na(coding_transcript_loeuf) & 
               nchar(ru) %in% c(3, 6) &
               (interruption_coding_annotation == "synonymous_variant" |
                  levenshtein_distance == 0)) %>%
      mutate(., is_constrained = coding_transcript_loeuf < median_coding_loeuf),
    aes(x = levenshtein_distance / len)
  ) +
    geom_density(aes(color = is_constrained), adjust = 5) +
    theme_minimal()
)




summary(
  glm(
    max_len_homology ~ nchar(ru) + fraction_gc + offset(log(len)) + coding_transcript_loeuf + levenshtein_distance,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
      filter(., maf == max(maf) & purity > 0.75 & !is.na(coding_transcript_loeuf) & nchar(ru) %in% c(3, 6)),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    max_len_homology ~ nchar(ru) + fraction_gc + offset(log(len)) + coding_transcript_loeuf + levenshtein_distance,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               !is.na(coding_transcript_loeuf) & 
               nchar(ru) %in% c(3, 6) &
               (interruption_coding_annotation == "synonymous_variant" |
                  levenshtein_distance == 0)),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    max_len_homology ~ nchar(ru) + fraction_gc + offset(log(len)) + levenshtein_distance + overlap_coding,
    data = jam_genotypes_enhancers %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
      filter(., maf == max(maf) & purity > 0.75 & nchar(ru) %in% c(3, 6)),
    family = poisson(link = "log")
  )
)

jam_genotypes_enhancers %>%
  group_by(., chrom, start, end, ru, fraction_gc, coding_transcript_loeuf) %>%
  count()
