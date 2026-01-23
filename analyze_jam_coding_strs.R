# playing around with coding STRs only

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Documents/Quinlan/Interruptions/")

coding_strs <- read.delim("./coding_interruptions_full_annotations.txt", 
                          header = F)

colnames(coding_strs) <-
  c("chrom", "start", "end", "ru", "len", "max_homo_ru", "levenshtein_distance", "maf", 
    "in_vista", "gh_element", "elite", "min_elite_loeuf", "genehancer_gene_is_ad", "overlap_coding_genes", 
    "overlap_coding_transcripts", "coding_gene_loeuf", "coding_gene_is_ad", "n_codons", "n_missense", "n_silent",
    "n_possible_missense", "n_possible_silent", "interruption_coding_annotation", "coding_pfam_domain")

coding_strs$"fraction_gc" <- 
  sapply(coding_strs$ru,
         function(x)
           nchar(gsub("[AT]+", "", x)) / 
           nchar(x))
coding_strs$"in_genehancer" <- coding_strs$gh_element != "."
coding_strs$"max_len_homology" <-
  coding_strs$max_homo_ru * nchar(coding_strs$ru)
coding_strs$"purity" <-
  1 - (coding_strs$levenshtein_distance / coding_strs$len)
coding_strs$elite <-
  sapply(
    coding_strs$elite,
    function(x)
      if(x == "0,1") "1" else x
  )
coding_strs$"overlap_coding" <- coding_strs$overlap_coding_transcripts != "."
coding_strs$overlap_coding_transcripts <- NULL
coding_strs$coding_gene_loeuf <- 
  sapply(
    coding_strs$coding_gene_loeuf,
    function(x)
      if(x == -1) NA else x
  )
coding_strs$"contains_cpg" <-
  sapply(
    coding_strs$ru,
    function(x)
      grepl("CG", paste0(x, x))
  )
coding_strs <-
  coding_strs %>%
  mutate(., n_codons = ifelse(n_codons == ".", NA, as.integer(n_codons))) %>%
  mutate(., n_missense = ifelse(n_missense == ".", NA, as.integer(n_missense))) %>%
  mutate(., n_silent = ifelse(n_silent == ".", NA, as.integer(n_silent))) %>%
  mutate(., n_possible_missense = ifelse(n_possible_missense == ".", NA, as.integer(n_possible_missense))) %>%
  mutate(., n_possible_silent = ifelse(n_possible_silent == ".", NA, as.integer(n_possible_silent)))

### summarizing coding STR loci
coding_strs %>%
  group_by(., chrom, start) %>%
  count()
# 2838

View(
  coding_strs %>%
  group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
  filter(., maf == max(maf) & !is.na(coding_gene_loeuf)))

# How interrupted are they in general
ggsave(
  "./figures/coding_strs_interruption_density_silent_vs_missense_20250314.pdf",
  ggplot(
    data = coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      pivot_longer(., cols = c(n_silent, n_missense), names_to = "interruption_consequence", values_to = "n") %>%
      mutate(., reweighted_n = if_else(interruption_consequence == "n_missense", n * n_possible_silent / n_possible_missense, as.double(n))),
    aes(x = reweighted_n / n_codons)
  ) +
    geom_density(aes(color = interruption_consequence), adjust = 3) +
    theme_minimal()
)

# how many are interrupted?
coding_strs %>%
  group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
  filter(., maf == max(maf) & 
           purity > 0.75 & 
           coding_gene_loeuf != -1 & 
           nchar(ru) %in% c(3, 6) &
           n_codons != ".") %>%
  mutate(., has_interruption_silent = n_silent >= 1,
         has_interruption_missense = n_missense >= 1) %>%
  group_by(., has_interruption_silent, has_interruption_missense) %>%
  count()

# can i do an offset(log(n_possible_missense * n_codons)) with a dependent variable of log(n * n_silent)

summary(
  glm(
    n_silent ~ nchar(ru) + fraction_gc + offset(log(n_codons)) + coding_gene_loeuf,
    data = coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != "."),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    n_missense ~ nchar(ru) + fraction_gc + offset(log(n_codons)) + coding_gene_loeuf,
    data = coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != "."),
    family = poisson(link = "log")
  )
)

# put them together?
summary(
  glm(
    n ~ nchar(ru) + fraction_gc + coding_gene_loeuf * interruption_consequence + offset(log(n_codons)),
    data = 
      coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      pivot_longer(., cols = c(n_silent, n_missense), names_to = "interruption_consequence", values_to = "n") %>%
      mutate(., consequence_offset = ifelse(interruption_consequence == "n_missense", n_possible_missense, n_possible_silent)),
    family = poisson(link = "log")
  )
) # significant interaction

# let's try it with a fancy offset
summary(
  glm(
    n ~ nchar(ru) + fraction_gc + coding_gene_loeuf * interruption_consequence + 
      offset(log(n_codons) + log(consequence_offset) - log(n_possible_silent)),
    data = 
      coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      pivot_longer(., cols = c(n_silent, n_missense), names_to = "interruption_consequence", values_to = "n") %>%
      mutate(., consequence_offset = ifelse(interruption_consequence == "n_missense", n_possible_missense, n_possible_silent)),
    family = poisson(link = "log")
  )
)

# without LOEUF:
summary(
  glm(
    n ~ nchar(ru) + fraction_gc + interruption_consequence + 
      offset(log(n_codons) + log(consequence_offset) - log(n_possible_silent)),
    data = 
      coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      pivot_longer(., cols = c(n_silent, n_missense), names_to = "interruption_consequence", values_to = "n") %>%
      mutate(., consequence_offset = ifelse(interruption_consequence == "n_missense", n_possible_missense, n_possible_silent)),
    family = poisson(link = "log")
  )
)

# imagining a pdf of fraction codons interrupted, broken out by consequence of interruption
ggsave(
  "./figures/coding_strs_loeuf_synonymous_vs_missense_20241017.pdf",
  ggplot(
    data = 
      coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      pivot_longer(., cols = c(n_silent, n_missense), names_to = "interruption_consequence", values_to = "n") %>%
      mutate(., reweighted_n = if_else(interruption_consequence == "n_missense", n * n_possible_silent / n_possible_missense, as.double(n))),
    aes(y = reweighted_n / n_codons, x = coding_gene_loeuf)
  ) +
    geom_point(aes(color = interruption_consequence), alpha = 0.1) +
    theme_minimal() +
    stat_smooth(aes(color = interruption_consequence)),
  width = 11, height = 7
)


# Autosomal dominant genes?
summary(
  glm(
    n_silent ~ nchar(ru) + fraction_gc + as.factor(coding_gene_is_ad) + offset(log(n_codons)),
    data = 
      coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != "."),
    family = poisson(link = "log")
  )
)
summary(
  glm(
    n_missense ~ nchar(ru) + fraction_gc + as.factor(coding_gene_is_ad) + offset(log(n_codons)),
    data = 
      coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != "."),
    family = poisson(link = "log")
  )
)
# interruptions are not associated one way or another with AD genes
summary(
  glm(
    n ~ nchar(ru) + fraction_gc + as.factor(coding_gene_is_ad) * interruption_consequence + 
      offset(log(n_codons) + log(consequence_offset) - log(n_possible_silent)),
    data = 
      coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      pivot_longer(., cols = c(n_silent, n_missense), names_to = "interruption_consequence", values_to = "n") %>%
      mutate(., consequence_offset = ifelse(interruption_consequence == "n_missense", n_possible_missense, n_possible_silent)),
    family = poisson(link = "log")
  )
)
# how many are there?
coding_strs %>%
  group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
  filter(., maf == max(maf) & 
           purity > 0.75 & 
           nchar(ru) %in% c(3, 6) &
           n_codons != ".") %>%
  group_by(., coding_gene_is_ad) %>%
  count(.)

ggplot(
  data = 
    coding_strs %>%
    group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
    filter(., maf == max(maf) & 
             purity > 0.75 & 
             nchar(ru) %in% c(3, 6) &
             n_codons != "."),
  aes(x = as.factor(nchar(ru)), y = n_silent/n_codons)
) +
  geom_violin(aes(color = as.factor(coding_gene_is_ad)))

# very interesting, they're super GC rich
ggplot(
  data = 
    coding_strs %>%
    group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
    filter(., maf == max(maf) & 
             purity > 0.75 & 
             nchar(ru) %in% c(3, 6) &
             n_codons != "."),
  aes(x = fraction_gc)
) +
  geom_density(aes(color = coding_pfam_domain != "."))

# in protein coding domain?
ggsave(
  "./figures/coding_strs_loeuf_synonymous_vs_missense_in_domain_20241122.pdf",
  ggplot(
    data = 
      coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      pivot_longer(., cols = c(n_silent, n_missense), names_to = "interruption_consequence", values_to = "n") %>%
      mutate(., reweighted_n = if_else(interruption_consequence == "n_missense", n * n_possible_silent / n_possible_missense, as.double(n))),
    aes(y = reweighted_n / n_codons, x = coding_gene_loeuf)
  ) +
    geom_point(aes(color = interruption_consequence), alpha = 0.1) +
    stat_smooth(aes(color = interruption_consequence)) +
    facet_grid(cols = vars(coding_pfam_domain != "."))
)

summary(
  glm(
    n_silent ~ nchar(ru) + fraction_gc + offset(log(n_codons)) + in_domain,
    data = coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      mutate(in_domain = coding_pfam_domain != "."),
    family = poisson(link = "log")
  )
)

summary(
  glm(
    n_silent ~ nchar(ru) + fraction_gc + offset(log(n_codons)) + in_domain * coding_gene_loeuf,
    data = coding_strs %>%
      group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
      filter(., maf == max(maf) & 
               purity > 0.75 & 
               coding_gene_loeuf != -1 & 
               nchar(ru) %in% c(3, 6) &
               n_codons != ".") %>%
      mutate(in_domain = coding_pfam_domain != "."),
    family = poisson(link = "log")
  )
)

coding_strs %>%
  group_by(., chrom, start, end, ru, fraction_gc, coding_gene_loeuf) %>%
  filter(., maf == max(maf) & 
           purity > 0.75 & 
           coding_gene_loeuf != -1 & 
           nchar(ru) %in% c(3, 6) &
           n_codons != ".") %>%
  mutate(in_domain = coding_pfam_domain != ".") %>%
  group_by(., in_domain) %>%
  count()
