#!/usr/bin/env Rscript

library(tidyverse)


# merge masked superstring stats ---------------------------------------------------------------

df.camel.ms <- read_tsv("input/camel.masked_superstrings.tsv")
df.tigs.ms <- read_tsv("input/tigs.masked_superstrings.tsv")

df.ms0 <- df.tigs.ms %>%
    bind_rows(df.camel.ms) %>%
    mutate(genome = str_replace(genome, "spneumo-pangenome", "spneumoniae-pangenome"))


# camel stats ---------------------------------------------------------------

df.ms <- df.ms0 %>%
    mutate(pref0 = str_replace(pref, '(.+)\\..+', '\\1')) %>%
    mutate(
        algorithm = case_when(
            grepl(regex("^pseudosimplitigs(|AC)$"), S_alg) ~ "loc-greedy",
            grepl(regex("^greedy(|AC)$"), S_alg) ~ "glob-greedy",
            grepl("streaming", S_alg) ~ "streaming",
            TRUE ~ S_alg
        )
    ) %>%
    mutate(`camel_DS` = case_when(grepl(
        regex("^(greedy|pseudosimplitigs)AC$"), S_alg
    ) ~ "AC",
    grepl(
        regex("^(greedy|pseudosimplitigs)$"), S_alg
    ) ~ "HT",
    TRUE ~ "na")) %>%
    relocate("genome", "k", "S_alg", "algorithm", "camel_DS")

df.mem0 <- read_tsv("input/camel.memtime.tsv")
df.mem <- df.mem0 %>%
    transmute(
        pref0 = pref,
        camel_mem_bytes = `max_RAM(kb)` * 1000,
        camel_real_s = `real(s)`,
        camel_usr_s = `user(s)`,
        camel_sys_s = `sys(s)`,
    )


df <- df.ms %>%
    full_join(df.mem) %>%
    relocate(pref0, pref, .after = last_col()) %>%
    arrange(genome, k, algorithm)

df %>%
    write_tsv("all_data.tsv",  na = "na")
