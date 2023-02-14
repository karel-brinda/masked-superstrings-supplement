#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

h <- 7
w <- 20
u <- "cm"


theme_set(theme_bw() +
              theme(
                  axis.text.x = element_text(
                      size = 8,
                      angle = 0,
                      hjust = 0.5,
                      vjust = 0.2

                  ),
                  panel.grid.major.x = element_blank()
              ))

scfill <- scale_fill_npg


# Longer format -----------------------------------------------------------

df0 <- read_tsv("all_data.tsv", na = "na")

df1 <- df0 %>%
    pivot_longer(c("S_xz_bytes",
                   "M_xz_bytes"),
                 names_to = "key",
                 values_to = "bytes") %>%
    mutate(key = str_replace(key, "_xz_bytes", "")) %>%
    mutate(mb = bytes / 1e6)

facet_labeller <-
    labeller(
        k = function(k) {
            paste0("k=", k)
        },
        genome = label_value
    )



# Filter  -----------------------------------------------------------

df2 <- df1 %>%
    filter(08 <= k & k <= 18) %>%
    filter(M_alg == "default") %>%
    #mutate(d = ifelse(is.na(d), 0, d)) %>%
    mutate(d = ifelse(is.na(d), 7, d)) %>%
    filter(! is.na(camel_DS)) %>%
    filter(S_alg != "streaming") %>%
    filter(
        genome == "spneumoniae-pangenome" |
            genome == "spneumoniae" #| genome == "sars-cov-2-pangenome"
    ) %>%
    mutate(genome = str_replace(genome, "spneumoniae-pangenome", "pan-genome")) %>%
    mutate(genome = str_replace(genome, "spneumoniae", "genome"))

df.ht <- df2 %>% filter(!grepl('AC', S_alg))
df.ac <- df2 %>% filter(grepl('AC', S_alg))




# Plot xz bits per k-mer -----------------------------------------------------------

# trick - copies of the same tick disappear <- using the same
# name can be used for tick removal
ggplot(df.ht) +
    aes(x = d,
        y = 8 * bytes / kmer_count,
        fill = key) +
    geom_bar(stat = "identity") +
    #, scales = "free"
    facet_grid(genome ~ k, labeller = facet_labeller, scales = "free") +
    scfill(labels = c("mask", "superstring")) +
    scale_x_discrete(name = "",
                     limits = c("L1 ", " .", " . ", ". ", "5", " . ", "G")) +
    scale_y_continuous(expand = c(0, 0),
                       #name = "xz bits per distinct k-mers",
                       name = "") +
	   theme(legend.title= element_blank())


ggsave(
    "xz_bits_per_kmer_HT.pdf",
    height = h,
    width = w,
    unit = u
)




ggplot(df.ac) +
    aes(x = d,
        y = 8 * bytes / kmer_count,
        fill = key) +
    geom_bar(stat = "identity") +
    #, scales = "free"
    facet_grid(genome ~ k, labeller = facet_labeller) +
    scfill() +
    scale_x_discrete(name = "",
                     limits = c("L1 ", " .", " . ", ". ", "5", " . ", "G")) +
    scale_y_continuous(expand = c(0, 0),
                       #name = "xz bits per distinct k-mers",
                       name = "") +
	   theme(legend.title= element_blank())


ggsave(
    "xz_bits_per_kmer_AC.pdf",
    height = h,
    width = w,
    unit = u
)







# Plot superstring chars per k-mer -----------------------------------------------------------

ggplot(df.ht %>% filter(key == "S")) +
    aes(x = d,
        y = l / kmer_count,
        shape = algorithm) +
    scfill() +
    geom_point() +
    geom_line(aes(group = algorithm)) +
    facet_grid(genome ~ k, labeller = facet_labeller) +
    scale_x_discrete(name = "",
                     limits = c("L1 ", " .", " . ", ". ", "5", " . ", "G")) +
    scale_y_continuous(
        #name = "characters per distinct k-mer",
        name = "",
        breaks = seq(1, 2, 0.2),
        lim = c(1.0, 2.05),
        expand = c(0, 0)
    ) +
    theme(plot.margin = margin(0, 3.3, 0, 0, "cm")) +
    guides(shape = "none")

ggsave(
    "chars_per_kmer_HT.pdf",
    height = h,
    width = w,
    unit = u
)




ggplot(df.ac %>% filter(key == "S")) +
    aes(x = d,
        y = l / kmer_count,
        shape = algorithm) +
    scfill() +
    geom_point() +
    geom_line(aes(group = algorithm)) +
    facet_grid(genome ~ k, labeller = facet_labeller) +
    scale_x_discrete(#name = "d",
        name = "",
        limits = c("L1 ", " .", " . ", ". ", "5", " . ", "G")) +
    scale_y_continuous(
        #name = "characters per distinct k-mer",
        name = "",
        breaks = seq(1, 2, 0.2),
        lim = c(1.0, 2.05),
        expand = c(0, 0)
    ) +
    theme(plot.margin = margin(0, 3.3, 0, 0, "cm")) +
    guides(shape = "none")

ggsave(
    "chars_per_kmer_AC.pdf",
    height = h,
    width = w,
    unit = u
)
