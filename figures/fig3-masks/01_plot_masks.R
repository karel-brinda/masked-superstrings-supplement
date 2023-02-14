#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

h <- 10
w <- 16
u <- "cm"

set.seed(42)

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

facet_labeller <-
    labeller(
        k = function(k) {
            paste0("k=", k)
        },
        genome = label_value
    )



# Filter  -----------------------------------------------------------

df <- df0 %>%
    filter(08 <= k & k <= 18) %>%
    filter(!is.na(camel_DS)) %>%
    #filter(M_alg == "default") %>%
    #mutate(d = ifelse(is.na(d), 0, d)) %>%
    mutate(d = paste0(" L", d)) %>%
    mutate(d = str_replace(d, " LNA", "G")) %>%
    mutate(
        M_alg = case_when(
            M_alg == "default" ~ "D",
            M_alg == "ones" ~ "O",
            M_alg == "zeros" ~ "Z",
            M_alg == "runs" ~ "R"
        )
    ) %>%
    filter(S_alg != "streaming") %>%
    mutate(desc = paste0(genome, "-", d))

#df.ht <- df %>% filter(camel_DS == "HT")
#df.ac <- df %>% filter(camel_DS == "AC")



# Plot xz bits per k-mer -----------------------------------------------------------


#df %>% filter(camel_DS == "HT") %>% filter(genome == "spneumoniae-pangenome")

for (g in c("spneumoniae",
            "spneumoniae-pangenome",
            "yeast",
            "sars-cov-2-pangenome")) {
    for (D in c("AC", "HT")) {
        dff <- df %>%
            filter(camel_DS == D) %>%
            filter(genome == g)

        #
        # Mxz_bits_per_superstring_char
        #
        print(paste(g, D, 1))
        ggplot(dff) +
            aes(x = M_alg,
                y = 8 * M_xz_bytes / l,
                fill = M_alg) +
            geom_bar(stat = "identity") +
            scale_x_discrete(name = "") +
            scale_y_continuous(
                name = "",
                breaks = seq(0, 2, 0.2),
                lim = c(0, 0.77),
                expand = c(0, NA)
            ) +
            facet_grid(d ~ k,
                       #, scales = "free") +
                       labeller = facet_labeller) +
            scfill()

        ggsave(
            paste("Mxz_bits_per_superstring_char", g , D, "pdf", sep = "."),
            height = h,
            width = w,
            unit = u
        )

        #
        # M_runs
        #

        print(paste(g, D, 2))
        ggplot(dff %>%
                   filter(camel_DS == D) %>%
                   filter(genome == g)) +
            aes(x = M_alg,
                y = r / l,
                fill = M_alg) +
            geom_bar(stat = "identity") +
            scale_x_discrete(name = "") +
            scale_y_continuous(
                name = "",
                breaks = seq(0, 2, 0.04),
                lim = c(0, NA),
                expand = c(0, NA)
            ) +
            facet_grid(d ~ k,
                       #, scales = "free") +
                       labeller = facet_labeller) +
            scfill()

        ggsave(
            paste("M_runs_per_superstring_char", g , D, "pdf", sep = "."),
            height = h,
            width = w,
            unit = u
        )


    }
}
