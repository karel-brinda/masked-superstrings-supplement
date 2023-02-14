# Masked superstrings – supplementary materials

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Citation](#citation)
* [KmerCamel🐫](#kmercamel)
* [Data used for experimental results](#data-used-for-experimental-results)
* [Reproducing experimental results](#reproducing-experimental-results)
* [Remarks](#remarks)

<!-- vim-markdown-toc -->

## Introduction

Here we provide supplementary materials for the paper [Masked superstrings as a unified framework for textual *k*-mer set representations](https://doi.org/10.1101/2023.02.01.526717), including the used data and pipelines.

## Citation

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Masked superstrings as a unified framework for textual *k*-mer set representations. *bioRxiv* 2023.02.01.526717, 2023.
[https://doi.org/10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

```
@article{sladky2023-masked-superstrings,
  title   = { Masked superstrings as a unified framework for textual $k$-mer set representations },
  author  = { Sladk{\'y}, Ond{\v r}ej and Vesel{\'y}, Pavel and B{\v r}inda, Karel },
  journal = { bioRxiv },
  volume  = { 2023.02.01.526717 },
  year    = { 2023 },
  doi     = { 10.1101/2023.02.01.526717 }
}
```

## KmerCamel🐫

The analyses were performed using [KmerCamel🐫](https://github.com/GordonHoklinder/kmercamel),
which experimentally implements local and global greedy heuristics for masked superstring
computation using hash tables and the Aho-Corasick automaton.

## Data used for experimental results
The results in Figures 2 and 3 and Tables 1 and 2 were obtained using data in `experiments/11_kmer_camel_comparison_v3/99_results/masked_superstrings_properties.kamenac.tsv`.

Additionally, Tables 1 and 2 use data from `experiments/12_tigs_stats/99_results/masked_superstrings_properties.kamenik.tsv`.

## Reproducing experimental results
After cloning this repository, run the following (besides standard Linux programs, it requires [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [seqtk](https://github.com/lh3/seqtk)):
```
git submodule update --init # to download KmerCamel
cd experiments/11_kmer_camel_comparison_v3/
make test # a fast test that everything works
make
cd ../12_tigs_stats/
make
```

## Remarks

- 11_kmer_camel_comparison_v3: Makefile specifies the number of cores for Snakemake, via argument -j. To test whether everything works, run `make test` first. To make it faster, limit the range of *k* in the Snakefile.
- 12_tigs_stats: uses compressed FASTA files computed by the [matchtigs](https://github.com/algbio/matchtigs) program
