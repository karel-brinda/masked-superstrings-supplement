# Masked superstrings – supplementary materials

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
  * [Citation](#citation)
* [Programs](#programs)
  * [Superstring computation - KmerCamel🐫](#superstring-computation---kmercamel)
  * [Mask optimization](#mask-optimization)
* [Data used for experimental results](#data-used-for-experimental-results)
* [Reproducing experimental results](#reproducing-experimental-results)
* [Remarks](#remarks)

<!-- vim-markdown-toc -->

## Introduction

Here we provide supplementary materials for the paper [Masked superstrings as a unified framework for textual *k*-mer set representations](https://doi.org/10.1101/2023.02.01.526717), including the used data and pipelines.

### Citation

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


## Programs

### Superstring computation - KmerCamel🐫

Superstrings were computed using the
[KmerCamel🐫](https://github.com/GordonHoklinder/kmercamel) program, which
experimentally implements local and global greedy heuristics for masked
superstring computation using hash tables and the [Aho-Corasick
automaton](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm).

The computed superstrings are provided with the default masks provided by
KmerCamel🐫. Such masks contain the minimum possible number of 1's (i.e., every
*k*-mer masked on only once). The specific pattern of 1's and 0's reflect the
orders in which individual *k*-mers are added to the superstrings, therefore,
changes of the underlying data structures (hash-table vs. AC automaton), as
well as changing machines or compilers, results or may result in different
superstrings and masks with different mask compressibility.

The default masks are denoted by `D` in the paper.


### Mask optimization

The masks were optimized using scripts in the
[experiments/08_optimize_masks](experiments/08_optimize_masks/)
directory, which implement individual mask optimization strategies.

* [maskMinNumRuns.py](experiments/08_optimize_masks/maskMinNumRuns.py).
  Minimization of the number of runs of 1's using [integer
  programming](https://en.wikipedia.org/wiki/Integer_programming) using the
  [PuLP solver](https://github.com/coin-or/pulp/), as described in Appendix H.
  (Denoted by `R` in the paper.)
* [maskMaxNumOnes.py](experiments/08_optimize_masks/maskMaxNumOnes.py).
  Maximization of the number of 1's in the mask, which is done by 2 passes
  through the superstring - collection of *k*-mers and masking on all of their
  occurrences. (Denoted by `O` in the paper.)
* [maskMaxNumZeros.py](experiments/08_optimize_masks/maskMaxNumZeros.py).
  Greedy minimization of the number of 1's (maximization of the number of 0's),
  which is done by 1 pass through the data, masking on the first occurrence of
  each *k*-mer and masking off the following ones. (Denoted by `Z` in the paper.)


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
