# Masked superstrings – supplementary materials

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
  * [Citation](#citation)
* [Programs](#programs)
  * [Superstring computation - KmerCamel🐫](#superstring-computation---kmercamel)
  * [Mask optimization](#mask-optimization)
* [Experimental evaluation](#experimental-evaluation)
  * [Input data](#input-data)
  * [Reproducing experimental results](#reproducing-experimental-results)
* [Figures](#figures)
  * [Fig. 1](#fig-1)
  * [Fig. 2](#fig-2)
  * [Fig. 3](#fig-3)
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

By default, the superstrings computed by KmerCamel🐫 come with *default* masks;
these contain the minimal possible number of 1's (i.e., every *k*-mer masked on
only once) and the patterns of 1's and 0's reflect the orders in which
individual *k*-mers were added to the superstrings. These default masks are
denoted by `D` in the paper.

Importantly, changes in the underlying data structures (hash-table vs. AC
automaton), as well as changing machines or compilers, results/may result in
different superstrings and their mask, and the specific choices can affect mask
compressibility. For instance, hash-table-based approaches tend to produce more
regular masks that are better compressible (e.g., for nearly complete de Bruijn
graphs).


### Mask optimization

The masks were optimized using the scripts in the
[experiments/08_optimize_masks](experiments/08_optimize_masks/)
directory, which implement individual mask optimization strategies.

* [maskMinNumRuns.py](experiments/08_optimize_masks/maskMinNumRuns.py).
  Minimization of the number of runs of 1's by [integer
  programming](https://en.wikipedia.org/wiki/Integer_programming) using the
  [PuLP solver](https://github.com/coin-or/pulp/), as described in Appendix H.
  (Denoted by `R` in the paper.)
* [maskMaxNumOnes.py](experiments/08_optimize_masks/maskMaxNumOnes.py).
  Maximization of the number of 1's in the mask, which is achieved by
  performing two passes through the superstring. In the first pass, the
  underlying *k*-mer *K* is computed, and in the second pass, all occurrences
  of *k*-mers from *K* are masked on. (Denoted by `O` in the paper.)
* [maskMaxNumZeros.py](experiments/08_optimize_masks/maskMaxNumZeros.py).
  Greedy minimization of the number of 1's (equivalent to maximization of the
  number of 0's), which is achieved by by performing one pass through the data,
  masking on the first occurrence of each *k*-mer while masking off all other
  occurrences. (Denoted by `Z` in the paper.)


## Experimental evaluation

### Input data

The results in Figures 2 and 3 and Tables 1 and 2 were obtained using data in `experiments/11_kmer_camel_comparison_v3/99_results/masked_superstrings_properties.kamenac.tsv`.

Additionally, Tables 1 and 2 use data from `experiments/12_tigs_stats/99_results/masked_superstrings_properties.kamenik.tsv`.

### Reproducing experimental results

After cloning this repository, run the following (besides standard Linux programs, it requires [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [seqtk](https://github.com/lh3/seqtk)):
```
git submodule update --init # to download KmerCamel
cd experiments/11_kmer_camel_comparison_v3/
make test # a fast test that everything works
make
cd ../12_tigs_stats/
make
```


## Figures

Figures are created using Adobe Illustrator and combine individual subfigures generated using R by scripts provided in the respective directories.


### Fig. 1

* [Directory](figures/fig1-overview/)
* [Fig. 1](figures/fig1-overview/fig1.pdf)


### Fig. 2

* [Directory](fig2-camel-comp/)
* [Fig. 2](figures/fig2-camel-comp/fig_camel_comp.pdf)


### Fig. 3

* [Directory](figures/fig3-masks/)
* [Fig. 3](figures/fig3-masks/fig_masks.pdf)



## Remarks

- 11_kmer_camel_comparison_v3: Makefile specifies the number of cores for Snakemake, via argument -j. To test whether everything works, run `make test` first. To make it faster, limit the range of *k* in the Snakefile.
- 12_tigs_stats: uses compressed FASTA files computed by the [matchtigs](https://github.com/algbio/matchtigs) program
