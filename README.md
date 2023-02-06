# Supplementary repository for paper on masked superstrings

Paper preprint at bioRxiv: 
  
&nbsp;&nbsp;&nbsp;&nbsp;Ondřej Sladký, Pavel Veselý, and Karel Břinda: *Masked superstrings as a unified framework for textual 𝑘-mer set representations.* 
[doi:10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

The [KmerCamel🐫](https://github.com/GordonHoklinder/kmercamel) program provides an implementation of the heuristics proposed in the paper (note: this is still an experimental implementation and can be substantially optimized w.r.t. running time or memory requirements).
  
**Particular data used for experimental results:**
The results in Figures 2 and 3 and Tables 1 and 2 were obtained using data in `experiments/11_kmer_camel_comparison_v3/99_results/masked_superstrings_properties.kamenac.tsv`.

Additionally, Tables 1 and 2 use data from `experiments/12_tigs_stats/99_results/masked_superstrings_properties.kamenik.tsv`.

**Reproducing experimental results:**
After cloning this repository, run the following (besides standard Linux programs, it requires [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [seqtk](https://github.com/lh3/seqtk)):
```
git submodule update --init # to download KmerCamel
cd experiments/11_kmer_camel_comparison_v3/
make test # a fast test that everything works
make
cd ../12_tigs_stats/
make
```

Notes:
- 11_kmer_camel_comparison_v3: Makefile specifies the number of cores for Snakemake, via argument -j. To test whether everything works, run `make test` first. To make it faster, limit the range of *k* in the Snakefile.
- 12_tigs_stats: uses compressed FASTA files computed by the [matchtigs](https://github.com/algbio/matchtigs) program
