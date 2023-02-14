#! /usr/bin/env bash

set -e
set -o pipefail
set -u

curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz \
	| seqtk seq -U \
	| xz -9 -T1 \
	> yeast.fa.xz
