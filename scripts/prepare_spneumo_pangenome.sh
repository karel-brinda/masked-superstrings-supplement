#! /usr/bin/env bash

set -e
set -o pipefail
set -u

pv ~/github/my/rase-db-spneumoniae-sparc/isolates/*.fa \
	| prophasm -i - -k 32 -o - \
	| seqtk seq -U \
	| xz -9 -T1 \
	> spneumo_pangenome_k32.fa.xz
