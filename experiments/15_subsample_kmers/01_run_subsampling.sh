#! /usr/bin/env bash

set -e
set -o pipefail
set -u


function subsample {
	k="$1"
	r=0.1
	echo "k=$k"
	echo "r=$r"
	./subsample_kmers.py -k $k -r $r ../../data/spneumo_pangenome_k32.fa.xz \
		| pv -l \
		| xz -9 -T1 \
		> spneumo_pangenome_subsampled_k${k}_r${r}.fa.xz
	echo
}

export -f subsample

printf '12\n14\n16\n18\n20\n' \
	| parallel -j10 subsample
