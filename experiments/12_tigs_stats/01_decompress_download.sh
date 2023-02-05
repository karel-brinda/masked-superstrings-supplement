#! /usr/bin/env bash

set -e
set -o pipefail
set -u

(
	cd 01_download
	for x in *.xz; do
		xz -vd "$x"
	done
)
