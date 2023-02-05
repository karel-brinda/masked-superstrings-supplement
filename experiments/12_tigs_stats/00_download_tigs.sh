#! /usr/bin/env bash

set -e
set -o pipefail
set -u

(
	cd 01_download
	curl -L https://iuuk.mff.cuni.cz/~vesely/download/compressed_tigs.tar.xz | tar xvf -
)
