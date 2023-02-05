.PHONY: all help clean cleanall format test test1 cleantest

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

.SUFFIXES:

h=$(shell hostname)
$(info Hostname $h)

all:
	$(MAKE) -C kmercamel
	snakemake -p -j all --rerun-incomplete --scheduler greedy

# Run test
test: cleantest
	$(MAKE) -C kmercamel
	snakemake -p -j 80 --rerun-incomplete --config testing=True

# Run test with 1 thread
test1: cleantest
	$(MAKE) -C kmercamel
	snakemake -p -j 1 --rerun-incomplete --config testing=True

format:
	snakefmt Snakefile

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -frv 10_masked_superstrings/*
	rm -frv 11_maskopt_memtime/$h*
	rm -frv 12_masked_superstrings_properties/$h*
	rm -fr 99_results/*$h*

cleantest:
	rm -fr 99_results/*$h__test.tsv

cleanall: clean ## Clean all
	rm -frv 00_hw_info/$h*
	rm -frv 01_download/$h*

