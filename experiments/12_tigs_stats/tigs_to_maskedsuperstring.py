#! /usr/bin/env python3

import argparse
import collections
import functools
import itertools
import logging as G
import os
import re
import sys
import subprocess

from pathlib import Path
from xopen import xopen
from pprint import pprint

c = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


def rc(kmer):
    rc_kmer = "".join([c[x] for x in kmer[::-1]])
    return rc_kmer


def readfq(fp):  # this is a generator function
    # From https://github.com/lh3/readfq/blob/master/readfq.py
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last: break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def tigs_to_masked_superstring(fn, k):
    sm = []
    with xopen(fn) as fo:
        for name, seq, _ in readfq(fo):
            s = seq[:-k + 1].upper() + seq[-k + 1:].lower()
            sm.append(s)
    return "".join(sm)


def main():

    parser = argparse.ArgumentParser(
        description="TIGS SPSS to masked superstring")

    parser.add_argument(
        '-k',
        metavar='int',
        dest='k',
        type=int,
        required=True,
        help='kmer size',
    )

    parser.add_argument(
        '-p',
        metavar='superstring.fa[.xz,.gz]',
        dest='fn',
        required=True,
        help='FASTA file with tig sequences',
    )

    args = parser.parse_args()

    G.basicConfig(
        level=G.INFO,
        format='[%(asctime)s.%(msecs)03d %(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    newMaskedSuperstring = tigs_to_masked_superstring(args.fn, k=args.k)
    print(">masked_superstring")
    print(newMaskedSuperstring)

    G.info(f"Finished")


if __name__ == "__main__":
    main()
