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


def encode_kmer(kmer):  # canonical representation
    #kmer=kmer.upper()
    rc_kmer = rc(kmer)
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
    return canonical_kmer


# Function to maximize # of 1s in the mask of a given superstring
# k: k
# maskedSuperstring: superstring, masked using upper- and lower-case letters; also defines the k-mer set
def maxNumOnes(maskedSuperstring, k):
    K = set()  # set of canonical k-mers

    oldNumOnes = 0
    newNumOnes = 0

    # first pass - collecting k-mers
    for i in range(len(maskedSuperstring) - k + 1):
        if maskedSuperstring[i].isupper():
            # masked as present
            oldNumOnes += 1
            Q = maskedSuperstring[i:i + k].upper()
            Qenc = encode_kmer(Q)
            K.add(Qenc)
    G.info(f"superstring length: {len(maskedSuperstring)}")
    G.info(f"number of kmers: {len(K)}")
    G.info(f"the original mask had {oldNumOnes} ones")

    # second pass - remasking
    nms = []
    for i, x in enumerate(maskedSuperstring):
        Q = maskedSuperstring[i:i + k].upper()
        Qenc = encode_kmer(Q)
        if Qenc in K:
            nms.append(maskedSuperstring[i].upper())
            newNumOnes += 1
        else:
            #this works even for tail, short k-mers aren't in K
            nms.append(maskedSuperstring[i].lower())
    G.info(f"the new mask has {oldNumOnes} ones")

    newMaskedSuperstring = "".join(nms)
    return newMaskedSuperstring


def main():

    parser = argparse.ArgumentParser(
        description=
        "Program for maximizing the number of 1s in a mask of the given superstring."
    )

    parser.add_argument(
        '-k',
        metavar='int',
        dest='k',
        type=int,
        required=True,
        help='kmer size',
    )

    parser.add_argument(
        '-t',
        action='store_true',
        dest='textFile',
        default=False,
        help=
        'input and output are text files with the masked superstring only (input can be split into more lines)',
    )

    parser.add_argument(
        '-p',
        metavar='superstring.fa[.xz,.gz]',
        dest='fs',
        required=True,
        help='FASTA file with superstring (masked by upper/lower-case letters)',
    )

    args = parser.parse_args()

    G.basicConfig(
        level=G.INFO,
        format='[%(asctime)s.%(msecs)03d %(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    # load superstring
    with xopen(args.fs) as fo:
        if args.textFile:  # assuming text file with the masked superstring only (can be split into more lines)
            s = fo.read().replace('\n', '')  # note: assuming UNIX line endings
        else:  # load fasta format
            ss = []
            for qname, seq, _ in readfq(fo):
                G.info(f"Appending superstring {qname}")
                ss.append(seq)
            s = "".join(ss)
    G.info(f"done loading superstring from {args.fs}")
    #G.info(f"superstring = {s}")

    newMaskedSuperstring = maxNumOnes(maskedSuperstring=s, k=args.k)

    if not args.textFile:  # output fasta
        print(f">masked_superstring")
        print(newMaskedSuperstring)
    else:  # text file with a single line
        print(newMaskedSuperstring, end='')  # no new line at the end
    G.info(f"Finished")


if __name__ == "__main__":
    main()
