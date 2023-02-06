#! /usr/bin/env python3

import argparse
import collections
import functools
import itertools
#import mmh3
import logging as G
import os
import re
import sys
import subprocess

from pulp import LpMaximize, LpMinimize, LpProblem, LpStatus, lpSum, LpVariable, PULP_CBC_CMD  # ILP solver

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


# Function to count the # of runs of 1s in the mask of a given superstring, masked using upper- and lower-case letters
def countNumRuns(maskedSuperstring):
    runs = 0
    lastOne = False
    for a in maskedSuperstring:
        if a.isupper():
            lastOne = True
        else:
            if lastOne:
                runs = runs + 1
            lastOne = False
    if lastOne:
        runs = runs + 1
    return runs


# Function to minimize # of runs of 1s in the mask of a given superstring
# k: k
# maskedSuperstring: superstring, masked using upper- and lower-case letters; also defines the k-mer set
# useMask: if set we only use intervals that contain 1 in the mask
def minNumRuns(k, maskedSuperstring, useMask=False):
    K = {}  # will count occurrences in the superstring
    superstring = maskedSuperstring.upper()  # unmask masked superstring
    # first collect k-mers as keys in K
    for i in range(len(maskedSuperstring) - k + 1):
        if maskedSuperstring[i].isupper():
            Q = superstring[i:i + k]  # computed repeatedly -> cache?
            Qenc = encode_kmer(Q)
            if Qenc not in K:
                K[Qenc] = 0  # only add to dictionary as key
    # find intervals of consecutive non-ghost k-mers
    intervals = []
    currentStart = 0
    intervalContainsOne = False  # does current interval contain 1 in the mask?
    for i in range(len(superstring)):
        Q = superstring[i:i + k].upper()
        Qenc = encode_kmer(Q)
        if Qenc in K:
            if not useMask:  # with useMask, we increase the counter only for intervals that contain 1 in the mask
                K[Qenc] += 1
            intervalContainsOne = True
        else:  # current k-mer is ghost
            if currentStart < i:  # interval ends here
                if not useMask or intervalContainsOne:
                    intervals.append((currentStart, i))
                    if useMask:
                        for j in range(currentStart, i):
                            R = superstring[j:j + k]
                            Renc = encode_kmer(R)
                            K[Renc] += 1
                if useMask and not intervalContainsOne:
                    G.info(
                        f"FOUND interval w/o any 1 in the mask"
                    )  # this is just out of curiosity (in theory, such intervals may exist, but do they appear in the actual output? How often?
            currentStart = i + 1
            intervalContainsOne = False
    G.info(
        f"found {len(intervals)} intervals (these are maximal intervals of possible ones)"
    )
    if useMask:
        G.info(
            f"useMask is set, so these intervals are only maximal extensions of runs of 1s in the mask"
        )
    # after this point, useMask does not make any difference
    numIntervals = 0
    taken = [False for i in range(len(intervals))]
    # go over intervals -- if the current interval contains a k-mer not represented elsewhere (actually, just once in the whole superstring), take it
    for j in range(len(intervals)):
        for i in range(intervals[j][0], intervals[j][1]):
            Q = superstring[i:i + k]
            Qenc = encode_kmer(Q)
            if K[Qenc] == 1:  # TODO: it should be sufficient that the k-mer appears only in the current interval (possibly more than once) and not anywhere else
                taken[j] = True
                break
        if taken[j]:  # we take this interval
            numIntervals = numIntervals + 1
            for i in range(intervals[j][0], intervals[j][1]):
                Q = superstring[i:i + k]
                Qenc = encode_kmer(Q)
                if Qenc in K:
                    K[Qenc] = -1  # k-mer will be represented
    # go over intervals again -- find intervals that are not needed (all of its k-mers are represented in taken intervals)
    numNotneeded = 0
    undecided = []
    remKmers = {
    }  # for each k-mer not in taken intervals, we'll have a list of undecided intervals containing it
    for j in range(len(intervals)):
        if taken[j]:
            continue
        isneeded = False
        for i in range(intervals[j][0], intervals[j][1]):
            Q = superstring[i:i + k]
            Qenc = encode_kmer(Q)
            if K[Qenc] > 0:
                isneeded = True
                break
        if not isneeded:
            numNotneeded += 1
        else:
            undecided.append(j)
            # collect k-mers
            for i in range(intervals[j][0], intervals[j][1]):
                Q = superstring[i:i + k]
                Qenc = encode_kmer(Q)
                if K[Qenc] > 0:  # only k=mers not represented by taken intervals
                    if Qenc not in remKmers:
                        remKmers[Qenc] = []
                    remKmers[Qenc].append(j)
    G.info(
        f"taken {numIntervals} intervals (these are intervals which have a k-mer represented just once)"
    )
    G.info(
        f"may skip {numNotneeded} intervals (all k-mers in these intervals represented in taken intervals)"
    )
    G.info(
        f"remanining {len(intervals) - numIntervals - numNotneeded} intervals (undecided)"
    )
    G.info(
        f"{len(remKmers)} k-mers in undecided intervals that are not represented in taken intervals"
    )
    takeUndecided = 0
    if len(intervals) - numIntervals - numNotneeded > 0:
        if len(remKmers) == 0:
            G.info(
                f"ERROR: no remaming k-mer, but {len(intervals) - numIntervals - numNotneeded} undecided intervals"
            )
        G.info(f"building ILP model")
        model = LpProblem(name="min-num-runs", sense=LpMinimize)
        xvars = {j: LpVariable(name=f"x{j}", cat="Binary") for j in undecided}
        obj_func = lpSum(xvars.values())
        model += obj_func
        for Qenc in remKmers:
            #print(f"k-mer {Qenc} constr {lpSum(xvars[j] for j in remKmers[Qenc]) >= 1}")
            model += (lpSum(xvars[j]
                            for j in remKmers[Qenc]) >= 1, f"k-mer {Qenc}")
        G.info(f"runningt ILP solver")
        status = model.solve(PULP_CBC_CMD(msg=False))
        G.info(f"ILP solver finished, status = {status}")
        #optsol = {}
        for v in model.variables():
            #optsol[v.name] = v.varValue
            taken[int(v.name[1:])] = (v.varValue == 1.0)
            #G.info(f"taking {int(v.name[1:])} with {v.varValue == 1.0}")
        takeUndecided = model.objective.value()
        G.info(f"taking {int(takeUndecided)} undecided intervals"
               )  #, solution = {optsol}")
    G.info(f"constructing optimal mask")
    result = numIntervals + takeUndecided
    newmask = ""
    j = 0  # current interval
    inTakenInterval = False
    for i in range(len(superstring)):
        if j < len(intervals) and i == intervals[j][1]:
            j += 1
        if j >= len(intervals) or i < intervals[j][0] or not taken[j]:
            newmask += superstring[i].lower()
        else:
            newmask += superstring[i].upper()
    if not result == countNumRuns(newmask):  # check for errors
        G.info(
            f"ERROR: result is {result}, while new mask contains {countNumRuns(newmask)} intervals"
        )
    return result, newmask


def main():

    parser = argparse.ArgumentParser(
        description=
        "Program for minimizing the number of runs of 1s in the mask")

    parser.add_argument(
        '-k',
        metavar='int',
        dest='k',
        type=int,
        required=True,
        help='kmer size',
    )

    parser.add_argument(
        '-m',
        action='store_true',
        dest='useMask',
        default=False,
        help=
        'if set, we only use intervals that have a "1" in the mask of the superstring (upper-case letter = 1, lower-case = 0)',
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
    s = ""
    with xopen(args.fs) as fo:
        if args.textFile:  # assuming text file with the masked superstring only (can be split into more lines)
            s = fo.read().replace('\n', '')  # note: assuming UNIX line endings
        else:  # load fasta format
            for qname, seq, _ in readfq(fo):
                G.info(f"Adding superstring {qname}")
                s += seq
    G.info(f"done loading superstring from {args.fs}")
    #G.info(f"superstring = {s}")

    origRuns = countNumRuns(s)
    G.info(f"Number of runs of 1s = {origRuns}")
    res, newmaskedsuperstring = minNumRuns(k=args.k,
                                           maskedSuperstring=s,
                                           useMask=args.useMask)
    G.info(
        f"masked superstring with optimal # of runs of 1s has {int(res)} runs")

    if not args.textFile:  # output fasta
        print(
            f"> masked superstring with optimal # of runs of 1s ({int(res)} runs)"
        )
        print(newmaskedsuperstring)
    else:  # text file with a single line
        print(newmaskedsuperstring, end='')  # no new line at the end
    G.info(f"Finished")


if __name__ == "__main__":
    main()
