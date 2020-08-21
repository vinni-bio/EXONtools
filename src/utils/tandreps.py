# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.


from __future__ import print_function, division
import logging

from mains.EXT_errors import EXONtoolsError


def mischeck(seq1, seq2):
    """Check number of mismatches between two sequences of same length"""

    counter = 0
    for i, j in zip(seq1, seq2):
        if i != j:
            counter += 1
    return counter


def microsat(seq, minlen, maxlen, minrep, minseqlen=0, mismatch=0.0, single=False):
    """Search for tandem repeats in provided sequence
        seq = sequence
        minlen = minimum unit length
        maxlen = maximum unit length
        minrep = minimum number of repeats
        minseqlen = minimal length of STR
        mismatch = percentage of allowed mismatches
        single = allow single nucleotide repeats

    """

    if not isinstance(seq, str):
        logging.error("Sequence must be in string format")
        raise EXONtoolsError("Microsat sequence type error")

    msranges = [x for x in range(minlen, maxlen + 1)]
    results = []
    blocked = []

    seqlen = len(seq)

    for j in msranges:
        for i, s in enumerate(seq):
            STR = ""
            if i + j > seqlen:
                break
            pattern = seq[i:i + j]
            if i not in blocked and (single or not len(set(pattern)) == 1):
                a = i + j
                b = i + j + j
                counter = 1
                while(b <= seqlen and seq[a] == pattern[0] and mischeck(seq[a:b], pattern) / j <= mismatch):
                    STR = STR + seq[a:b]
                    a = a + j
                    b = b + j
                    counter += 1
                if STR and counter >= minrep and len(STR) >= minseqlen:
                    STR = pattern + STR
                    results.append((i + 1, j, counter, pattern, STR))
                    blocked = blocked + [x for x in range(i, a)]
    return results
