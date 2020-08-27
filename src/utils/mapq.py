# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2020

# The following mapping quality is similar to Bowtie2/Tophat
# function for single reads with global alignment
# The code was originally proposed here:
# https://github.com/BenLangmead/bowtie2/blob/master/unique.h
# under GNU General Public License v3.0

# Here it is used to correct '255' mapping scores in corrected sam files

from __future__ import print_function, division
import logging
from mains.EXT_errors import EXONtoolsError


def MAPQ(AS, XS, seqlen):
    """Calculates aka bowtie2 mapping quality score for a single read"""

    if not isinstance(AS, int):
        logging.error("Alignment score (AS) should have an integer format")
        raise EXONtoolsError("MAPQ error")

    if not isinstance(seqlen, int):
        logging.error("Sequence length should have an integer format")
        raise EXONtoolsError("MAPQ error")

    scmin = -0.6 - 0.6 * (seqlen)

    if not XS:
        XS = scmin - 1
    if XS > AS:
        return 0

    diff = abs(scmin)
    bestover = AS - scmin
    bestdiff = abs(abs(AS) - abs(XS))
    if XS < scmin:
        if bestover >= diff * 0.8:
            return 42
        elif bestover >= diff * 0.7:
            return 40
        elif bestover >= diff * 0.6:
            return 24
        elif bestover >= diff * 0.5:
            return 23
        elif bestover >= diff * 0.4:
            return 8
        elif bestover >= diff * 0.3:
            return 3
        else:
            return 0
    else:
        if bestdiff >= diff * 0.9:
            if bestover == diff:
                return 39
            else:
                return 33
        elif bestdiff >= diff * 0.8:
            if bestover == diff:
                return 38
            else:
                return 27
        elif bestdiff >= diff * 0.7:
            if bestover == diff:
                return 37
            else:
                return 26
        elif bestdiff >= diff * 0.6:
            if bestover == diff:
                return 36
            else:
                return 22
        elif bestdiff >= diff * 0.5:
            if bestover == diff:
                return 35
            elif bestover >= diff * 0.84:
                return 25
            elif bestover >= diff * 0.68:
                return 16
            else:
                return 5
        elif bestdiff >= diff * 0.4:
            if bestover == diff:
                return 34
            elif bestover >= diff * 0.84:
                return 21
            elif bestover >= diff * 0.68:
                return 14
            else:
                return 4
        elif bestdiff >= diff * 0.3:
            if bestover == diff:
                return 32
            elif bestover >= diff * 0.88:
                return 18
            elif bestover >= diff * 0.67:
                return 15
            else:
                return 3
        elif bestdiff >= diff * 0.2:
            if bestover == diff:
                return 31
            elif bestover >= diff * 0.88:
                return 17
            elif bestover >= diff * 0.67:
                return 11
            else:
                return 0
        elif bestdiff >= diff * 0.1:
            if bestover == diff:
                return 30
            elif bestover >= diff * 0.88:
                return 12
            elif bestover >= diff * 0.67:
                return 7
            else:
                return 0
        elif bestdiff > 0:
            if bestover >= diff * 0.67:
                return 6
            else:
                return 2
        else:
            if bestover >= diff * 0.67:
                return 1
            else:
                return 0
