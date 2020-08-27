# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from __future__ import print_function, division
from mains.EXT_errors import EXONtoolsError
import logging
from utils.seqIO import SeqIO


class phred(object):
    """Phred score numeric sequence"""

    strict = False
    sanger = (33, 73)
    solexa = (59, 104)
    illumina13 = (64, 104)
    illumina15 = (67, 105)
    illumina18 = (33, 74)
    modes = {"S": sanger, "X": solexa, "I": illumina13, "J": illumina15, "L": illumina18}

    def __init__(self, qseq, mode="L"):

        self.mode = mode.upper()

        if isinstance(qseq, str):
            if self.mode in ["S", "L"]:
                adjust = 33
            else:
                adjust = 64
            self.qual = [ord(x) - adjust for x in qseq]
        else:
            logging.error("Phred format error. Please check that quality sequence is a string")
            raise EXONtoolsError("Phred format error")

        if phred.strict and (min(self.qual) + adjust < phred.modes[mode][0] or max(self.qual) + adjust > phred.modes[mode][1]):
            logging.error("Phred format error. Please check your input settings")
            raise EXONtoolsError("Phred format error")

    def average(self):
        try:
            return int(round(sum(self.qual) / len(self.qual), 0))
        except ZeroDivisionError:
            logging.error("Phred values error. Zero length")
            raise EXONtoolsError("Phred format error")

    @staticmethod
    def modecheck(qualseq):
        """Check quality sequence to predict its encoding system"""

        if isinstance(qualseq, str):
            for mode in ["L", "S", "J", "X", "I"]:
                if not (set(qualseq) - set([chr(x) for x in range(phred.modes[mode][0], phred.modes[mode][1] + 1)])):
                    return mode
            logging.error("Unable to determine quality score encoding system")
            raise EXONtoolsError("Phred score error")
        else:
            logging.error("Phred format error. Please check that quality sequence is a string")
            raise EXONtoolsError("Phred format error")

    @classmethod
    def makestrict(cls):
        cls.strict = True


def phredmode(filepath, filetype="FASTQ", nchecks=1000000, error=0.01):
    """Predicts phred encoding for fastq file"""

    checkfile = SeqIO(filepath, filetype)
    modes = []
    count = 0

    for seq in checkfile.read():
        count += 1
        modes.append(phred.modecheck(seq.qual))
        if count > nchecks:
            break

    results = {}
    for mode in ["L", "S", "J", "X", "I"]:
        results[mode] = modes.count(mode)

    maxmode = sorted(results.keys(), key=lambda x: results[x], reverse=True)[0]
    if not SeqIO.dryrun:
        try:
            if results[maxmode] / count > error:
                logging.debug("Phred score is succesfully estimated: '{0:s}'".format(maxmode))
                return maxmode
            else:
                logging.error("Phred prediction is failed due to ambiguous results")
                raise EXONtoolsError("Phred mode prediction error")
        except ZeroDivisionError:
            logging.error("Phred prediction is failed due to empty file provided for evaluation")
            raise EXONtoolsError("Phred mode prediction error")
