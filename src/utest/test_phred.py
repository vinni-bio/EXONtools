# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import absolute_import, print_function
import unittest
import sys
import os
from uwarnings import ignore_warnings
import env
import logging

# TESTING PROGRAMS
from mains.EXT_errors import EXONtoolsError
from utils.seqIO import SeqIO
from utils.phred import phred, phredmode

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)


class test_phredIO(unittest.TestCase):
    """Testing Fasta processing functions"""

    def setUp(self):
        """Runs before each test"""

        self.inpath1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba8_R1.fq")
        self.inpath2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba9.sam")
        self.solexa = """;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefgh"""
        self.illumina13 = """@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefgh"""  # included in Solexa
        self.illumina15 = """CDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghi"""

    def tearDown(self):
        """Runs after each test"""
        pass

    def test_phred(self):
        """Testing automatic PHRED score estimation"""

        pmode = phredmode(self.inpath1, "FASTQ")
        self.assertEqual(pmode, "L")

        pmode = phredmode(self.inpath2, "SAM")
        self.assertEqual(pmode, "L")

        phred.makestrict()

        self.assertEqual(phred.modecheck(self.solexa), "X")
        self.assertEqual(phred.modecheck(self.illumina15), "J")
        self.assertEqual(phred.modecheck(self.illumina13), "X")

        fastqfile = SeqIO(self.inpath1, "FASTQ")
        for read in fastqfile.read():
            self.assertRaises(EXONtoolsError, phred, read.qual, "J")
            break


if __name__ == "__main__":
    unittest.main()
