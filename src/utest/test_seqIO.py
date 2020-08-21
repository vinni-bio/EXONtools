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

from mains.EXT_errors import EXONtoolsError

# TESTING PROGRAMS
from utils.seqIO import SeqIO

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)


class test_seqIO(unittest.TestCase):
    """Testing Fasta processing functions"""

    def setUp(self):
        """Runs before each test"""

        self.inpath1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba1.fasta")
        self.inpath2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba2.fasta")
        self.inpath3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba3.fasta")
        self.inpath4 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba4.fasta")
        self.inpath5 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba5.fasta")
        self.inpath6 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba6.fastaq")
        self.inpath7 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba7.fasta")
        self.inpath8a = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba8_R1.fq")
        self.inpath8b = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba8_R2.fq")
        self.inpath9 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba9.sam")

    def tearDown(self):
        """Runs after each test"""
        pass

    def test_fasta(self):
        """Testing FASTA input"""

        test1 = SeqIO(self.inpath1)
        test2 = SeqIO(self.inpath2)
        test3 = SeqIO(self.inpath3)
        test4 = SeqIO(self.inpath4)
        test5 = SeqIO(self.inpath5)
        test7 = SeqIO(self.inpath7)

        for x in [test1, test2, test3, test4, test5, test7]:
            self.assertEqual(x.type, "fasta")

        test1.readall()
        test3.readall()

        for x in test5.read():
            raise Exception

        for x in test7.read():
            raise Exception

        test5.readall()

        with self.assertRaises(EXONtoolsError):
            test2.readall()

        with self.assertRaises(EXONtoolsError):
            test4.readall()

        with self.assertRaises(EXONtoolsError):
            test4.readall()

        with self.assertRaises(EXONtoolsError):
            test7.readall()

        self.assertRaises(EXONtoolsError, SeqIO, self.inpath6)

        test2.nostrict()
        test7.nostrict()
        test2.readall()
        test7.readall()

        self.assertEqual(len(test1.SEQS), 14)
        self.assertEqual(len(test2.SEQS), 13)
        self.assertEqual(len(test3.SEQS), 14)
        self.assertEqual(len(set(test1.SEQS) ^ set(test3.SEQS)), 0)
        self.assertEqual(len(set(test1.SEQS) ^ set(test2.SEQS)), 1)
        self.assertEqual(sum(map(lambda x: len(x.seq), sum(test1.SEQS.values(), []))), 49479)
        self.assertEqual(sum(map(lambda x: len(x.seq), sum(test2.SEQS.values(), []))), 49479)
        self.assertEqual(sum(map(lambda x: len(x.seq), sum(test3.SEQS.values(), []))), 49479)
        self.assertEqual(len(test5.SEQS), 0)
        self.assertEqual(len(test7.SEQS), 0)

    def test_fastq(self):
        """Testing FASTQ input"""

        self.assertRaises(EXONtoolsError, SeqIO, self.inpath8a)
        self.assertRaises(EXONtoolsError, SeqIO, self.inpath8b)
        test8a = SeqIO(self.inpath8a, "FASTQ")
        test8b = SeqIO(self.inpath8b, "FASTQ")

        test8a.readall()
        test8b.readall()

        self.assertEqual(len(test8a.SEQS), 1000)
        self.assertEqual(len(test8b.SEQS), 1000)

        reads = test8b.callread("61DFRAAXX100204:1:32:12852:17173")[0]
        self.assertEqual(reads.name, "61DFRAAXX100204:1:32:12852:17173")
        self.assertEqual(len(reads.seq), 76)
        self.assertEqual(len(reads.qual), 76)
        self.assertEqual(len(reads.qual), 76)
        self.assertEqual(len(reads.info), 0)
        self.assertEqual(reads.extra, "CHECK")
        self.assertEqual(reads.pair, 2)

        reads = test8b.callread("61DFRAAXX100204:1:32:12852:17174")
        self.assertIsNone(reads)

    def test_sam(self):
        """Testing SAM input"""

        test9 = SeqIO(self.inpath9, "SAM")
        test9.readall()
        reads = sum(test9.SEQS.values(), [])
        self.assertEqual(len(reads), 150)
        self.assertEqual(len([x for x in reads if x.target == "*"]), 50)


if __name__ == "__main__":
    unittest.main()
