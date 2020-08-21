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
from progs.assasses import metrics
from utils.gccheck import GC, Ncheck

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)


class testASSASSES(unittest.TestCase):
    """Testing assembly assessment class"""

    def setUp(self):
        """QUAST results"""

        self.asspath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set3", "proba1.fasta")

        self.lib = "proba1"
        self.ID = "MMMM"
        self.totnum = 14
        self.long500 = 14
        self.long1000 = 10
        self.long5000 = 5
        self.long10000 = 0
        self.max = 8724
        self.min = 503
        self.totlen = 49479
        self.N25 = 7144
        self.N50 = 6567
        self.N75 = 3880
        self.L50 = 4
        self.L75 = 6
        self.mean = 3534
        self.esize = 5733
        self.NG50 = 6567
        self.nbases = 201
        self.nbas100 = 406
        self.GC = 0.49

    def tearDown(self):
        """Runs after each test"""
        pass

    def test_metrics(self):
        """Testing assembly quality metrics"""

        self.assertIsNone(metrics.NG50)

        asseval = metrics(self.asspath)

        self.assertEqual(metrics.NG50, self.totlen)
        self.assertEqual(asseval.NG50, self.totlen)
        self.assertEqual(asseval.NG, self.N50)
        self.assertEqual(asseval.N50, self.N50)
        self.assertEqual(asseval.totnum, self.totnum)
        self.assertEqual(asseval.longest, self.ID)
        self.assertEqual(asseval.ID, self.lib)
        self.assertEqual(asseval.long500, self.long500)
        self.assertEqual(asseval.long1000, self.long1000)
        self.assertEqual(asseval.long5000, self.long5000)
        self.assertEqual(asseval.long10000, self.long10000)
        self.assertEqual(asseval.esize, self.esize)
        self.assertEqual(asseval.mean, self.mean)
        self.assertEqual(asseval.N25, self.N25)
        self.assertEqual(asseval.N50, self.N50)
        self.assertEqual(asseval.N75, self.N75)
        self.assertEqual(asseval.L50, self.L50)
        self.assertEqual(asseval.L75, self.L75)
        self.assertEqual(asseval.nbases, self.nbases)
        self.assertEqual(asseval.n100, self.nbas100)
        self.assertEqual(asseval.gc, self.GC)

    def test_Ncheck(self):
        """Testing N base count in sequences"""

        count = 0

        with open(self.asspath, 'r') as infile:
            for line in infile:
                if not line.startswith(">"):
                    count += Ncheck(line)

        self.assertEqual(count, self.nbases)
        self.assertEqual(count * 100000 // self.totlen, self.nbas100)

    def test_GCcont(self):
        """Testing GC assessment in sequences"""

        GCcont = []

        with open(self.asspath, 'r') as infile:
            for line in infile:
                if not line.startswith(">"):
                    GCcont.append(GC(line))

        self.assertEqual(round(sum(GCcont) / len(GCcont), 2), self.GC)


if __name__ == "__main__":
    unittest.main()
