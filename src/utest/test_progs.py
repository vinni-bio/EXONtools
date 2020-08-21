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
from progs import abyss
from progs import transabyss
from progs import trinity
from progs import spades

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)


class testPROGS(unittest.TestCase):
    """Testing EXONtools programs"""

    def setUp(self):
        """Runs before each test"""
        self.inpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set2")
        self.outdir = "OUTPUT"
        self.kmers = 'program default'
        self.coverage = 'program default'
        self.npar = 10
        self.mpirun = ""
        self.cores = 1
        self.parallel = 2
        self.memory = 4

    def tearDown(self):
        """Runs after each test"""
        pass

    def test_abyss(self):
        """Testing abyss input parameters"""

        # Testing k-mers
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, [9], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, [], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, ['9'], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, ['proba'], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, ['program default'], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, [1, 10], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, [2.5], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, [-1], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)

        # Testing coverage
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, ['4'], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, [], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, ['4', '10'], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, ['proba'], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, 'proba', self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, ['program default'], self.npar, self.mpirun, self.cores, self.parallel)

        # testing paired coverage
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, self.coverage, 4, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, self.coverage, "10", self.mpirun, self.cores, self.parallel)

        # testing mpi
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, self.coverage, self.npar, "proba", self.cores, self.parallel)

        # testing N cores per job
        self.assertRaises(EXONtoolsError, abyss.abyss_pars, self.kmers, self.coverage, self.npar, self.mpirun, 2, self.parallel)

        # testing function output
        mpi, kmers, coverage = abyss.abyss_pars(self.kmers, self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertFalse(mpi)
        self.assertIsNone(coverage)
        self.assertEqual(kmers, [25])

        # testing mpirun function
        mpi, kmers, coverage = abyss.abyss_pars(self.kmers, self.coverage, self.npar, "echo", self.cores, self.parallel)
        self.assertEqual(str(mpi), " mpirun=echo np=2")
        self.assertNotEqual(str(mpi), " mpirun=echo np=0")

    def test_transabyss(self):
        """Testing transabyss input parameters"""

        # Testing k-mers
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, [9], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, [], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, ['9'], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, ['proba'], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, ['program default'], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, [1, 10], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, [2.5], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, [-1], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)

        # Testing coverage
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, ['4'], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, [], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, ['4', '10'], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, ['proba'], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, 'proba', self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, ['program default'], self.npar, self.mpirun, self.cores, self.parallel)

        # testing paired coverage
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, self.coverage, 4, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, self.coverage, "10", self.mpirun, self.cores, self.parallel)

        # testing mpi
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, self.coverage, self.npar, "proba", self.cores, self.parallel)

        # testing N cores per job
        self.assertRaises(EXONtoolsError, transabyss.transabyss_pars, self.kmers, self.coverage, self.npar, self.mpirun, 2, self.parallel)

        # testing function output
        mpi, kmers, coverage = transabyss.transabyss_pars(self.kmers, self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertEqual(mpi, 0)
        self.assertIsNone(coverage)
        self.assertEqual(kmers, [25])

        # testing mpirun function
        mpi, kmers, coverage = transabyss.transabyss_pars(self.kmers, self.coverage, self.npar, "echo", self.cores, self.parallel)
        self.assertEqual(mpi, 2)
        self.assertNotEqual(mpi, 0)

    def test_trinity(self):
        """Testing trinity input parameters"""

        self.assertRaises(EXONtoolsError, trinity.trinity_pars, ["25"], self.coverage, self.npar, self.mpirun, 2, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, ["24", "25"], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, ["33"], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, ["9"], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, ["proba"], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, [False], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, [25], [1, 2], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, [25], ["1", "2"], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, [25], ["proba"], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, [25], [False], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, [25], self.coverage, "proba", self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, [25], self.coverage, 4.0, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, trinity.trinity_pars, [25], self.coverage, 1, self.mpirun, self.cores, self.parallel)

        params = trinity.trinity_pars(self.kmers, self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertEqual(params, "--min_glue 10")

    def test_spades(self):
        """Testing spades input parameters"""

        self.assertRaises(EXONtoolsError, spades.spades_pars, ["25"], self.coverage, self.npar, self.mpirun, 2, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, ["24", "25"], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, ["128"], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, ["9"], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, ["proba"], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, [False], self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, [25], [1, 2], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, [25], [0], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, [25], ["proba"], self.npar, self.mpirun, self.cores, self.parallel)
        self.assertRaises(EXONtoolsError, spades.spades_pars, [25], [False], self.npar, self.mpirun, self.cores, self.parallel)

        params = spades.spades_pars(self.kmers, self.coverage, self.npar, self.mpirun, self.cores, self.parallel)
        self.assertEqual(params, "--only-assembler")


if __name__ == "__main__":
    unittest.main()
