# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import absolute_import
import unittest
import sys
import os
import env
import logging

from mains.EXT_IO import output, getinput, parseinput
from mains.EXT_errors import EXONtoolsError

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)


class testIO(unittest.TestCase):
    """Testing IO commands"""

    def setUp(self):
        """Runs before each test"""
        self.t1 = "./OUTPUT"
        self.t2 = os.path.join(os.path.abspath(os.curdir), "OUTPUT")
        self.t3 = "OUTPUT"
        self.t4 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set1")
        self.t5 = os.path.join(self.t4, "ST-1_test.fastq")
        self.t6 = os.path.dirname(self.t4)
        self.t7 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "set2")

    def tearDown(self):
        """Runs after each test"""
        pass

    def test_output(self):
        """Testing the output class methods"""

        output.setdry()

        result1 = output.test_outpath(self.t1)
        result2 = output.test_outpath(self.t2)
        result3 = output.test_outpath(self.t3)

        self.assertEqual(result1, self.t2)
        self.assertEqual(result2, self.t2)
        self.assertEqual(result3, self.t2)

    def test_inpath(self):
        """Testing the input class method - inpath"""

        self.assertTrue(os.path.exists(self.t4))

        getinput.format(".fasta")
        self.assertRaises(EXONtoolsError, getinput.test_inpath, self.t4)
        self.assertRaises(EXONtoolsError, getinput.test_inpath, self.t5)
        self.assertRaises(EXONtoolsError, getinput.test_inpath, self.t6)
        self.assertRaises(EXONtoolsError, getinput.test_inpath, "./abracadabra.fasta")

        getinput.format("default")

        result4 = getinput.test_inpath(self.t4)
        result5 = getinput.test_inpath(self.t5)
        getinput.setdry()
        result6 = getinput.test_inpath(self.t6)
        getinput.undry()

        self.assertIsInstance(result4, list)
        self.assertIsInstance(result5, list)
        self.assertIsInstance(result6, list)

        self.assertEqual(len(result4), 5)
        self.assertIn(self.t5, result4)
        self.assertEqual(len(result5), 1)
        self.assertIn(self.t5, result5)
        self.assertEqual(len(result6), 0)

    def test_paired(self):
        """Testing the input class method - paired"""

        result = getinput(self.t4)
        self.assertRaises(EXONtoolsError, result.find_paired)
        result.setdry()
        result4_1, result4_2 = result.find_paired()
        self.assertIn("ST-1_paired", result4_1)
        self.assertIn("ST-1_unpaired", result4_2)
        self.assertIn("ST-1_unpaired_", result4_2)
        self.assertIn("ST-2_unpaired", result4_2)
        self.assertEqual(len(result4_1), 1)
        self.assertEqual(len(result4_2), 3)
        result.undry()

        result = getinput(self.t7)
        result7_1, result7_2 = result.find_paired()
        self.assertIn("ST-1_paired", result7_1)
        self.assertIn("ST-1_unpaired", result7_2)
        self.assertIn("ST-2_unpaired", result7_2)
        self.assertEqual(len(result7_1), 1)
        self.assertEqual(len(result7_2), 2)

        result.files = [result.files[0]]
        result8_1, result8_2 = result.find_paired()
        self.assertIn("ST-1_unpaired", result8_2)
        self.assertEqual(len(result8_1), 0)
        self.assertEqual(len(result8_2), 1)

        result.files = [os.path.join(self.t4, 'abracadabra_R2.fastq')]
        result9_1, result9_2 = result.find_paired()
        self.assertIn("abracadabra_unpaired", result9_2)
        self.assertEqual(result.files, result9_2["abracadabra_unpaired"])
        self.assertEqual(len(result9_1), 0)
        self.assertEqual(len(result9_2), 1)

        result.files = [os.path.join(self.t4, 'abracadabra.fastq')]
        result10_1, result10_2 = result.find_paired()
        self.assertIn("abracadabra_unpaired", result10_2)
        self.assertEqual(result.files, result10_2["abracadabra_unpaired"])
        self.assertEqual(len(result10_1), 0)
        self.assertEqual(len(result10_2), 1)

    def test_parse_input(self):
        """Testing the input parsing function"""
        self.assertRaises(EXONtoolsError, parseinput, self.t4, "", "")

        result11_1, result11_2 = parseinput(self.t7, "", "")
        self.assertEqual(len(result11_1), 1)
        self.assertEqual(len(result11_1["ST-1_paired"]), 2)
        self.assertEqual(len(result11_2), 2)
        self.assertRaises(EXONtoolsError, parseinput, inpath=self.t7, forward='abracadabra.fastq', reverse="")
        self.assertRaises(EXONtoolsError, parseinput, inpath="", reverse='abracadabra.fastq', forward="")
        self.assertRaises(EXONtoolsError, parseinput, inpath=self.t7, reverse='abracadabra.fastq', forward="")

        result = getinput(self.t7)
        paired, unpaired = result.find_paired()
        result12_1, result12_2 = parseinput(inpath="", forward=paired["ST-1_paired"][0], reverse=paired["ST-1_paired"][1])
        self.assertEqual(len(result12_1), 1)
        self.assertEqual(len(result12_1["ST-1_paired"]), 2)
        self.assertEqual(len(result12_2), 0)

        result13_1, result13_2 = parseinput(inpath="", forward=paired["ST-1_paired"][0], reverse="")
        self.assertEqual(len(result13_1), 0)
        self.assertEqual(len(result13_2["ST-1_unpaired"]), 1)
        self.assertEqual(len(result13_2), 1)
        self.assertEqual(result13_2["ST-1_unpaired"], [paired["ST-1_paired"][0]])

        self.assertRaises(EXONtoolsError, parseinput, inpath="", reverse=paired["ST-1_paired"][1], forward="")


if __name__ == "__main__":
    unittest.main()
