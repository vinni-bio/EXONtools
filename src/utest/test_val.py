# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, absolute_import
import unittest
import env
import logging
import os

from mains.EXT_validator import memory_check, positive, pathexists
from mains.EXT_errors import EXONtoolsError

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)


class testValidator(unittest.TestCase):
    """Testing IO commands"""

    def setUp(self):
        """Runs before each test"""
        self.t1 = "1"
        self.t2 = None
        self.t3 = 0.999
        self.t4 = 1
        self.t5 = 7
        self.t6 = 8000

    def tearDown(self):
        """Runs after each test.. even when fails"""
        pass

    def test_memoryval(self, msg="Memory tests"):
        """Determines the correct input for memory settings"""

        self.assertRaises(TypeError, memory_check, self.t1)
        self.assertRaises(TypeError, memory_check, self.t2)
        self.assertRaises(EXONtoolsError, memory_check, self.t3)
        self.assertTrue(memory_check(self.t4))
        self.assertTrue(memory_check(self.t5))
        self.assertRaises(EXONtoolsError, memory_check, self.t6)

    def test_positive(self):
        """Checks if the positive validation works correctly"""

        self.assertRaises(EXONtoolsError, positive, [])
        self.assertRaises(EXONtoolsError, positive, [1, -1, 2])
        self.assertRaises(EXONtoolsError, positive, [self.t1])
        self.assertRaises(EXONtoolsError, positive, self.t1)
        self.assertRaises(EXONtoolsError, positive, [self.t2])
        self.assertRaises(EXONtoolsError, positive, self.t2)
        self.assertRaises(EXONtoolsError, positive, self.t3)
        self.assertTrue(positive([self.t3]))
        self.assertRaises(EXONtoolsError, positive, self.t4)
        self.assertTrue(positive([self.t4]))

    def test_pathexists(self):
        """Checks if the path existence check works correctly"""

        self.assertTrue(pathexists, [__file__])
        self.assertTrue(pathexists, [os.path.curdir])
        self.assertRaises(EXONtoolsError, pathexists, __file__)
        self.assertRaises(EXONtoolsError, pathexists, os.path.curdir)
        self.assertRaises(EXONtoolsError, pathexists, "abracadabra")


if __name__ == "__main__":
    unittest.main()
