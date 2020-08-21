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
from mains.EXT_executor import executor

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)


class testEXEC(unittest.TestCase):
    """Testing EXONtools executor commands"""

    def setUp(self):
        """Runs before each test"""
        pass

    def tearDown(self):
        """Runs after each test"""
        pass

    def test_configload(self):
        """Executor: Testing config load method"""

        executor.DEFAULT_CONFIG_FILEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utestdata", "proba.ini")

        self.assertFalse(executor.config_set)
        self.assertEqual(len(executor.program_paths), 0)
        executor.loadConfigs()
        self.assertTrue(executor.config_set)

        self.assertEqual(len(executor.program_paths), 1)
        self.assertIn("uname", executor.program_paths)

    def test_versioncheck(self):
        """Executor: Testing version check method"""

        self.assertTrue(executor.check_version("uname", "0.0.1"))
        self.assertIsNone(executor.check_version("uname", "NA"))
        self.assertRaises(EXONtoolsError, executor.check_version, "uname", "9.9.9.9.9")
        self.assertRaises(EXONtoolsError, executor.check_version, "abracadabra", "0.0.1")

    def test_progcheck(self):
        """Executor: Testing program check method"""
        self.assertIsNone(executor.progcheck("uname"))
        self.assertRaises(EXONtoolsError, executor.progcheck, "pname")
        self.assertRaises(EXONtoolsError, executor.progcheck, "uname", "pname")


if __name__ == "__main__":
    unittest.main()
