# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function
import argparse
import sys
import os


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """SETTING THE OUTPUT CLASS FOR ARGPARSE HELP MENUS"""
    pass


class LICENSE(argparse.Action):
    """SETTING THE ACTION CLASS FOR LICENSE OPTION"""

    LCNS_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), "LICENSE.txt")

    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help="Show EXONtools license notice and exit"):
        if os.path.exists(LICENSE.LCNS_PATH):
            with open(LICENSE.LCNS_PATH, 'r') as lcns_file:
                self.lcns_lines = "".join(lcns_file.readlines())
        else:
            self.lcns_lines = "Sorry... The EXONtools license file is missing from its default location."
        super(LICENSE, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)


    def __call__(self, parser, namespace, values, option_string=None):
        formatter = parser._get_formatter()
        formatter.add_text(self.lcns_lines)
        parser._print_message(formatter.format_help(), sys.stdout)
        parser.exit()
