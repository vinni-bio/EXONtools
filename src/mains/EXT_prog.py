# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

import logging
from mains.EXT_validator import memory_check


class EXTprogram:
    """Represents a way of completing a EXONtools program using a particular program (with some combination of external
        programs, and custom EXONtools scripts. EXONtools program is responsbile for carrying out their associated
        EXONcommand.
            """
    name = ""
    dryrun = False
    debug = False
    stats = False
    memory = 4
    threads = 1
    keeptmp = False
    suffix = ""
    extra = ""

    def __init__(self, args):
        if args.debugmode:
            EXTprogram.debug = True
        if args.rundry:
            EXTprogram.dryrun = True
        if args.stats:
            EXTprogram.stats = True
        if args.keeptmp:
            EXTprogram.keeptmp = True
        if args.threads > 1:
            EXTprogram.threads = args.threads
        if args.memory:
            memory_check(args.memory, minmem=1, maxpercmem=0.85)
            EXTprogram.memory = args.memory
        if args.extra:
            EXTprogram.extra = args.extra
        if "suffix" in vars(args):
            EXTprogram.suffix = args.suffix
        self.args = args

    @classmethod
    def run_dry(cls, *args):
        if args and cls.dryrun:
            for arg in args:
                arg.setdry()
            logging.debug("Setting up the dry run for selected classes: OK")

    @classmethod
    def set_debug(cls, *args):
        if args and cls.debug:
            for arg in args:
                arg.setdebug()

    def execute_program(self):
        """Completes EXONcommand using the specified program.
        """
        pass
