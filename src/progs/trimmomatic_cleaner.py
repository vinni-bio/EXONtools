# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import os
import logging
import csv
import pdb
import shutil

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import parseinput, getinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, run_executor, set_threads
from mains.EXT_errors import EXONtoolsError
from mains.EXT_validator import positive
from utils.sorting import natural_sort


class trimmomatic(EXTprogram):
    """
    Trimming and filtering FASTQ reads, including adapter removal
    """

    name = "trimmomatic"

    def execute_program(self):
        args = self.args
        self.clean_reads(args.inpath, args.outdir, args.forward, args.reverse, args.unpaired, args.adapters, args.select, args.minlen, args.minqual, args.leading, args.trailing, args.cut5end, args.cut3end, args.rqc, args.phred)

    def clean_reads(self, inpath, outdir, forward, reverse, unpaired, adapters, select, minlen, minqual, leading, trailing, cut5end, cut3end, rqc, phred):

        print(self.__dict__)
