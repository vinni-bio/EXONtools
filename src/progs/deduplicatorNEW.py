# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.


from __future__ import print_function, division
import os
import re
import logging
import csv
import pdb
import gzip
import sys
import shutil

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, parseinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_parallel import run_executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, set_threads
from mains.EXT_errors import EXONtoolsError
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from utils.rqc import rqc_test


class deduplicator(EXTprogram):
    """This command finds and removes reads with PCR duplicates in fastq files leaving
	   just one read with the highest sequencing quality or largest size.
	   PCR duplicates are recognized as identical sequences
	"""

    name = "deduplicator"

    def execute_program(self):
        args = self.args
        self.deduplicate(args.inpath, args.forward, args.reverse, args.unpaired, args.outdir, args.gzoutput, args.skip, args.rqc)

    def deduplicate(self, inpath, forward, reverse, unpaired_in, outdir, gzoutput, skip, rqc):

        debug()

