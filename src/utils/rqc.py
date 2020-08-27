# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from __future__ import print_function, division
import os
import logging

from mains.EXT_IO import makenewdir
from mains.EXT_executor import executor
from mains.EXT_parallel import run_executor


def rqc_test(indir, outdir, threads, extension="*.fq", filetype="fastq", comment=" "):
    """Runs FastQC analysis for all FASTQ files in the indir directory"""

    logging.debug("Running read quality analysis")

    RQCdir = makenewdir(name=os.path.join(outdir, "RQC"), fullname="RQC")

    logging.info("All read quality results will be saved to '{0:s}' directory".format(indir))

    run_executor(executor(
        program="fastqc",
        params=[threads, RQCdir.path, filetype, os.path.join(indir, extension)],
        conditions={"positive": [threads], "pathexists": [indir]},
        custom_arg_string=" > /dev/null 2>&1"
    ))

    try:
        [os.remove(os.path.join(RQCdir.path, x)) for x in os.listdir(RQCdir.path) if x.endswith(".zip")]
    except IOError:
        pass

    run_executor(executor(
        program="multiqc",
        params=[comment, outdir, "fastqc", RQCdir.path],
        conditions={"positive": [threads], "pathexists": [indir]},
        custom_arg_string=" > /dev/null 2>&1"
    ))

    RQCdir.delete()

    logging.debug("Read quality analysis completed: OK")
