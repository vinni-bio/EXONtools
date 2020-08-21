# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import os
import logging

from mains.EXT_IO import makenewdir
from mains.EXT_executor import executor
from mains.EXT_parallel import run_executor


def fastqc_test(indir, outdir, threads, extension="*.fq", filetype="fastq", comment=" "):
    """Runs FastQC analysis for all FASTQ files in the indir folder"""

    logging.debug("Running FastQC analysis")

    FQdir = makenewdir(name=os.path.join(outdir, "FastQC"), fullname="FastQC")

    logging.info("The FastQC results will be saved to '{0:s}' folder".format(os.path.basename(FQdir.path)))

    run_executor(executor(
        program="fastqc",
        params=[threads, FQdir.path, filetype, os.path.join(indir, extension)],
        conditions={"positive": [threads], "pathexists": [indir]},
        custom_arg_string=" > /dev/null 2>&1"
    ))

    try:
        [os.remove(os.path.join(FQdir.path, x)) for x in os.listdir(FQdir.path) if x.endswith(".zip")]
    except IOError:
        pass

    run_executor(executor(
        program="multiqc",
        params=[comment, outdir, "fastqc", FQdir.path],
        conditions={"positive": [threads], "pathexists": [indir]},
        custom_arg_string=" > /dev/null 2>&1"
    ))

    FQdir.delete()

    logging.debug("FastQC analysis completed: OK")
