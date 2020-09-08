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

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        deduplicator.run_dry(SeqIO, getinput, output, makenewdir, executor)
        deduplicator.set_debug(executor)
        if rqc:
            executor.setconfig("fastqc", "multiqc")

        # IO SETTINGS
        paired, unpaired = fixunpaired(inpath, forward, reverse, unpaired_in)
        librs = findlibrs(paired, unpaired)
        output(outdir)
        tmpdir = makenewdir(name="tmp", fullname="temporary")
        logging.debug("IO settings: OK")
        debug()

        # SET DEDUP PARAMETERS
        deduppars(skip, gzoutput)
        
        TASKS = []
        if paired:
            for lib in sorted(paired.keys(), key=natural_sort):
                TASKS.append(worker(dedup, [lib, tmpdir.path, skip, suffix, paired[lib][0], paired[lib][1]]))
        if unpaired:
            for lib in sorted(unpaired.keys(), key=natural_sort):
                TASKS.append(worker(dedup, [lib, tmpdir.path, skip, suffix, paired[lib][0], None]))

def dedup():
    """deduplication"""
    pass

def debug():
    """debuger"""
    if deduplicator.debug:
        pdb.set_trace()


def savestats(results, outpath):
    """Save all stats"""
    if deduplicator.stats and not deduplicator.dryrun:
        statdir = makenewdir(name=os.path.join(outpath, "STATS"), fullname="STATS")
        logging.info("Read stats will be saved to 'STATS/correct_stats.csv'")
        header = ["No", "FILE", "LIBRARY", "#READS", "#TRIMMED", "#CORRECTED"]
        with open(os.path.join(statdir.path, "correct_stats.csv"), 'w') as statfile:
            csv_writer = csv.writer(statfile)
            csv_writer.writerow(header)
            counter = 0
            for lib in sorted(results.keys(), key=natural_sort):
                for direct in ['forward', 'reverse', 'unpaired']:
                    if results[lib][direct]["outpath"]:
                        counter += 1
                        csv_writer.writerow([
                            counter,
                            os.path.basename(results[lib][direct]["outpath"]),
                            lib,
                            results[lib][direct]["counts"],
                            results[lib][direct]["trimmed"],
                            results[lib][direct]["changed"]
                        ])
        logging.debug("Read stats were successfully written to the file: OK")
        debug()
        return statdir


def deduppars(skip, gzoutput):
    """set dedup parameters"""
    logging.debug("Setting up DEDUP parameters")
    if deduplicator.extra:
        logging.warning("Extra string argument cannot be used in 'deduplicator' program and therefore will be omitted")
    if gzoutput:
        logging.warning("All output FASTQ files will be compressed with gzip")
    direction = ['forward' , 'forward', 'reverse', 'reverse']
    for i,x in enumerate(skip):
        if x != 0:
            logging.info("The first {0:s} bases will be masked in {1:s} reads".format(x,direction[i]))
    logging.debug("DEDUP parameters: OK")
    return SPADES_params


def runqc(rqc, outpath, message, gzout):
    """Run QC analysis"""
    if rqc and not deduplicator.stats and not deduplicator.dryrun:
        statdir = makenewdir(name=os.path.join(outpath, "STATS"), fullname="STATS").path
    else:
        statdir = os.path.join(outpath, "STATS")
    if rqc and gzout and not deduplicator.dryrun:
        rqc_test(outpath, statdir, deduplicator.threads, extension="*.gz", comment=message)
        debug()
    elif rqc and not deduplicator.dryrun:
        rqc_test(outpath, statdir, deduplicator.threads, comment=message)
        debug()
    else:
        pass


def findlibrs(paired, unpaired):
    """Parse IO to find librs"""
    librs = list(set([x.split("_")[0] for x in list(paired.keys()) + list(unpaired.keys())]))
    librs.sort(key=natural_sort)
    return librs


def makelibtemp(sample, tmppath):
    """make temporary file for provided library"""
    libouttmp = os.path.join(tmppath, sample + "_hammer")
    makenewdir(libouttmp)
    logging.debug("Temporary directory: OK")
    return libouttmp


def fixunpaired(inpath, forward, reverse, unpaired):
    """parse input files"""
    logging.debug("Parsing input files")
    getinput.format(['.fq', '.fastq', '.fastq.gz', '.fq.gz'])
    paired, newunpaired = parseinput(inpath, forward, reverse)
    if unpaired and not inpath:
        newunpaired = {os.path.basename(unpaired).split("_")[0].split(".")[0] + "_unpaired": getinput(unpaired).files}
        logging.debug("Input file parsing: OK")
        return paired, newunpaired
    elif unpaired and inpath:
        logging.error("Input arguments '-i' and '-U' cannot be used simultaneously")
        raise EXONtoolsError("Input path error")
    else:
        logging.debug("Input file parsing: OK")
        return paired, newunpaired


def count_reads(inpath):
    """Count reads by number of lines"""
    logging.debug("Count reads by line")
    infileseq = SeqIO(inpath, fileformat="FASTQ")
    infileseq.totalcount(lines=True)
    if infileseq.total % 4 != 0:
        logging.error("FASTQ format error. Please check your input file.")
        raise EXONtoolsError("Wrong fastq format")
    logging.debug("Count reads by line: OK")
    return infileseq.total // 4


def truncreads(skip, inpathR2):
    """skip conformation"""
    logging.debug("Changing skip vector")
    trunc = skip[:]
    for i, x in enumerate(trunc):
        if x == 0:
            trunc[i] = None
        elif i > 1 and not inpathR2:
            trunc[i] = None
        elif i % 2 == 1:
            trunc[i] = -x
    logging.debug("Changing skip vector: OK")
    return trunc[0], trunc[1], trunc[2], trunc[3]
