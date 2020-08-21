# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import os
import re
import logging
import csv
import pdb
import gzip

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, parseinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, set_threads
from mains.EXT_errors import EXONtoolsError
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from utils.fastqc import fastqc_test


class seqformatter(EXTprogram):
    """This program verifies names and structure of raw reads"""

    name = "seqformatter"

    def execute_program(self):
        args = self.args
        self.format_seqs(args.inpath, args.forward, args.reverse, args.unpaired, args.outdir, args.gzoutput, args.rename, args.pattern, args.customgrep, args.skipcheck, args.fastqc)

    def format_seqs(self, inpath, forward, reverse, unpaired_in, outdir, gzoutput, rename, pattern, customgrep, skipcheck, fastqc):

        if seqformatter.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        seqformatter.run_dry(SeqIO, getinput, output, makenewdir, executor)
        seqformatter.set_debug(executor)

        if fastqc:
            executor.setconfig("fastqc")

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        seqformatter.run_dry(getinput, output, makenewdir, executor)
        seqformatter.set_debug(executor)

        if seqformatter.debug:
            pdb.set_trace()

        if seqformatter.extra:
            logging.warning("Extra string argument cannot be used in 'seqformatter' program and therefore will be ommited")

        if skipcheck:
            logging.warning("Skipping FASTQ read identifier checkup... Only parsing output will be performed")

        if customgrep and pattern.lower() != "custom":
            logging.error("You can use any custom grep pattern only after you select '--type CUSTOM' option")
            raise EXONtoolsError("Missing custom pattern option")
        elif not customgrep and pattern.lower() == "custom":
            logging.error("Provide custom grep pattern using '--custom' option")
            raise EXONtoolsError("Missing custom pattern option")
        else:
            if pattern.lower() == "illumina":
                logging.info("Reads will be verified for conforming the Illumina raw read format")
                patterncheck = "^(.*) [1|2]*:N:\d*:.*$"
            elif pattern.lower() == "exontools":
                logging.info("Reads will be verified for conforming the EXONtools read name format")
                patterncheck = "^([^\s|_]+_[un]*paired_ID\d+)/\d$"
            else:
                patterncheck = customgrep

        # GET READ INPUT FILES
        getinput.format(['.fq', '.fastq', '.fastq.gz', '.fq.gz'])
        paired, unpaired = parseinput(inpath, forward, reverse)
        if unpaired_in and not inpath:
            unpaired = {os.path.basename(unpaired_in).split("_")[0].split(".")[0] + "_unpaired": getinput(unpaired_in).files}
        elif unpaired_in and inpath:
            logging.error("Input arguments '-i' and '-U' cannot be used simultaneously")
            raise EXONtoolsError("Input path error")
        else:
            pass

        # checks if reads should be renamed and informs user
        if rename:
            logging.warning("All reads in FASTQ files will be renamed according to library identifiers (--rename)")

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        if gzoutput:
            logging.warning("All output FASTQ files will be compressed with gzip")

        # # MAKE TMP DIRECTORIES
        # tmpdir = makenewdir(name="tmp", fullname="temporary")

        logging.debug("Test and IO settings: OK")

        if seqformatter.debug:
            pdb.set_trace()

        # Collects all tasks
        TASKS = []
        if rename:
            datnames = makenewdir(name="CATALOG", fullname="CATALOG")
        else:
            datnames = None
        if paired:
            for lib in sorted(paired.keys(), key=natural_sort):
                TASKS.append(worker(validate_reads, [paired[lib][0], paired[lib][1], output.path, lib, patterncheck, datnames, rename, gzoutput, skipcheck]))
        if unpaired:
            for lib in sorted(unpaired.keys(), key=natural_sort):
                TASKS.append(worker(validate_reads, [unpaired[lib][0], None, output.path, lib, patterncheck, datnames, rename, gzoutput, skipcheck]))

        logging.info("Running the FASTQ read checkup analysis for all files")
        if TASKS:
            processes_requested = set_threads("FASTQ file formatter", len(TASKS), seqformatter.threads)
            pool = create_pool(processes_requested)
            jobs = hard_worker(run_instance, TASKS, pool)
            close_pool(pool)

        stats_collector = {}
        if jobs:
            for result in jobs:
                stats_collector.update(result)
            del jobs
        else:
            logging.error("Multiprocessing error in FASTQ file analysis")
            raise EXONtoolsError("Multiprocessing error in FASTQ file analysis")
        logging.debug("Checkup analysis succesfully finished: OK")

        if seqformatter.debug:
            pdb.set_trace()

        # SAVES STATS TO *.CSV TABLE
        if seqformatter.stats and not seqformatter.dryrun:
            logging.info("Read stats will be saved to 'read_stats.csv'")
            header = ["No", "Library", "BARCODE", "#READS", "#FILTERED"]
            with open(os.path.join(output.path, "read_stats.csv"), 'w') as statfile:
                csv_writer = csv.writer(statfile)
                csv_writer.writerow(header)
                for i, lib in enumerate(sorted(stats_collector.keys(), key=natural_sort)):
                    csv_writer.writerow([i + 1] + stats_collector[lib])
            logging.debug("Preclean stats were successfully written to the file: OK")

            if seqformatter.debug:
                pdb.set_trace()

        # PERFORMS FASTQC TESTS
        if fastqc and gzoutput and not seqformatter.dryrun:
            fastqc_test(output.path, seqformatter.threads, extension="*.gz")
        elif fastqc and not seqformatter.dryrun:
            fastqc_test(output.path, seqformatter.threads)
        else:
            pass


def validate_reads(inpathR1, inpathR2, outdir, libname, pattern, datout, rename, gzout, skipcheck):
    """Perform all FASTQ checkups and save the output"""

    if seqformatter.dryrun:
        def write_seq_name(ffile, lline="", nname="", ssuffix=0, ccount=0):
            pass

    elif rename:
        datpath = os.path.join(datout.path, libname + seqformatter.suffix + ".dat")
        datfile = open(datpath, 'w')

        def write_seq_name(ffile, lline="", nname="", ssuffix=0, ccount=0):
            ffile.write("@{0:s}_ID{1:d}/{2:d}\n".format(nname, ccount, ssuffix))
            if ssuffix == 1:
                datfile.write("{0:s}_ID{1:d}\t{2:s}\n".format(nname, ccount, lline.split()[0]))
    else:
        # Switch the rename function off.
        def write_seq_name(ffile, lline="", nname="", ssuffix=0, ccount=0):
            ffile.write("@{0:s}\n".format(lline))

    if inpathR1 and inpathR2:
        logging.info("Running checkup analysis on paired reads in {0:s} library".format(libname.split()[0]))
        infileR1 = SeqIO(inpathR1, fileformat="FASTQ")
        infileR2 = SeqIO(inpathR2, fileformat="FASTQ")

        # CREATE AND OPEN OUTPUT FILES
        if seqformatter.dryrun:
            pass
        elif gzout:
            outpathR1 = os.path.join(outdir, libname + seqformatter.suffix + "_R1.fq.gz")
            outpathR2 = os.path.join(outdir, libname + seqformatter.suffix + "_R2.fq.gz")
            outfileR1 = gzip.open(outpathR1, 'wt')
            outfileR2 = gzip.open(outpathR2, 'wt')
        else:
            outpathR1 = os.path.join(outdir, libname + seqformatter.suffix + "_R1.fq")
            outpathR2 = os.path.join(outdir, libname + seqformatter.suffix + "_R2.fq")
            outfileR1 = open(outpathR1, 'w')
            outfileR2 = open(outpathR2, 'w')

    elif inpathR1 and not inpathR2:
        logging.info("Running checkup analysis on unpaired reads in '{0:s}' library".format(libname.split()[0]))
        infileR1 = SeqIO(inpathR1, fileformat="FASTQ")

        # CREATE AND OPEN OUTPUT FILES
        if seqformatter.dryrun:
            pass
        elif gzout:
            outpathR1 = os.path.join(outdir, libname + seqformatter.suffix + ".fq.gz")
            outfileR1 = gzip.open(outpathR1, 'wt')
        else:
            outpathR1 = os.path.join(outdir, libname + seqformatter.suffix + ".fq")
            outfileR1 = open(outpathR1, 'w')

    else:
        logging.error("No input files provided for FASTQ validation")
        raise EXONtoolsError("FASTQ file validation error")

    if inpathR2:
        count = 0
        filtered = 0
        for read1, read2 in zip(infileR1.read(), infileR2.read()):
            if skipcheck or (not read1.filtered and not read1.filtered and read1.name == read2.name and read1.pair == 1 and read2.pair == 2 and re.search(pattern, read1.identifier)):
                count += 1
                write_seq_name(ffile=outfileR1, lline=read1.identifier, nname=libname, ssuffix=1, ccount=count)
                outfileR1.write("{0:s}\n+{1:s}\n{2:s}\n".format(read1.seq, read1.info, read1.qual))
                write_seq_name(ffile=outfileR2, lline=read2.identifier, nname=libname, ssuffix=2, ccount=count)
                outfileR2.write("{0:s}\n+{1:s}\n{2:s}\n".format(read2.seq, read2.info, read2.qual))
            else:
                filtered += 1
        if not seqformatter.dryrun:
            [x.close() for x in [outfileR1, outfileR2]]

    else:
        count = 0
        filtered = 0
        for read1 in infileR1.read():
            if skipcheck or re.search(pattern, read1.identifier):
                count += 1
                seqname = read1.identifier.replace("/2", "/1").replace("_R2", "_R1")
                write_seq_name(ffile=outfileR1, lline=seqname, nname=libname, ssuffix=1, ccount=count)
                outfileR1.write("{0:s}\n+{1:s}\n{2:s}\n".format(read1.seq, read1.info, read1.qual))
            else:
                filtered += 1
        if not seqformatter.dryrun:
            outfileR1.close()

    if rename and not seqformatter.dryrun:
        datfile.close()

    if not seqformatter.dryrun:
        barcode = read1.barcode
    else:
        barcode = "NA"

    if not barcode:
        barcode = "NA"

    logging.info("Total {0:d} FASTQ reads passed EXONtools checkup in '{1:s}' library".format(count, libname))

    return {libname: [libname, barcode, count, filtered]}
