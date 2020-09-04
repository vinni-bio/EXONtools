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

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, parseinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, set_threads
from mains.EXT_errors import EXONtoolsError
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from utils.rqc import rqc_test


class readformatter(EXTprogram):
    """This program verifies names and structure of raw reads"""

    name = "readformatter"

    def execute_program(self):
        args = self.args
        self.format_seqs(args.inpath, args.forward, args.reverse, args.outdir, args.gzoutput, args.fq2fa, args.rename, args.pattern, args.customgrep, args.skipcheck, args.rqc)

    def format_seqs(self, inpath, forward, reverse, outdir, gzoutput, fq2fa, rename, pattern, customgrep, skipcheck, rqc):

        if readformatter.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        readformatter.run_dry(SeqIO, getinput, output, makenewdir, executor)
        readformatter.set_debug(executor)

        if rqc:
            executor.setconfig("fastqc", "multiqc")

        if readformatter.debug:
            pdb.set_trace()

        if readformatter.extra:
            logging.warning("Extra string argument cannot be used in 'readformatter' program and therefore will be ommited")

        if skipcheck:
            logging.warning("Skipping FASTQ read identifier checkup... Only parsing output will be performed")

        if fq2fa:
            logging.warning("All reads will be transformed to FASTA format (--fasta)")

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
            elif pattern.lower() == "torrent":
                logging.info("Reads will be verified for conforming the Ion Torrent raw read format")
                patterncheck = "^\w+:\d+:\d+$"
            elif pattern.lower() == "exontools":
                logging.info("Reads will be verified for conforming the EXONtools read name format")
                patterncheck = "^([^\s|_]+_[un]*paired_ID\d+)/\d$"
            else:
                patterncheck = customgrep

        # GET READ INPUT FILES
        getinput.format(['.fq', '.fastq', '.fastq.gz', '.fq.gz'])
        paired, unpaired = parseinput(inpath, forward, reverse)

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

        if readformatter.debug:
            pdb.set_trace()

        # Collects all tasks
        TASKS = []
        if rename:
            logging.warning("Old and new read names will be saved to 'CATALOG' directory")
            datnames = makenewdir(name="CATALOG", fullname="CATALOG")
        else:
            datnames = None
        if paired:
            for lib in sorted(paired.keys(), key=natural_sort):
                TASKS.append(worker(validate_reads, [paired[lib][0], paired[lib][1], output.path, lib, patterncheck, datnames, rename, gzoutput, skipcheck, fq2fa]))
        if unpaired:
            for lib in sorted(unpaired.keys(), key=natural_sort):
                TASKS.append(worker(validate_reads, [unpaired[lib][0], None, output.path, lib, patterncheck, datnames, rename, gzoutput, skipcheck, fq2fa]))

        logging.info("Running the FASTQ read checkup analysis for all files")
        if TASKS:
            processes_requested = set_threads("FASTQ file formatter", len(TASKS), readformatter.threads)
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

        if readformatter.debug:
            pdb.set_trace()

        # SAVES STATS TO *.CSV TABLE
        if readformatter.stats and not readformatter.dryrun:
            statdir = makenewdir(name=os.path.join(output.path, "STATS"), fullname="STATS")
            logging.info("Read stats will be saved to 'STATS/read_stats.csv'")
            header = ["No", "Library", "BARCODE", "#READS", "#FILTERED"]
            with open(os.path.join(statdir.path, "read_stats.csv"), 'w') as statfile:
                csv_writer = csv.writer(statfile)
                csv_writer.writerow(header)
                for i, lib in enumerate(sorted(stats_collector.keys(), key=natural_sort)):
                    csv_writer.writerow([i + 1] + stats_collector[lib])
            logging.debug("Read stats were successfully written to the file: OK")

            if readformatter.debug:
                pdb.set_trace()

        # PERFORMS READ QUALITY TESTS
        if rqc and fq2fa:
            logging.warning("Quality analysis cannot be done on FASTA files. Please provide them in FASTQ format")
        else:
            if rqc and not readformatter.stats:
                statdir = makenewdir(name=os.path.join(output.path, "STATS"), fullname="STATS")
            if rqc and gzoutput and not readformatter.dryrun:
                rqc_test(output.path, statdir.path, readformatter.threads, extension="*.gz", comment="Read formatting analysis")
            elif rqc and not readformatter.dryrun:
                rqc_test(output.path, statdir.path, readformatter.threads, comment="Read formatting analysis")
            else:
                pass


def validate_reads(inpathR1, inpathR2, outdir, libname, pattern, datout, rename, gzout, skipcheck, fq2fa):
    """Perform all FASTQ checkups and save the output"""

    if readformatter.dryrun:
        def write_seq_name(ffile, lline="", nname="", pprefix ='@', ssuffix=0, ccount=0):
            pass

    elif rename:
        datpath = os.path.join(datout.path, libname + readformatter.suffix + ".dat")
        datfile = open(datpath, 'w')

        def write_seq_name(ffile, lline="", nname="", pprefix ='@', ssuffix=0, ccount=0):
            ffile.write("{0:s}{1:s}_ID{2:d}/{3:d}\n".format(pprefix,nname, ccount, ssuffix))
            if ssuffix == 1:
                datfile.write("{0:s}_ID{1:d}\t{2:s}\n".format(nname, ccount, lline.split()[0]))
    else:
        # Switch the rename function off.
        def write_seq_name(ffile, lline="", nname="", pprefix ='@', ssuffix=0, ccount=0):
            ffile.write("{0:s}{1:s}\n".format(pprefix,lline))

    if inpathR1 and inpathR2:
        logging.info("Running checkup analysis on paired reads in {0:s} library".format(libname.split()[0]))
        infileR1 = SeqIO(inpathR1, fileformat="FASTQ")
        infileR2 = SeqIO(inpathR2, fileformat="FASTQ")

        # CREATE AND OPEN OUTPUT FILES
        if readformatter.dryrun:
            pass
        elif gzout and fq2fa:
            outpathR1 = os.path.join(outdir, libname + readformatter.suffix + "_R1.fa.gz")
            outpathR2 = os.path.join(outdir, libname + readformatter.suffix + "_R2.fa.gz")
            outfileR1 = gzip.open(outpathR1, 'wt')
            outfileR2 = gzip.open(outpathR2, 'wt')
        elif gzout:
            outpathR1 = os.path.join(outdir, libname + readformatter.suffix + "_R1.fq.gz")
            outpathR2 = os.path.join(outdir, libname + readformatter.suffix + "_R2.fq.gz")
            outfileR1 = gzip.open(outpathR1, 'wt')
            outfileR2 = gzip.open(outpathR2, 'wt')
        elif fq2fa:
            outpathR1 = os.path.join(outdir, libname + readformatter.suffix + "_R1.fa")
            outpathR2 = os.path.join(outdir, libname + readformatter.suffix + "_R2.fa")
            outfileR1 = open(outpathR1, 'w')
            outfileR2 = open(outpathR2, 'w')
        else:
            outpathR1 = os.path.join(outdir, libname + readformatter.suffix + "_R1.fq")
            outpathR2 = os.path.join(outdir, libname + readformatter.suffix + "_R2.fq")
            outfileR1 = open(outpathR1, 'w')
            outfileR2 = open(outpathR2, 'w')

    elif inpathR1 and not inpathR2:
        logging.info("Running checkup analysis on unpaired reads in '{0:s}' library".format(libname.split()[0]))
        infileR1 = SeqIO(inpathR1, fileformat="FASTQ")

        # CREATE AND OPEN OUTPUT FILES
        if readformatter.dryrun:
            pass
        elif gzout and fq2fa:
            outpathR1 = os.path.join(outdir, libname + readformatter.suffix + ".fa.gz")
            outfileR1 = gzip.open(outpathR1, 'wt')
        elif gzout:
            outpathR1 = os.path.join(outdir, libname + readformatter.suffix + ".fq.gz")
            outfileR1 = gzip.open(outpathR1, 'wt')
        elif fq2fa:
            outpathR1 = os.path.join(outdir, libname + readformatter.suffix + ".fa")
            outfileR1 = open(outpathR1, 'w')
        else:
            outpathR1 = os.path.join(outdir, libname + readformatter.suffix + ".fq")
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
                if fq2fa:
                    write_seq_name(ffile=outfileR1, lline=read1.identifier, nname=libname, pprefix = '>', ssuffix=1, ccount=count)
                    outfileR1.write("{0:s}\n".format(read1.seq))
                    write_seq_name(ffile=outfileR2, lline=read2.identifier, nname=libname, pprefix = '>', ssuffix=2, ccount=count)
                    outfileR2.write("{0:s}\n".format(read2.seq))
                else:
                    write_seq_name(ffile=outfileR1, lline=read1.identifier, nname=libname, ssuffix=1, ccount=count)
                    outfileR1.write("{0:s}\n+{1:s}\n{2:s}\n".format(read1.seq, read1.info, read1.qual))
                    write_seq_name(ffile=outfileR2, lline=read2.identifier, nname=libname, ssuffix=2, ccount=count)
                    outfileR2.write("{0:s}\n+{1:s}\n{2:s}\n".format(read2.seq, read2.info, read2.qual))
            else:
                filtered += 1
        if not readformatter.dryrun:
            [x.close() for x in [outfileR1, outfileR2]]

    else:
        count = 0
        filtered = 0
        for read1 in infileR1.read():
            if skipcheck or re.search(pattern, read1.identifier):
                count += 1
                if fq2fa:
                    seqname = read1.identifier.replace("/2", "/1").replace("_R2", "_R1")
                    write_seq_name(ffile=outfileR1, lline=seqname, nname=libname, pprefix=">", ssuffix=1, ccount=count)
                    outfileR1.write("{0:s}\n".format(read1.seq))
                else:
                    seqname = read1.identifier.replace("/2", "/1").replace("_R2", "_R1")
                    write_seq_name(ffile=outfileR1, lline=seqname, nname=libname, ssuffix=1, ccount=count)
                    outfileR1.write("{0:s}\n+{1:s}\n{2:s}\n".format(read1.seq, read1.info, read1.qual))
            else:
                filtered += 1
        if not readformatter.dryrun:
            outfileR1.close()

    if rename and not readformatter.dryrun:
        datfile.close()

    if not readformatter.dryrun:
        barcode = read1.barcode
    else:
        barcode = "NA"

    if not barcode:
        barcode = "NA"

    if fq2fa:
        logging.info("Total {0:d} FASTA sequences passed EXONtools checkup in '{1:s}' library".format(count, libname))
    else:
        logging.info("Total {0:d} FASTQ sequences passed EXONtools checkup in '{1:s}' library".format(count, libname))

    return {libname: [libname, barcode, count, filtered]}
