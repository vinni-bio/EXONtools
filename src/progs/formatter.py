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

        debug()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        readformatter.run_dry(SeqIO, getinput, output, makenewdir, executor)
        readformatter.set_debug(executor)

        if rqc:
            executor.setconfig("fastqc", "multiqc")
        debug()

        # GET INPUT PARAMETERS
        patterncheck = format_checks(skipcheck, fq2fa, customgrep, pattern, rename, gzoutput)
        debug()

        # IO SETTINGS
        getinput.format(['.fq', '.fastq', '.fastq.gz', '.fq.gz'])
        paired, unpaired = parseinput(inpath, forward, reverse)
        output(outdir)
        # tmpdir = makenewdir(name="tmp", fullname="temporary")
        logging.debug("IO settings: OK")
        debug()

        # COLLECT AND RUN ALL TASKS
        logging.info("Running the FASTQ read checkup analysis for all files")
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

        results = runtask(TASKS, "FASTQ file formatter")
        logging.debug("Checkup analysis succesfully finished: OK")
        debug()

        # SAVES STATS TO *.CSV TABLE
        savestats(results, output.path)

        # READ QUALITY TESTS
        if rqc and fq2fa:
            logging.warning("Quality analysis cannot be done on FASTA files. Please provide them in FASTQ format")
        else:
            runqc(rqc, output.path, message="Read formatting analysis", gzout=gzoutput)


def debug():
    """debuger"""
    if readformatter.debug:
        pdb.set_trace()


def format_checks(skipcheck, fq2fa, customgrep, pattern, rename, gzout):
    """PERFORM ALL CHECKS"""
    if rename:
        logging.warning("All reads in FASTQ files will be renamed according to library identifiers (--rename)")
    if gzout:
        logging.warning("All output FASTQ files will be compressed with gzip")
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
    logging.debug("Read formatter checks: OK")
    return patterncheck


def savestats(results, outpath):
    """Save stats for readformatter"""
    if readformatter.stats and not readformatter.dryrun:
        statdir = makenewdir(name=os.path.join(output.path, "STATS"), fullname="STATS")
        logging.info("Read stats will be saved to 'STATS/read_stats.csv'")
        header = ["No", "Library", "BARCODE", "#READS", "#FILTERED"]
        with open(os.path.join(statdir.path, "read_stats.csv"), 'w') as statfile:
            csv_writer = csv.writer(statfile)
            csv_writer.writerow(header)
            for i, lib in enumerate(sorted(results.keys(), key=natural_sort)):
                csv_writer.writerow([i + 1] + results[lib])
        logging.debug("Read stats were successfully written to the file: OK")
        debug()


def runqc(rqc, outpath, message, gzout):
    """Run QC analysis"""
    if rqc and not readformatter.stats and not readformatter.dryrun:
        statdir = makenewdir(name=os.path.join(outpath, "STATS"), fullname="STATS")
    else:
        statdir = os.path.join(outpath, "STATS")
    if rqc and gzout and not readformatter.dryrun:
        rqc_test(outpath, statdir.path, readformatter.threads, extension="*.gz", comment=message)
        debug()
    elif rqc and not readformatter.dryrun:
        rqc_test(outpath, statdir.path, readformatter.threads, comment=message)
        debug()
    else:
        pass


def parsejobs(jobs):
    """Convert jobs from list to dictionary"""
    job_collector = {}
    if jobs:
        for result in jobs:
            job_collector.update(result)
    else:
        logging.error("Multiprocessing produced no results")
        raise EXONtoolsError("Multiprocessing error")
    logging.debug("Checkup analysis succesfully finished: OK")
    return job_collector


def runtask(TASKS, message):
    """Run tasks"""
    if TASKS:
        processes_requested = set_threads(message, len(TASKS), readformatter.threads)
        pool = create_pool(processes_requested)
        jobs = hard_worker(run_instance, TASKS, pool)
        close_pool(pool)
        logging.debug("TASKs running process: OK")
        return parsejobs(jobs)


def gzopenfile(outpath, libname, gzout, suffix, fasta, mode):
    """Open file"""
    if fasta:
        ext = ".fa"
    else:
        ext = ".fq"
    if gzout:
        outpath = os.path.join(outpath, libname + suffix + ext + ".gz")
        return gzip.open(outpath, mode + 't')
    else:
        outpath = os.path.join(outpath, libname + suffix + ext)
        return gzip.open(outpath, mode)


def validate_reads(inpathR1, inpathR2, outdir, libname, pattern, datout, rename, gzout, skipcheck, fq2fa):
    """Perform all FASTQ checkups and save the output"""

    if inpathR2:
        suffix1 = '_R1' + readformatter.suffix
        suffix2 = '_R2' + readformatter.suffix
    else:
        suffix1 = readformatter.suffix

    if readformatter.dryrun:
        def write_seq_name(ffile, lline="", nname="", pprefix='@', ssuffix=0, ccount=0):
            pass
    elif rename:
        datpath = os.path.join(datout.path, libname + readformatter.suffix + ".dat")
        datfile = open(datpath, 'w')

        def write_seq_name(ffile, lline="", nname="", pprefix='@', ssuffix=0, ccount=0):
            ffile.write("{0:s}{1:s}_ID{2:d}/{3:d}\n".format(pprefix, nname, ccount, ssuffix))
            if ssuffix == 1:
                datfile.write("{0:s}_ID{1:d}\t{2:s}\n".format(nname, ccount, lline.split()[0]))
    else:
        # Switch off rename function
        def write_seq_name(ffile, lline="", nname="", pprefix='@', ssuffix=0, ccount=0):
            ffile.write("{0:s}{1:s}\n".format(pprefix, lline))

    if inpathR1 and inpathR2:
        logging.debug("Running checkup analysis on paired reads in {0:s} library".format(libname.split()[0]))
        infileR1 = SeqIO(inpathR1, fileformat="FASTQ")
        infileR2 = SeqIO(inpathR2, fileformat="FASTQ")

        # CREATE AND OPEN OUTPUT FILES
        if not readformatter.dryrun:
            outfileR1 = gzopenfile(outdir, libname, gzout=gzout, suffix=suffix1, fasta=fq2fa, mode="w")
            outfileR2 = gzopenfile(outdir, libname, gzout=gzout, suffix=suffix2, fasta=fq2fa, mode="w")
    elif inpathR1 and not inpathR2:
        logging.debug("Running checkup analysis on unpaired reads in '{0:s}' library".format(libname.split()[0]))
        infileR1 = SeqIO(inpathR1, fileformat="FASTQ")
        if not readformatter.dryrun:
            outfileR1 = gzopenfile(outdir, libname, gzout, suffix1, fasta=fq2fa, mode="w")
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
                    write_seq_name(ffile=outfileR1, lline=read1.identifier, nname=libname, pprefix='>', ssuffix=1, ccount=count)
                    outfileR1.write("{0:s}\n".format(read1.seq))
                    write_seq_name(ffile=outfileR2, lline=read2.identifier, nname=libname, pprefix='>', ssuffix=2, ccount=count)
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
