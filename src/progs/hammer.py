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


class hammer(EXTprogram):
    """This program finds and fixes sequencing errors within fastq files using SPAdes"""

    name = "hammer"

    def execute_program(self):
        args = self.args
        self.correct_errors(args.inpath, args.forward, args.reverse, args.unpaired, args.outdir, args.gzoutput, args.torrent, args.rqc)

    def correct_errors(self, inpath, forward, reverse, unpaired_in, outdir, gzoutput, torrent, rqc):

        debug()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        hammer.run_dry(SeqIO, getinput, output, makenewdir, executor)
        hammer.set_debug(executor)
        if rqc:
            executor.setconfig("fastqc", "multiqc", "spades")
        else:
            executor.setconfig("spades")
        debug()

        # IO SETTINGS
        paired, unpaired = fixunpaired(inpath, forward, reverse, unpaired_in)
        librs = findlibrs(paired, unpaired)
        output(outdir)
        tmpdir = makenewdir(name="tmp", fullname="temporary")
        log_dir = makenewdir(name="spades_logs", fullname="LOG")
        logging.info("All SPAdes log files will be saved in '{0:s}' directory".format(os.path.basename(log_dir.path)))
        logging.debug("IO settings: OK")
        debug()

        # SET SPADES PARAMETERS
        SPADES_params = spadespars(torrent, gzoutput)
        debug()

        # RUN SPADES ANALYSIS
        logging.debug("Start SPAdes analysis")
        libtmps = {}
        for lib in librs:
            infiles = makeinfiles(lib, paired, unpaired)
            libouttmp = makelibtemp(lib, tmpdir.path)
            libtmps[lib] = libouttmp

            run_executor(executor(
                program="SPAdes",
                params=[SPADES_params, libouttmp, hammer.threads, int(hammer.memory), infiles],
                conditions={"positive": [hammer.threads], "pathexists": [libouttmp]},
                custom_arg_string=hammer.extra + " > /dev/null"
            ))

            copylog(lib, log_dir.path, libouttmp)
            logging.info("SPAdes read correction for the library '{0:s}' is done".format(lib))

        logging.debug("Successfully finished all corrections: OK")
        debug()

        # CREATE AND RUN TASKS
        TASKS = []
        logging.debug("Preparing files for multiproc analaysis")
        for lib in sorted(libtmps.keys()):
            corrdir = os.path.join(libtmps[lib], "corrected")
            FileList = [os.path.join(corrdir, x) for x in os.listdir(corrdir) if x.endswith(".fastq")]
            TASKS.append(worker(process_stats, [FileList, output.path, lib, hammer.suffix, gzoutput]))
        logging.info("Processing files with corrected sequencing data")
        results = runtask(TASKS, "SPAdes stat collector")
        logging.debug("Read correction analysis succesfully finished: OK")
        debug()

        # SAVE STATS
        savestats(results, output.path)

        # READ QUALITY TESTS
        runqc(rqc, outpath=output.path, message="Corrected reads", gzout=gzoutput)

        if not hammer.keeptmp:
            tmpdir.delete()


def debug():
    """debuger"""
    if hammer.debug:
        pdb.set_trace()


def savestats(results, outpath):
    """Save all stats"""
    if hammer.stats and not hammer.dryrun:
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


def runqc(rqc, outpath, message, gzout):
    """Run QC analysis"""
    if rqc and not hammer.stats and not hammer.dryrun:
        statdir = makenewdir(name=os.path.join(outpath, "STATS"), fullname="STATS").path
    else:
        statdir = os.path.join(outpath, "STATS")
    if rqc and gzout and not hammer.dryrun:
        rqc_test(outpath, statdir, hammer.threads, extension="*.gz", comment=message)
        debug()
    elif rqc and not hammer.dryrun:
        rqc_test(outpath, statdir, hammer.threads, comment=message)
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


def spadespars(torrent, gzoutput):
    """set spades parameters"""
    logging.debug("Setting up SPAdes parameters")
    SPADES_params = '--only-error-correction'
    if gzoutput:
        logging.warning("All output FASTQ files will be compressed with gzip")
    if hammer.extra:
        logging.warning("The following arguments will be added to 'SPAdes' command line:")
        logging.warning(hammer.extra)
    if torrent:
        SPADES_params = SPADES_params + " --iontorrent"
        logging.warning("Ion Torrent data type is set for SPAdes")
        logging.info("Correcting sequencing errors with SPAdes IonHammer")
    else:
        logging.info("Correcting sequencing errors with SPAdes BayesHammer")
    logging.debug("SPAdes parameters: OK")
    return SPADES_params


def makeinfiles(lib, paired, unpaired):
    """Create spades input files"""
    if lib + "_paired" in paired and lib + "_unpaired" in unpaired:
        infiles = "-1 {0:s} -2 {1:s} -s {2:s}".format(paired[lib + "_paired"][0], paired[lib + "_paired"][1], unpaired[lib + "_unpaired"][0])
    elif lib + "_paired" in paired and not lib + "_unpaired" in unpaired:
        infiles = "-1 {0:s} -2 {1:s}".format(paired[lib + "_paired"][0], paired[lib + "_paired"][1])
    elif not lib + "_paired" in paired and lib + "_unpaired" in unpaired:
        infiles = "-s {0:s}".format(unpaired[lib + "_unpaired"][0])
    else:
        logging.error("You found the bug during parsing of SPAdes input files")
        raise EXONtoolsError("Input files command line error")
    logging.debug("SPAdes input files: OK")
    return infiles


def parsejobs(jobs):
    """Convert jobs from list to dictionary"""
    job_collector = {}
    if jobs:
        for result in jobs:
            job_collector.update(result)
    else:
        logging.error("Multiprocessing produced no results")
        raise EXONtoolsError("Multiprocessing error")
    return job_collector


def runtask(TASKS, message):
    """Run tasks"""
    if TASKS:
        processes_requested = set_threads(message, len(TASKS), hammer.threads)
        pool = create_pool(processes_requested)
        jobs = hard_worker(run_instance, TASKS, pool)
        close_pool(pool)
        logging.debug("TASKs running process: OK")
        return parsejobs(jobs)


def copylog(sample, logdir, libouttmp):
    """Copy SPAdes output"""
    logpath = os.path.join(logdir, sample + "_spades.log")
    if not hammer.dryrun:
        spadelogpath = os.path.join(libouttmp, "spades.log")
    try:
        shutil.move(spadelogpath, logpath)
        logging.debug("Moved SPAdes logs for '{0:s}' library: OK".format(sample))
    except IOError:
        logging.error("EXONtools cannot find 'spades.log' file for '{0:s}' sample".format(sample))
        raise EXONtoolsError


def check_spades_output(infile):
    """Find corrected files"""
    check = os.path.basename(infile).split('.')[0]
    if "_R1" in check and 'unpaired' not in check:
        return "forward"
    elif "_R2" in check and 'unpaired' not in check:
        return "reverse"
    else:
        return "unpaired"


def make_outpath(sample, direction, infile, outdir, suffix):
    """Make outpath for spades output file"""
    if direction == 'forward' or direction == 'reverse':
        outfile_name = os.path.basename(infile).split(".")[0] + suffix + ".fq"
        outfilepath = os.path.join(outdir, outfile_name)
    else:
        outfile_name = os.path.basename(infile).split(".")[0] + suffix + ".fq"
        outfilepath = os.path.join(outdir, sample + "_unpaired" + suffix + ".fq")
    logging.debug("Made output path for '{0:s}' library: OK".format(sample))
    return outfilepath


def opengzfile(inpath, gzout=False, mode='a'):
    if gzout:
        return gzip.open(inpath + ".gz", mode + 't')
    else:
        return open(inpath, mode)


def process_stats(inpathlist, outpath, sample, suffix, gzoutput):
    """Collecting stats from read correction results"""

    stats = {sample: {
        "forward": {"outpath": None, "trimmed": 0, "changed": 0, "counts": 0},
        "reverse": {"outpath": None, "trimmed": 0, "changed": 0, "counts": 0},
        "unpaired": {"outpath": None, "trimmed": 0, "changed": 0, "counts": 0}
    }}
    logging.info("Cleaning SPAdes output for '{0:s}' library".format(sample))
    for inpath in inpathlist:
        direction = check_spades_output(inpath)
        outfilepath = make_outpath(sample, direction, inpath, outpath, suffix)
        stats[sample][direction]['outpath'] = outfilepath

        if direction == 'unpaired':
            def writefile(outfile, name, seq, qual):
                outfile.write("@{0:s}\n{1:s}\n+\n{2:s}\n".format(name.replace("_paired", "_unpaired"), seq, qual))
        else:
            def writefile(outfile, name, seq, qual):
                outfile.write("@{0:s}\n{1:s}\n+\n{2:s}\n".format(name, seq, qual))
        outfile = opengzfile(outfilepath, gzout=gzoutput)
        infile = SeqIO(inpath, fileformat="FASTQ")
        for read in infile.read():
            if "ltrim=" in read.info:
                stats[sample][direction]["trimmed"] += 1
            if "changed" in read.info:
                changed_data = int(read.info.split("changed:")[-1].split(" ")[0])
                stats[sample][direction]["changed"] = stats[sample][direction]["changed"] + changed_data
            writefile(outfile, read.name, read.seq, read.qual)
        stats[sample][direction]['counts'] = stats[sample][direction]['counts'] + infile.total
        outfile.close()
        logging.debug("SPAdes output for '{0:s}' library: OK".format(sample))
    return stats
