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
import itertools
import sys
import time

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, set_threads
from mains.EXT_errors import EXONtoolsError
from utils.revcomp import DNArevcomp
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from utils.fastqc import fastqc_test
from utils.progbar import printProgressBar


class demultiplexer(EXTprogram):
    """This program subdivides file with sequencing reads by library names or indexes"""

    name = "demultiplexer"

    def execute_program(self):
        args = self.args
        self.demultiplex_reads(args.inpath, args.outdir, args.barcode, args.inseqsearch, args.start, args.trim, args.pattern, args.tolerance, args.gzoutput, args.fastqc)

    def demultiplex_reads(self, inpath, outdir, barcode, inseqsearch, start, trim, pattern, tolerance, gzoutput, fastqc):

        if demultiplexer.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        demultiplexer.run_dry(SeqIO, getinput, output, makenewdir, executor)
        demultiplexer.set_debug(executor)

        if fastqc:
            executor.setconfig("fastqc")

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        demultiplexer.run_dry(getinput, output, makenewdir, executor)
        demultiplexer.set_debug(executor)

        if demultiplexer.debug:
            pdb.set_trace()

        if demultiplexer.extra:
            logging.warning("Extra string argument cannot be used in 'demultiplexer' program and therefore will be ommited")

        if barcode and pattern:
            logging.error("'--pattern' and '--barcode' flags are mutually exclusive. You cannot use them together")
            raise EXONtoolsError("Demultiplexer command line error")
        elif barcode:
            logging.info("Barcode option was enabled for demultiplexing")
            barcodes = get_barcodes(barcode)
            if inseqsearch:
                logging.warning("Library barcodes will be searched within read sequences")
        elif pattern:
            logging.info("Library name pattern option was enabled")
            if pattern == "default":
                pattern = "^@([^_:]+)[_:].*$"
            if pattern.count("(") != 1 or pattern.count(")") != 1:
                logging.error("Include round parentheses into your RegExp pattern to capture specific fragment")
                raise EXONtoolsError("GREP pattern error")
            logging.info("GREP pattern '{0}' will be used for demultiplexing".format(pattern))
        else:
            logging.error("You must provide either '--barcode' or '--pattern' option to start demultiplexer")
            raise EXONtoolsError("Demultiplexer command line error")

        # GET READ INPUT FILES
        logging.info("The following files with sequecing reads will be demultiplexed")
        getinput.format(['.fq', '.fastq', '.fastq.gz', '.fq.gz'])
        FileList = getinput(inpath).files

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE TMP DIRECTORIES
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        if gzoutput:
            logging.warning("All output FASTQ files will be compressed with gzip")

        logging.debug("Test and IO settings: OK")

        if demultiplexer.debug:
            pdb.set_trace()

        # Count total reads in FASTQ files
        TASKS = []
        for inpathfile in FileList:
            TASKS.append(worker(count_reads, [inpathfile]))

        logging.info("Counting reads in FASTQ files")
        if TASKS:
            processes_requested = set_threads("Read counter", len(TASKS), demultiplexer.threads)
            pool = create_pool(processes_requested)
            jobs = hard_worker(run_instance, TASKS, pool)
            close_pool(pool)

        totalreads = 0
        for result in jobs:
            for infile in result:
                totalreads += result[infile]

        logging.info("{0:d} reads are found in input file(s)".format(totalreads))

        ncluster = int(totalreads / demultiplexer.threads)
        ncluster_orig = int(totalreads / demultiplexer.threads) + 1

        counter = 0
        inFileList = []
        for inpathseq in FileList:
            infileseq = SeqIO(inpathseq, fileformat="FASTQ")
            outfilename = "clust1_" + os.path.basename(inpathseq)
            inFileList = [os.path.join(tmpdir.path, outfilename)]
            for inread in infileseq.read():
                counter += 1
                with open(os.path.join(tmpdir.path, outfilename), 'a') as outfile:
                    outfile.write("@{0:s} {1:s}\n{2:s}\n+{3:s}\n{4:s}\n".format(inread.name, inread.extra, inread.seq, inread.info, inread.qual))
                if counter > ncluster:
                    ncluster += ncluster_orig
                    outfilename = "clust" + str(ncluster) + "_" + os.path.basename(inpathseq)
                    inFileList.append(os.path.join(tmpdir.path, outfilename))

        # Collects all tasks
        TASKS = []

        logging.debug("Preparing tasks for your demultiplex analysis")

        outmpdirs = []
        for inpathfile in inFileList:
            outtmpdir = makenewdir(os.path.join(tmpdir.path, inpathfile.split("_")[0]))
            outmpdirs.append(outtmpdir.path)
            if barcode:
                TASKS.append(worker(debarcode, [barcodes, inpathfile, outtmpdir.path, "_demux", tolerance, inseqsearch, trim, start, gzoutput, ncluster_orig]))
            else:
                TASKS.append(worker(demux, [inpathfile, outtmpdir.path, pattern, "_unpaired", gzoutput, ncluster_orig]))

        if demultiplexer.debug:
            pdb.set_trace()

        logging.info("Running demultiplexing analysis for all selected files")
        if TASKS:
            processes_requested = set_threads("Read demultiplexer", len(TASKS), demultiplexer.threads)
            pool = create_pool(processes_requested)
            jobs = hard_worker(run_instance, TASKS, pool)
            close_pool(pool)

        mergedirs(outmpdirs, output.path)

        tmpdir.delete()

        stats_collector = {}
        final_stats = {}
        if jobs:
            for result in jobs:
                stats_collector.update(result)
            del jobs
        else:
            logging.error("Multiprocessing error in demultiplexing analysis")
            raise EXONtoolsError("Demultiplex error")

        for filename in stats_collector:
            for libname in stats_collector[filename]:
                try:
                    final_stats[libname] += stats_collector[filename][libname]
                except KeyError:
                    final_stats[libname] = stats_collector[filename][libname]

        logging.debug("Demultiplex analysis succesfully finished: OK")

        if demultiplexer.debug:
            pdb.set_trace()

        if demultiplexer.stats and not demultiplexer.dryrun:
            statdir = makenewdir(name=os.path.join(output.path, "STATS"), fullname="STATS")
            logging.info("Demultiplexing stats will be saved to 'STATS/demux_stats.csv'")
            header = ["No", "Library", "#READS"]
            with open(os.path.join(statdir.path, "demux_stats.csv"), 'w') as statfile:
                csv_writer = csv.writer(statfile)
                csv_writer.writerow(header)
                for i, libname in enumerate(sorted(final_stats.keys(), key=natural_sort)):
                    csv_writer.writerow([i + 1, libname, final_stats[libname]])
            logging.debug("Demultiplexing stats were successfully written to the file: OK")

            if demultiplexer.debug:
                pdb.set_trace()

        # PERFORMS FASTQC TESTS
        if not os.path.exists(statdir.path):
            statdir = makenewdir(name=os.path.join(output.path, "STATS"), fullname="STATS")
        if fastqc and gzoutput and not demultiplexer.dryrun:
            fastqc_test(output.path, statdir.path, demultiplexer.threads, extension="*.gz", comment="Demultiplexing analysis")
        elif fastqc and not demultiplexer.dryrun:
            fastqc_test(output.path, statdir.path, demultiplexer.threads, comment="Demultiplexing analysis")
        else:
            pass


def get_barcodes(barcpath):
    """This function parses the file with barcodes and returns the dictionary with barcodes as keys"""

    logging.debug("Reading library barcodes")
    barcpath = os.path.abspath(barcpath)

    try:
        with open(barcpath, 'r') as barcfile:
            logging.info("Parsing file with barcodes from:")
            logging.info(barcpath)
            barcodes = {}
            for line in barcfile:
                line = line.strip().upper()
                if line:
                    dat = re.split('\s|\t', line)
                    if len(dat) == 3:
                        dat[2] = DNArevcomp(dat[2])
                    elif len(dat) > 3:
                        logging.error("Your file with barcodes contains more than 3 columns!")
                        raise EXONtoolsError("Demultiplexer barcode file format error")
                    elif len(dat) == 1:
                        logging.error("Your file with barcodes contains less than 2 columns!")
                        raise EXONtoolsError("Demultiplexer barcode file format error")
                    barcodes["+".join(dat[1:])] = dat[0]
    except IOError:
        logging.error("Your file with library indexes does not exist")
        raise EXONtoolsError("Demultiplexer wrong input path for barcode file")

    if len(barcodes) == 0:
        logging.error("No barcodes were found in the barcode file")
        raise EXONtoolsError("Demultiplexer barcode file input error")
    logging.debug("Read library indexes: OK")
    return barcodes


def debarcode(barcodes, filepath, outpath, suffix, tolerance, inseq, trim, start, gzout, total):
    """Demultiplex reads by barcodes sequences"""

    total_defined = 0
    undefined_count = 0
    barc_stats = {x: 0 for x in list(barcodes.values())}
    if inseq:
        barcrange = (min(map(len, list(barcodes.keys()))), max(map(len, list(barcodes.keys()))))
        if len(barcrange) != len(set(barcrange)):
            barcrange = barcrange[0]
        infileseq = SeqIO(filepath, fileformat="FASTQ")
        for i, read in enumerate(infileseq.read()):
            printProgressBar(i, total, prefix="", suffix="", decimals=1, length=100, fill='█', printEnd="\r")
            barc_list = []
            selected_barc = {}
            for i in range(barcrange[0], barcrange[1] + 1):
                barc = read.seq[start - 1:start - 1 + i]
                result = tolerate_barcode(barc, tolerance)
                barc_list = barc_list + result

                barc = read.seq[len(read.seq) - start - 1 - i:len(read.seq) - start - 1]
                result = tolerate_barcode(barc, tolerance)
                barc_list = barc_list + result

            for barc in barc_list:
                if barc in barcodes:
                    selected_barc[barc] = barcodes[barc]            # checks corresponding barcodes

            # SAVES READS WITH A SINGLE MATCH
            if len(selected_barc) == 1:
                libname = list(selected_barc.values())[0]
                barcseq = list(selected_barc.keys())[0]
                barc_stats[libname] += 1
                total_defined += 1
                file_name = os.path.join(outpath, libname + suffix + ".fq")
                readseq = read.seq[len(barcseq) + trim:]
                qualseq = read.qual[len(barcseq) + trim:]
                if gzout:
                    with gzip.open(file_name + ".gz", "at") as outfile:
                        outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, readseq, read.info, qualseq))
                else:
                    with open(file_name, "a") as outfile:
                        outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, readseq, read.info, qualseq))

            # SAVES UNMATCHED READS or READS WITH MULTIPLE MATCHES
            else:
                file_name = os.path.join(outpath, "Undefined_reads" + suffix + ".fq")
                undefined_count += 1
                if gzout:
                    with gzip.open(file_name + ".gz", "at") as outfile:
                        outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, read.seq, read.info, read.qual))
                else:
                    with open(file_name, "a") as outfile:
                        outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, read.seq, read.info, read.qual))

    else:
        infileseq = SeqIO(filepath, fileformat="FASTQ")

        for read in infileseq.read():
            if read.barcode:
                selected_barc = {}                                      # initiates the dict of selected barcodes
                barc_list = tolerate_barcode(read.barcode, tolerance)        # makes the list with barcode variations

                for barc in barc_list:
                    if barc in barcodes:
                        selected_barc[barc] = barcodes[barc]            # checks corresponding barcodes

                # SAVES READS WITH A SINGLE MATCH
                if len(selected_barc) == 1:
                    libname = list(selected_barc.values())[0]
                    barc_stats[libname] += 1
                    total_defined += 1
                    file_name = os.path.join(outpath, libname + suffix + ".fq")
                    if gzout:
                        with gzip.open(file_name + ".gz", "at") as outfile:
                            outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, read.seq, read.info, read.qual))
                    else:
                        with open(file_name, "a") as outfile:
                            outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, read.seq, read.info, read.qual))

                # SAVES UNMATCHED READS or READS WITH MULTIPLE MATCHES
                else:
                    file_name = os.path.join(outpath, "Undefined_reads" + suffix + ".fq")
                    undefined_count += 1
                    if gzout:
                        with gzip.open(file_name + ".gz", "at") as outfile:
                            outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, read.seq, read.info, read.qual))
                    else:
                        with open(file_name, "a") as outfile:
                            outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, read.seq, read.info, read.qual))
            else:
                file_name = os.path.join(outpath, "Empty_barcodes" + suffix + ".fq")
                undefined_count += 1
                if gzout:
                    with gzip.open(file_name + ".gz", "at") as outfile:
                        outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, read.seq, read.info, read.qual))
                else:
                    with open(file_name, "a") as outfile:
                        outfile.write("@{0} {1}\n{2}\n+{3}\n{4}\n".format(read.name, read.extra, read.seq, read.info, read.qual))

    barc_stats["undefined_count"] = undefined_count
    os.remove(filepath)
    logging.info("Total {0} reads were demultiplexed in {1}".format(total_defined, os.path.basename(filepath)))
    logging.info("Total {0} reads were not identified in {1}".format(undefined_count, os.path.basename(filepath)))

    return {os.path.basename(filepath): barc_stats}


def tolerate_barcode(barcode, tolerance):
    """
    Takes barcode sequence (single or dual) and returns the list
    with all possible combindations allowing some tolerance for mismatches (max=2).
    """

    barcode_list = []                                   # initiates the list for barcodes
    dat = barcode.split("+")                            # split if dual barcode

    for seq in dat:                                     # for each barcode in the list do:
        bucket1 = [seq]                                 # bucket 1 includes barcodes for modification
        bucket2 = []                                    # empty bucket 2
        bucket3 = [seq]                                 # bucket 3 includes the original barcodes
        step = 0                                        # initiates step counter
        while(step < tolerance):                        # does operation for each tolerance level
            for x in bucket1:                           #
                for i, base in enumerate(list(x)):      # splits barcode by nucleotides
                    nuc_list = ["A", "T", "C", "G"]     # list of bases for replacement in barcode seq
                    try:
                        nuc_list.remove(base)           # removes current nucleotide
                    except ValueError:                  # returns replacements for uknown base
                        result = replaceN(barcode)
                        return result
                    for nuc in nuc_list:                # for each nucleotide in replacement list do:
                        seq_list = list(x)              # split original barcode
                        seq_list[i] = nuc               # replace each base with one nucleotide
                        final = "".join(seq_list)       # join all bases
                        bucket2.append(final)           # saves new barcode in bucket2
            bucket1 = []                                # empties bucket1
            for y in bucket2:                           # for each barcode in bucket 2 do:
                bucket1.append(y)                       # put barcode into bucket 1 to repeat all thing on them
                if y not in bucket3:                    # append all new barcodes to bucket 3
                    bucket3.append(y)
            bucket2 = []                                # empties bucket 2
            step += 1
        barcode_list.append(bucket3)                     # adds all barcodes to list
    final_list = list(itertools.product(*barcode_list))  # make all combination pairs
    result = ["+".join(x) for x in final_list]           # produce result with joint barcodes
    return result


def replaceN(barcode):
    dat = barcode.split("+")                            # split barcode if dual
    barcode_list = []                                   # create empty list for new barcodes
    nuc_list = ["A", "T", "C", "G"]                     # normal bases
    for seq in dat:                                     # check each barcode
        new_seqs = []                                   # empty list for barcode variations
        for nuc in list(seq):                           # check strange
            if nuc not in nuc_list:
                for k in nuc_list:
                    newseq = re.sub(nuc, k, seq)
                    if newseq not in new_seqs:
                        new_seqs.append(newseq)
        if new_seqs:                                     # append new seqs if found
            barcode_list.append(new_seqs)
        else:
            barcode_list.append([seq])                   # append old seq if OK
    final_list = list(itertools.product(*barcode_list))  # make all combination pairs
    result = ["+".join(x) for x in final_list]           # produce result with joint barcodes
    return result


def demux(filepath, outpath, pattern, suffix, gzout, total):
    """Demultiplex reads using GREP pattern in read name identifier"""

    total_defined = 0
    undefined_count = 0
    demux_stats = {}

    infileseq = SeqIO(filepath, fileformat="FASTQ")

    for read in infileseq.read():
        libname = check_pattern(pattern, "@" + read.name + " " + read.extra)

        if libname:
            total_defined += 1
            try:
                demux_stats[libname] += 1
            except KeyError:
                demux_stats[libname] = 1

            if gzout:
                with gzip.open(os.path.join(outpath, libname + suffix + ".fq.gz"), "at") as outfile:
                    outfile.write("@{0}\n{1}\n+{2}\n{3}\n".format(read.name + " " + read.extra.strip(), read.seq, read.info, read.qual))
            else:
                with open(os.path.join(outpath, libname + suffix + ".fq"), "a") as outfile:
                    outfile.write("@{0}\n{1}\n+{2}\n{3}\n".format(read.name + " " + read.extra.strip(), read.seq, read.info, read.qual))
        else:
            undefined_count += 1
            if gzout:
                with gzip.open(os.path.join(outpath, "Undefined_reads" + suffix + ".fq.gz"), "at") as outfile:
                    outfile.write("@{0}\n{1}\n+{2}\n{3}\n".format(read.name + " " + read.extra.strip(), read.seq, read.info, read.qual))
            else:
                with open(os.path.join(outpath, "Undefined_reads" + suffix + ".fq"), "a") as outfile:
                    outfile.write("@{0}\n{1}\n+{2}\n{3}\n".format(read.name + " " + read.extra.strip(), read.seq, read.info, read.qual))

    demux_stats["undefined_count"] = undefined_count

    logging.info("Total {0} reads were demultiplexed in {1}".format(total_defined, os.path.basename(filepath)))
    logging.info("Total {0} reads were not identified in {1}".format(undefined_count, os.path.basename(filepath)))
    return demux_stats


def check_pattern(pattern, line):
    """Check and return pattern result"""

    maskcheck = re.match(pattern, line)
    if maskcheck:
        result = maskcheck.groups()
        if len(result) > 1:
            logging.error("Number of matches is more than one")
            logging.error("Please check your pattern settings")
            raise EXONtoolsError("Search pattern error")
    else:
        return
        # logging.error("No matches found based on provided pattern. Please check your pattern settings")
        # raise EXONtoolsError("Search pattern error")

    return result[0]


def count_reads(inpath):

    counter = 0
    infileseq = SeqIO(inpath, fileformat="FASTQ")
    for r in infileseq.read():
        counter += 1

    return {os.path.basename(inpath): counter}


def mergedirs(DirList, outpath):
    """Merge files from all folders"""

    for indir in DirList:
        FileList = os.listdir(indir)
        for file in FileList:
            with open(os.path.join(indir, file), 'r') as infile:
                with open(os.path.join(outpath, file), 'a') as outfile:
                    for line in infile:
                        outfile.write(line)
            os.remove(os.path.join(indir, file))
