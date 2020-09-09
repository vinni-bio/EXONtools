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
import itertools
import time
import sys

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, parseinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_errors import EXONtoolsError
from utils.revcomp import DNArevcomp
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from utils.rqc import rqc_test
from utils.progbar import printProgressBar


class demultiplexer(EXTprogram):
    """This program subdivides file with sequencing reads by library names or indexes"""

    name = "demultiplexer"

    def execute_program(self):
        args = self.args
        self.demultiplex_reads(args.forward, args.reverse, args.outdir, args.barcode, args.inseqsearch, args.start, args.trim, args.pattern, args.tolerance, args.gzoutput, args.rqc)

    def demultiplex_reads(self, forward, reverse, outdir, barcode, inseqsearch, start, trim, pattern, tolerance, gzoutput, rqc):

        debug()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        logging.debug("Setting dry and debugging modes for subclasses")
        demultiplexer.run_dry(SeqIO, getinput, output, makenewdir, executor)
        demultiplexer.set_debug(executor)

        if rqc:
            executor.setconfig("fastqc", "multiqc")
        debug()

        # GET INPUT PARAMETERS
        barcodes, pattern = demux_checks(barcode, pattern, inseqsearch, gzoutput, start)
        debug()

        # IO SETTINGS
        FileList = demux_io(forward, reverse)
        output(outdir)
        # tmpdir = makenewdir(name="tmp", fullname="temporary")
        logging.debug("IO demultiplexer settings: OK")
        debug()

        logging.info("Counting the total number of reads")
        totalreads1 = count_reads(FileList[0])
        if reverse:
            totalreads2 = count_reads(FileList[1])
            if totalreads1 != totalreads2:
                logging.error("Number of reads in paired do not match")
                raise EXONtoolsError
            logging.info("{0:d} paired-end reads are found within input files".format(totalreads1))
        else:
            logging.info("{0:d} single-end reads are found within input file".format(totalreads1))
        debug()

        # RUN ANALYSIS
        if pattern:
            stats_collector = findpatterns(FileList, output.path, pattern, demultiplexer.suffix, gzoutput, totalreads1)
        if barcode and not inseqsearch:
            stats_collector = findindex(barcodes, FileList, output.path, demultiplexer.suffix, tolerance, gzoutput, totalreads1)
        if barcode and inseqsearch:
            stats_collector = findbarcode(barcodes, FileList, output.path, demultiplexer.suffix, tolerance, inseqsearch, trim, start, gzoutput, totalreads1)
        debug()

        # SAVE STATS
        savestats(stats_collector, output.path)

        # READ QUALITY TESTS
        runqc(rqc, outpath=output.path, message="Demultiplexing analysis", gzout=gzoutput)


def debug():
    """debuger"""
    if demultiplexer.debug:
        pdb.set_trace()


def demux_checks(barcode, pattern, inseqsearch, gzoutput, start):
    """PERFORM ALL CHECKS"""
    barcodes = None
    logging.debug("Performing input parameter checks")
    if demultiplexer.extra:
        logging.warning("Extra string argument cannot be used in 'demultiplexer' program and therefore will be ommited")
    if gzoutput:
        logging.warning("All output FASTQ files will be compressed with gzip")
    if barcode and pattern:
        logging.error("'--pattern' and '--barcode' flags are mutually exclusive. You cannot use them together in one run")
        raise EXONtoolsError("Demultiplexer command line error")
    elif barcode:
        logging.info("Barcode option was enabled for demultiplexing")
        barcodes = get_barcodes(barcode)
        if inseqsearch:
            logging.warning("Library barcodes will be searched within read sequences")
            if '+' in (list(barcodes.keys())[0]):
                logging.error("Dual indexes cannot be used for inread search")
                logging.error("Please correct your input command and try again")
                raise EXONtoolsError("Demultiplexer command line error")
        else:
            if start > 1:
                logging.warning("Start position cannot be changed if barcodes are not in sequence")
    elif pattern:
        logging.info("Library name pattern option was enabled")
        if pattern == "default":
            pattern = "^@([^_:]+)[_:].*$"
        if pattern.count("(") != 1 or pattern.count(")") != 1:
            logging.error("Capture a single pattern using round parentheses in your RegExp")
            raise EXONtoolsError("GREP pattern error")
        logging.info("GREP pattern '{0}' will be used for demultiplexing".format(pattern))
    else:
        logging.error("You must provide either '--barcode' or '--pattern' option to start the demultiplexer program")
        raise EXONtoolsError("Demultiplexer command line error")
    logging.debug("Demultiplexer checks: OK")
    return barcodes, pattern


def demux_io(forward, reverse):
    """Input path parsing"""
    logging.debug("Performinng input path parsing")
    getinput.format(['.fq', '.fastq', '.fastq.gz', '.fq.gz'])
    paired, unpaired = parseinput(inpath=None, forward=forward, reverse=reverse)
    if unpaired:
        FileList = list(unpaired.values())[0]
        logging.info("The following file with unpaired sequencing reads will be demultiplexed")
        logging.info(FileList[0])
    elif paired:
        FileList = list(paired.values())[0]
        logging.info("The following files with paired sequencing reads will be demultiplexed")
        logging.info(FileList[0])
        logging.info(FileList[1])
    else:
        logging.error("Demultiplexer program can only work with single FASTQ files (forward & reverse)")
        raise EXONtoolsError("I/O error")
    logging.debug("Input path parsing: OK")
    return FileList


def savestats(results, outpath):
    """Save stats for demultiplexer"""
    if demultiplexer.stats and not demultiplexer.dryrun:
        statdir = makenewdir(name=os.path.join(outpath, "STATS"), fullname="STATS")
        logging.info("Demultiplexing stats will be saved to 'STATS/demux_stats.csv'")
        header = ["No", "Library", "NReads"]
        undefined_count = results.pop('undefined_count')
        with open(os.path.join(statdir.path, "demux_stats.csv"), 'w') as statfile:
            csv_writer = csv.writer(statfile)
            csv_writer.writerow(header)
            for i, libname in enumerate(sorted(results.keys(), key=natural_sort)):
                csv_writer.writerow([i + 1, libname, results[libname]])
            csv_writer.writerow([i + 2, 'undefined_count', undefined_count])
        logging.debug("Demultiplexing stats were successfully written to the file: OK")
        debug()


def runqc(rqc, outpath, message, gzout):
    """Run QC analysis"""
    if rqc and not demultiplexer.stats and not demultiplexer.dryrun:
        statdir = makenewdir(name=os.path.join(outpath, "STATS"), fullname="STATS").path
    else:
        statdir = os.path.join(outpath, "STATS")
    if rqc and gzout and not demultiplexer.dryrun:
        rqc_test(outpath, statdir, demultiplexer.threads, extension="*.gz", comment=message)
        debug()
    elif rqc and not demultiplexer.dryrun:
        rqc_test(outpath, statdir, demultiplexer.threads, comment=message)
        debug()
    else:
        pass


def get_barcodes(barcpath):
    """
    This function parses input file with barcodes
    and returns the dictionary:
    keys = barcode sequences
    values = library names
    """
    logging.debug("Reading file with library barcodes")
    barcpath = os.path.abspath(barcpath)
    try:
        with open(barcpath, 'r') as barcfile:
            logging.info("Parsing the file with barcodes:")
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


def opengzfile(inpath, gzout=False, mode='a'):
    """Open file"""
    if gzout:
        return gzip.open(inpath + ".gz", mode + 't')
    else:
        return open(inpath, mode)


def writeindexread(outpath, libname, read, suffix, gzout):
    """Write read by index"""
    file_name = os.path.join(outpath, libname + suffix + ".fq")
    outfile = opengzfile(file_name, gzout)
    outfile.write("@{0}\n{1}\n+{2}\n{3}\n".format(read.identifier, read.seq, read.info, read.qual))
    outfile.close()


def writebacrread(outpath, libname, barcode, read, suffix, trim, gzout):
    """Write read by barcode"""
    readseq = read.seq[len(barcode) + trim:]
    qualseq = read.qual[len(barcode) + trim:]
    if read.identifier.endswith(':'):
        identifier = rread.identifier + barcode
    else:
        identifier = read.identifier
    file_name = os.path.join(outpath, libname + suffix + ".fq")
    outfile = opengzfile(file_name, gzout)
    outfile.write("@{0}\n{1}\n+{2}\n{3}\n".format(identifier, readseq, read.info, qualseq))
    outfile.close()


def findindex(barcodes, filepath, outpath, suffix, tolerance, gzout, total):
    """
    Search barcodes in read names
    Can work with dual indexes
    """
    logging.info("Demultiplexing reads by index search within read name")
    if len(filepath) == 2:
        suffix1 = '_R1' + suffix
        suffix2 = '_R2' + suffix
    else:
        suffix1 = suffix
    total_defined = 0
    undefined_count = 0
    barc_stats = {x: 0 for x in list(barcodes.values())}                 # list all libraries and add 0 counts

    if not demultiplexer.dryrun:
        infileseq1 = SeqIO(filepath[0], fileformat="FASTQ")                  # open file with forward reads
        gen1 = infileseq1.read()
        read1 = next(gen1)
        if len(filepath) > 1:                                                # open file with reverse reads
            infileseq2 = SeqIO(filepath[1], fileformat="FASTQ")
            gen2 = infileseq2.read()
            read2 = next(gen2)
        else:
            infileseq2 = None

        # PERFORM CHECKS FOR EACH READ
        while(infileseq1.total <= total):
            if read1.barcode:
                selected_barc = []
                if set(list(read1.barcode)) - {"A", "T", "G", "C", "+"}:            # solve Ns in barcodes
                    barc_list = replaceN(read1.barcode)
                    for barc in barc_list:
                        if barc in barcodes:
                            selected_barc.append(barcodes[barc])            # find corresponding barcodes
                else:
                    for tolr in range(tolerance + 1):
                        barc_list = tolerate_barcode(read1.barcode, tolr)   # alternatively tolerate barcodes
                        for barc in barc_list:
                            if barc in barcodes:
                                selected_barc.append(barcodes[barc])        # find corresponding barcodes
                        if selected_barc:
                            break

                # SAVE READS WITH A SINGLE MATCH
                if len(selected_barc) == 1:
                    barc_stats[selected_barc[0]] += 1
                    total_defined += 1
                    writeindexread(outpath, selected_barc[0], read1, suffix1, gzout)
                    if infileseq2:
                        writeindexread(outpath, selected_barc[0], read2, suffix2, gzout)

                # TRY TO SOLVE BARCODE BOTH WITH N and TOLERANCE
                elif len(selected_barc) == 0 and read1.barcode.count('N') == 1 and tolerance > 0:
                    for barc in barc_list:
                        for tolr in range(tolerance + 1):
                            new_barc_list = tolerate_barcode(barc, tolr)
                            for barc in new_barc_list:
                                if barc in barcodes:
                                    selected_barc.append(barcodes[barc])       # find corresponding barcodes
                            if selected_barc:
                                break
                    if len(selected_barc) == 1:
                        barc_stats[selected_barc[0]] += 1
                        total_defined += 1
                        writeindexread(outpath, selected_barc[0], read1, suffix1, gzout)
                        if infileseq2:
                            writeindexread(outpath, selected_barc[0], read2, suffix2, gzout)
                    else:
                        undefined_count += 1
                        writeindexread(outpath, "Undefined_reads", read1, suffix1, gzout)
                        if infileseq2:
                            writeindexread(outpath, "Undefined_reads", read2, suffix2, gzout)

                # TRY TO SOLVE BARCODE WITH A SECOND READ
                elif infileseq2:
                    selected_barc = []
                    if set(list(read2.barcode)) - {"A", "T", "G", "C", "+"}:       # solve Ns in barcodes
                        barc_list = replaceN(read2.barcode)
                        for barc in barc_list:
                            if barc in barcodes:
                                selected_barc.append(barcodes[barc])        # find corresponding barcodes
                    else:
                        for tolr in range(tolerance + 1):
                            barc_list = tolerate_barcode(read2.barcode, tolr)    # alternatively tolerate barcodes
                            for barc in barc_list:
                                if barc in barcodes:
                                    selected_barc.append(barcodes[barc])        # find corresponding barcodes
                            if selected_barc:
                                break
                    if len(selected_barc) == 1:
                        barc_stats[selected_barc[0]] += 1
                        total_defined += 1
                        writeindexread(outpath, selected_barc[0], read1, suffix1, gzout)
                        writeindexread(outpath, selected_barc[0], read2, suffix2, gzout)
                    else:
                        undefined_count += 1
                        writeindexread(outpath, "Undefined_reads", read1, suffix1, gzout)
                        writeindexread(outpath, "Undefined_reads", read2, suffix2, gzout)

                # SAVES UNMATCHED READS or READS WITH MULTIPLE MATCHES
                else:
                    undefined_count += 1
                    writeindexread(outpath, "Undefined_reads", read1, suffix1, gzout)
                    if infileseq2:
                        writeindexread(outpath, "Undefined_reads", read2, suffix2, gzout)
            else:
                undefined_count += 1
                writeindexread(outpath, "Empty_index", read1, suffix1, gzout)
                if infileseq2:
                    writeindexread(outpath, "Empty_index", read2, suffix2, gzout)

            printProgressBar(infileseq1.total, total, prefix="", suffix="", decimals=1, length=50, fill='█', printEnd="\r")
            try:
                read1 = next(gen1)
                if infileseq2:
                    read2 = next(gen2)
            except StopIteration:
                break
    logging.info("Read demultiplexing is done: ███████████████| 100.0 %")
    time.sleep(0.5)

    barc_stats["undefined_count"] = undefined_count
    try:
        logging.info("Total {0:d} reads were demultiplexed from total {1:d} reads ({2:0.2f} %) ".format(total_defined, total, total_defined / total * 100))
    except ZeroDivisionError:
        logging.info("Total 0 reads were demultiplexed from total 0 reads")
    return barc_stats


def findbarcode(barcodes, filepath, outpath, suffix, tolerance, inseq, trim, start, gzout, total):
    """
    Demultiplex reads by barcodes inserted within sequences
    Cannot work with dual indexes
    """
    logging.info("Demultiplexing reads by barcode search within sequence")

    if len(filepath) == 2:
        suffix1 = '_R1' + suffix
        suffix2 = '_R2' + suffix
    else:
        suffix1 = suffix
    total_defined = 0
    undefined_count = 0
    barc_stats = {x: 0 for x in list(barcodes.values())}
    barcrange = (min(map(len, list(barcodes.keys()))), max(map(len, list(barcodes.keys()))))

    if not demultiplexer.dryrun:
        infileseq1 = SeqIO(filepath[0], fileformat="FASTQ")                  # open file with forward reads
        gen1 = infileseq1.read()
        read1 = next(gen1)
        if len(filepath) > 1:                                                # open file with reverse reads
            infileseq2 = SeqIO(filepath[1], fileformat="FASTQ")
            gen2 = infileseq2.read()
            read2 = next(gen2)
        else:
            infileseq2 = None

        # PERFORM CHECKS FOR EACH READ
        while(infileseq1.total <= total):
            selected_barc = []
            for tolr in range(tolerance + 1):
                barc_list = []
                for i in range(barcrange[0], barcrange[1] + 1):
                    barc = read1.seq[start - 1 : start - 1 + i]
                    if barc.count('N') == 1:
                        newbarks = replaceN(barc)
                        for bb in newbarks:
                            barc_list = barc_list + tolerate_barcode(bb, tolr)
                    else:
                        barc_list = barc_list + tolerate_barcode(barc, tolr)

                for barc in barc_list:
                    if barc in barcodes:
                        selected_barc.append((barcodes[barc], barc))            # find corresponding barcodes
                if selected_barc:                                               # break tolerance loop search
                    break

            # SAVE READS WITH A SINGLE MATCH
            if len(selected_barc) == 1:
                libname = selected_barc[0][0]
                barcseq = selected_barc[0][1]
                barc_stats[libname] += 1
                total_defined += 1
                writebacrread(outpath, libname, barcseq, read1, suffix1, trim, gzout)
                if infileseq2:
                    writebacrread(outpath, libname, barcseq, read2, suffix2, trim, gzout)
            else:
                undefined_count += 1
                writeindexread(outpath, "Undefined_reads", read1, suffix1, gzout)
                if infileseq2:
                    writeindexread(outpath, "Undefined_reads", read2, suffix2, gzout)

            printProgressBar(infileseq1.total, total, prefix="", suffix="", decimals=1, length=50, fill='█', printEnd="\r")
            try:
                read1 = next(gen1)
                if infileseq2:
                    read2 = next(gen2)
            except StopIteration:
                break
    logging.info("Read demultiplexing is done: ███████████████| 100.0 %")
    time.sleep(0.5)

    barc_stats["undefined_count"] = undefined_count
    try:
        logging.info("Total {0:d} reads were demultiplexed from total {1:d} reads ({2:0.2f} %) ".format(total_defined, total, total_defined / total * 100))
    except ZeroDivisionError:
        logging.info("Total 0 reads were demultiplexed from total 0 reads")
    return barc_stats


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
        final_seqs = []
        for i, nuc in enumerate(seq):                    # check strange
            if nuc not in nuc_list:
                for k in nuc_list:
                    new_seqs.append(seq[:i] + k + seq[i + 1:])
        for newbarc in new_seqs:
            if 'N' in newbarc:
                final_seqs = final_seqs + replaceN(newbarc)
            else:
                final_seqs.append(newbarc)
        if final_seqs:                                  # append new seqs if found
            final_seqs = list(set(final_seqs))
            barcode_list.append(final_seqs)
        else:
            barcode_list.append([seq])                   # append old seq if OK
    final_list = list(itertools.product(*barcode_list))  # make all combination pairs
    result = ["+".join(x) for x in final_list]           # produce result with joint barcodes
    return result


def findpatterns(filepath, outpath, pattern, suffix, gzout, total):
    """Demultiplex reads using GREP pattern in read name identifier"""
    logging.info("Demultiplexing reads by pattern search within read names")
    if len(filepath) == 2:
        suffix1 = '_R1' + suffix
        suffix2 = '_R2' + suffix
    else:
        suffix1 = suffix
    total_defined = 0
    undefined_count = 0
    demux_stats = {}

    if not demultiplexer.dryrun:
        infileseq1 = SeqIO(filepath[0], fileformat="FASTQ")                  # open file with forward reads
        gen1 = infileseq1.read()
        read1 = next(gen1)
        if len(filepath) > 1:                                                # open file with reverse reads
            infileseq2 = SeqIO(filepath[1], fileformat="FASTQ")
            gen2 = infileseq2.read()
            read2 = next(gen2)
        else:
            infileseq2 = None

        # PERFORM CHECKS FOR EACH READ
        while(infileseq1.total <= total):
            libname = check_pattern(pattern, "@" + read1.name + " " + read1.extra)
            if libname:
                total_defined += 1
                try:
                    demux_stats[libname] += 1
                except KeyError:
                    demux_stats[libname] = 1
                writeindexread(outpath, libname, read1, suffix1, gzout)
                if infileseq2:
                    writeindexread(outpath, libname, read2, suffix2, gzout)
            else:
                undefined_count += 1
                writeindexread(outpath, "Undefined_reads", read1, suffix1, gzout)
                if infileseq2:
                    writeindexread(outpath, "Undefined_reads", read2, suffix2, gzout)

            printProgressBar(infileseq1.total, total, prefix="", suffix="", decimals=1, length=50, fill='█', printEnd="\r")
            try:
                read1 = next(gen1)
                if infileseq2:
                    read2 = next(gen2)
            except StopIteration:
                break
    logging.info("Read demultiplexing is done: ███████████████| 100.0 %")
    time.sleep(0.5)

    demux_stats["undefined_count"] = undefined_count
    try:
        logging.info("Total {0:d} reads were demultiplexed from total {1:d} reads ({2:0.2f} %) ".format(total_defined, total, total_defined / total * 100))
    except ZeroDivisionError:
        logging.info("Total 0 reads were demultiplexed from total 0 reads")
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
    return result[0]


def count_reads(inpath):
    """Count reads by number of lines"""
    infileseq = SeqIO(inpath, fileformat="FASTQ")
    infileseq.totalcount(lines=True)
    if infileseq.total % 4 != 0:
        logging.error("FASTQ format error. Please check your input file.")
        raise EXONtoolsError("Wrong fastq format")
    return infileseq.total // 4
