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
import subprocess

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
        self.deduplicate(args.inpath, args.forward, args.reverse, args.unpaired, args.outdir, args.gzoutput, args.skip, args.phred, args.rqc)

    def deduplicate(self, inpath, forward, reverse, unpaired_in, outdir, gzoutput, skip, phred, rqc):

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
        dupldir = makenewdir(name="DUPLICATES", fullname="DUPLICATES")
        tmpdir = makenewdir(name="tmp", fullname="temporary")
        logging.debug("IO settings: OK")
        debug()

        # SET DEDUP PARAMETERS
        deduppars(skip, gzoutput)
        debug()
 
        # RUN DEDUPLICATION        
        TASKS = []
        if paired:
            trunc = truncreads(skip,True)
            for lib in sorted(paired.keys(), key=natural_sort):
                TASKS.append(worker(dedup, [lib, tmpdir.path, trunc, gzoutput, deduplicator.suffix, paired[lib][0], paired[lib][1],phred]))
        if unpaired:
            trunc = truncreads(skip,False)
            for lib in sorted(unpaired.keys(), key=natural_sort):
                TASKS.append(worker(dedup, [lib, tmpdir.path, trunc, gzoutput, deduplicator.suffix, unpaired[lib][0],None,phred]))

        results = runtask(TASKS, "Deduplicator")
        logging.debug("Deduplication analysis succesfully finished: OK")
        debug()

        # SAVE STATS
        savestats(results, output.path)

        # READ QUALITY TESTS
        runqc(rqc, outpath=output.path, message="Deduplicated reads", gzout=gzoutput)

        if not deduplicator.keeptmp:
            tmpdir.delete()


def runtask(TASKS, message):
    """Run tasks"""
    if TASKS:
        processes_requested = set_threads(message, len(TASKS), deduplicator.threads)
        pool = create_pool(processes_requested)
        jobs = hard_worker(run_instance, TASKS, pool)
        close_pool(pool)
        logging.debug("TASKs running process: OK")
        return parsejobs(jobs)


def parsejobs(jobs):
    """Convert jobs from list to dictionary"""
    if not jobs:
        logging.error("Multiprocessing produced no results")
        raise EXONtoolsError("Multiprocessing error")
    job_collector = {}
    for result in jobs:
        if result:
            job_collector.update(result)
    if not job_collector:
        logging.warning("Multiprocessing produced no results")
    return job_collector


def meanphred(inline, phred = 33):
    """Calculate mean phred score for input quality line"""
    try:
        phredstr = [ord(x)- phred for x in inline]
        return sum(phredstr)*100//len(phredstr)
    except ZeroDivisionError:
        logging.error("No quality line is provided for PHRED conversion")
        raise EXONtoolsError


def external_sort(tmpfile):
    """Use bash sort function"""
    list_dup = subprocess.Popen("cat "+tmpfile+" | sort -k1,1r -k2nr -T " + os.path.dirname(tmpfile) +" | awk '{if (NR>1 && length($1)>length(f) && index($1,f)>0) {print name} else if (NR>1 && length($1)==length(f) && (index($1,f)>0)) {print $3}} {f=$1} {name =$3}'", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    (results,err) = list_dup.communicate()
    if err != '':
        logging.error("EXONtools ERROR (subprocess failed)")
        logging.error(err)
        raise EXONtoolsError
    results = results.split()
    duplicates = [int(x.replace("LINE","")) for x in results]
    duplicates.sort()
    return duplicates

        
def opengzfile(inpath, gzout=False, mode='w'):
    if gzout:
        return gzip.open(inpath + ".gz", mode + 't')
    else:
        return open(inpath, mode)


def dedup(sample, tmppath, skip, gzoutput, suffix, inpathR1, inpathR2=None, phred=33):
    """DEDUPLICATION"""

    # MAKE TEMP FILE
    if inpathR2:
        logging.info("Deduplicating paired reads for '{0:s}' library".format(sample))
        outpath = os.path.join(tmppath,sample+"_paired.out")
        infile1 = SeqIO(inpathR1,fileformat="FASTQ")
        infile2 = SeqIO(inpathR2,fileformat="FASTQ")
        with open(outpath,'w') as outfile:
            for read1,read2 in zip(infile1.read(),infile2.read()):
                outfile.write("{0:s} {1:d} LINE{2:d}\n".format(
                    read1.seq[skip[0]:skip[1]]+read2.seq[skip[2]:skip[3]],
                    meanphred(read1.qual[skip[0]:skip[1]]+read2.qual[skip[2]:skip[3]],phred=phred),
                    infile1.total))
    else:
        logging.info("Deduplicating unpaired reads for '{0:s}' library".format(sample))
        outpath = os.path.join(tmppath,sample+"_unpaired.out")
        infile1 = SeqIO(inpathR1,fileformat="FASTQ")
        with open(outpath,'w') as outfile:
            for read1 in infile1.read():
                outfile.write("{0:s} {1:d} LINE{2:d}\n".format(
                    read1.seq[skip[0]:skip[1]],
                    meanphred(read1.qual[skip[0]:skip[1]],phred=phred),
                    infile1.total))
     
    # FIND DUPLICATES FROM TEMPORARY FILE
    duplicates = external_sort(outpath)
    NDUPS= len(duplicates)
    logging.info("{0:d} duplicates are found in '{1:s}' library".format(NDUPS,sample))

    # REMOVE TEMPORARY FILE
    os.remove(outpath)

    # SAVE RESULTS
    outdir = os.path.dirname(tmppath)
    dupldir = os.path.join(outdir,"DUPLICATES")
    ndup = 0
    if inpathR2:
        with open(os.path.join(dupldir,sample+"_duplicates.dat"), 'w') as duplfile:
            outpath1 = os.path.join(outdir,sample+"_R1"+suffix+".fq")
            outpath2 = os.path.join(outdir,sample+"_R2"+suffix+".fq")
            outfile1 = opengzfile(outpath1,gzout=gzoutput)
            outfile2 = opengzfile(outpath2,gzout=gzoutput)
            infile1 = SeqIO(inpathR1,fileformat="FASTQ")
            infile2 = SeqIO(inpathR2,fileformat="FASTQ")
            for read1,read2 in zip(infile1.read(),infile2.read()):
                if ndup < NDUPS and infile1.total == duplicates[ndup]:
                    duplfile.write(read1.name+'\n')
                    ndup +=1
                else:
                    outfile1.write("@{0:s}\n{1:s}\n+\n{2:s}\n".format(read1.identifier,read1.seq,read1.qual))
                    outfile2.write("@{0:s}\n{1:s}\n+\n{2:s}\n".format(read2.identifier,read2.seq,read2.qual))
            return {sample: [infile1.total,NDUPS]}
    else:
        with open(os.path.join(dupldir,sample+"_duplicates.dat"), 'w') as duplfile:
            outpath1 = os.path.join(outdir,sample+suffix+".fq")
            outfile1 = opengzfile(outpath1,gzout=gzoutput)
            infile1 = SeqIO(inpathR1,fileformat="FASTQ")
            for read1 in infile1.read():
                if ndup < NDUPS and infile1.total == duplicates[ndup]:
                    duplfile.write(read1.name+'\n')
                    ndup +=1
                else:
                    outfile1.write("@{0:s}\n{1:s}\n+\n{2:s}\n".format(read1.identifier,read1.seq,read1.qual))
            return {sample: [infile1.total,NDUPS]}


def debug():
    """debuger"""
    if deduplicator.debug:
        pdb.set_trace()


def savestats(results, outpath):
    """Save all stats"""
    if deduplicator.stats and not deduplicator.dryrun:
        statdir = makenewdir(name=os.path.join(outpath, "STATS"), fullname="STATS")
        logging.info("Read stats will be saved to 'STATS/deduplicate_stats.csv'")
        header = ["No", "LIBRARY", "#INITIAL_READS", "#DUPLICATES",  "#PERCENTAGE", "#FINAL_READS"]
        with open(os.path.join(statdir.path, "deduplicate_stats.csv"), 'w') as statfile:
            csv_writer = csv.writer(statfile)
            csv_writer.writerow(header)
            for i, lib in enumerate(sorted(results.keys(), key=natural_sort)):
                csv_writer.writerow([
                    i,
                    lib,
                    results[lib][0],
                    results[lib][1],
                    round(results[lib][1]/results[lib][0]*100,2),
                    results[lib][0] - results[lib][1]
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


def truncreads(skip, pair=False):
    """skip conformation"""
    trunc = skip[:]
    logging.debug("Changing skip vector")
    for i, x in enumerate(trunc):
        if x == 0:
            trunc[i] = None
        elif i > 1 and not pair:
            trunc[i] = None
        elif i % 2 == 1:
            trunc[i] = -x
    logging.debug("Changing skip vector: OK")
    return trunc
