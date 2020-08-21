# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import os
import logging
import csv
import pdb

from mains.EXT_prog import EXTprogram
from utils.gccheck import GC, Ncheck
from mains.EXT_errors import EXONtoolsError
from mains.EXT_IO import getinput, output
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, set_threads
from mains.EXT_worker import worker
from utils.sorting import natural_sort
from mains.EXT_validator import positive
from utils.seqIO import SeqIO


class assass(EXTprogram):
    """Estimate metrics of all assemblies within provided path"""

    name = "assass"

    def execute_program(self):
        args = self.args
        self.evaluate(args.inpath, args.outdir, args.minlen)

    def evaluate(self, inpath, outdir, minlen):

        if assass.debug:
            pdb.set_trace()

        assass.run_dry(getinput, output)

        positive([minlen], False)

        logging.info("Collecting stats for the assemblies in '{0:s}'".format(inpath))

        # GET INPUT FILES
        getinput.format(['.fa', '.fasta'])
        FileList = getinput(inpath).files

        # GET OUTPUT PATH
        output.nostrict()
        output(outdir)
        outpath = os.path.join(output.path, "assembly_stats.csv")

        logging.info("The assembly metrics will be saved to '{0:s}'".format(os.path.basename(outpath)))

        if os.path.exists(outpath):
            logging.warning("The file 'assembly_stats.csv' exists and will be overwritten")

        logging.debug("Input files and checks: OK")

        if assass.debug:
            pdb.set_trace()

        # COLLECT STATS
        logging.info("Starting to collect the assembly stats")
        stat_collector = []
        if not assass.dryrun:
            stat_collector.append(metrics(FileList[0], minlen))
            if len(FileList) > 1:
                TASKS = [worker(metrics, [x, minlen]) for x in FileList[1:]]
                processes_requested = set_threads("assass", len(TASKS), assass.threads, userdef=False)
                pool = create_pool(processes_requested)
                results = hard_worker(run_instance, TASKS, pool)
                close_pool(pool)
                stat_collector = stat_collector + results
                stat_collector.sort(key=lambda x: natural_sort(x.ID))

        logging.debug("All metrics were successfully collected")

        header = ["Library", "N", "N>500", "N>1000", "N>5000", "N>10000", "LongestID", "Min", "N75", "N50", "N25", "Max", "L50", "L75", "NG50", "Mean", "E-size", "Total", "GC", "N's", "N's per 100kb"]

        logging.debug("Ready to write stats to the output file")

        if assass.debug:
            pdb.set_trace()

        if not assass.dryrun:
            with open(outpath, 'w') as statout:
                csv_writer = csv.writer(statout)
                csv_writer.writerow(header)

                for x in stat_collector:
                    csv_writer.writerow([x.ID, x.totnum, x.long500, x.long1000, x.long5000, x.long10000, x.longest, x.min, x.N75, x.N50, x.N25, x.max, x.L50, x.L75, x.NG, x.mean, x.esize, x.totlen, x.gc, x.nbases, x.n100])

        logging.info("All assemblies were successfully evaluated")


class metrics(object):
    """Calculates various assembly metrics"""

    NG50 = None

    def __init__(self, path, seqminlen=200):

        if os.path.exists(path) and os.path.isfile(path):

            logging.info("Collecting stats for assembly '{0:s}'".format(os.path.basename(path)))

            fastafile = SeqIO(path)

            lengths = []
            gccont = []
            nbases = []
            minlen = 100000000000000
            maxlen = 0
            longestID = None
            totnum = 0

            for x in fastafile.read():
                totnum += 1
                xlen = len(x.seq)
                if xlen >= seqminlen:
                    lengths.append(xlen)
                    minlen = min(xlen, minlen)
                    maxlen = max(xlen, maxlen)
                    if xlen == maxlen:
                        longestID = x.name
                    gccont.append(GC(x.seq))
                    nbases.append(Ncheck(x.seq))
            logging.debug("Finished file parsing for stats: OK")

            if not lengths:
                logging.error("Found not contigs in {0:s}".format(path))
                raise EXONtoolsError("Assembly file error")

            self.ID = os.path.basename(os.path.splitext(fastafile.path)[0])
            self.totnum = totnum
            self.gc = round(sum(gccont) / len(gccont), 2)
            self.totlen = sum(lengths)
            self.min = minlen
            self.max = maxlen
            self.longest = longestID.split()[0]
            self.long500 = len([x for x in lengths if x >= 500])
            self.long1000 = len([x for x in lengths if x >= 1000])
            self.long5000 = len([x for x in lengths if x >= 5000])
            self.long10000 = len([x for x in lengths if x >= 10000])
            self.mean = int(round(sum(lengths) / len(lengths), 0))
            self.N25, _ = metrics.NL(lengths, perc=25)
            self.N50, self.L50 = metrics.NL(lengths, perc=50)
            self.N75, self.L75 = metrics.NL(lengths, perc=75)
            self.NG, _ = metrics.NL(lengths, NG=True)
            self.esize = int(round(sum([x * x / self.totlen for x in lengths]), 0))
            self.nbases = sum(nbases)
            self.n100 = int(round(100000 * self.nbases / self.totlen, 0))
        else:
            logging.error("Provided path does not exist")
            raise EXONtoolsError("Non-existing path error")

    @classmethod
    def NL(cls, lengths, perc=50, NG=False):

        if perc not in range(1, 100):
            logging.error("The percentage value must be within the range from 1 to 99")
            raise EXONtoolsError("NL method error")

        if isinstance(lengths, list):
            if lengths:
                if NG is True:
                    try:
                        half = cls.NG50 / 2
                    except TypeError:
                        half = sum(lengths) * perc / 100
                        cls.setNG50(lengths)
                else:
                    try:
                        half = sum(lengths) * perc / 100
                    except TypeError:
                        logging.error("Length list must contain only integer values")
                        raise EXONtoolsError("NL method error")

                lengths.sort(reverse=True)
                check = 0
                count = 0
                for l in lengths:
                    check = check + l
                    count += 1
                    if check > half:
                        return l, count
                return l, count
            else:
                logging.error("No lengths provided for NL calculation")
                raise EXONtoolsError("NL method error")
        else:
            logging.error("Lengths must be provided as a list")
            raise EXONtoolsError("NL method error")

    @classmethod
    def setNG50(cls, lengths):
        """Assings class NG as N50 for the first assembly"""

        try:
            cls.NG50 = sum(lengths)
        except TypeError:
            logging.error("Length list must contain only integer values")
            raise EXONtoolsError("NL method error")
