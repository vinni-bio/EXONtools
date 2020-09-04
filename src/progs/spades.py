#! /usr/bin/env python
# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from __future__ import print_function, division
import pdb
import os
import logging
import shutil
import csv
import re

from mains.EXT_prog import EXTprogram
from mains.EXT_executor import executor
from mains.EXT_parallel import run_executor
from mains.EXT_errors import EXONtoolsError
from mains.EXT_IO import getinput, output, makenewdir, parseinput
from progs.assasses import metrics
from utils.sorting import natural_sort


class spades(EXTprogram):
    """SPADES program assembles sequencing reads into contigs for single cell sequencing data"""

    name = "spades"

    def execute_program(self):
        args = self.args
        self.assemble(args.inpath, args.forward, args.reverse, args.unpaired, args.outdir, args.kmers, args.coverage, args.npar, args.mpirun, args.cores, args.parallel)

    def assemble(self, inpath, forward, reverse, unpaired_in, outdir, kmers, coverage, npar, mpirun, cores, parallel):

        if spades.debug:
            pdb.set_trace()

        executor.setconfig("spades")

        if spades.debug:
            pdb.set_trace()

        # TESTING INPUT PARAMETER VALUES
        params = spades_pars(kmers, coverage, npar, mpirun, cores, parallel)

        if spades.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        spades.run_dry(getinput, output, makenewdir, executor)
        spades.set_debug(executor)

        # GET INPUT FILES
        getinput.format(['.fq', '.fastq'])
        paired, unpaired = parseinput(inpath, forward, reverse)
        if unpaired_in and not inpath:
            unpaired = {os.path.basename(unpaired_in).split("_")[0].split(".")[0] + "_unpaired": getinput(unpaired_in).files}
        elif unpaired_in and inpath:
            logging.error("Input arguments '-i' and '-U' cannot be used simultaneously")
            raise EXONtoolsError("Input path error")
        else:
            pass

        librs = list(set([x.split("_")[0] for x in list(paired.keys()) + list(unpaired.keys())]))
        librs.sort(key=natural_sort)

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE TMP DIRECTORIES
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        # MAKE LOG DIR
        log_dir = makenewdir(name="spades_logs", fullname="LOG")
        logging.info("Find all SPAdes log files in {0:s} directory".format(os.path.basename(log_dir.path)))

        logging.debug("IO settings: OK")

        # RUNNING THE ASSEMBLY ANALYSIS
        outpaths = []

        if spades.debug:
            pdb.set_trace()

        for lib in librs:

            if lib + "_paired" in paired and lib + "_unpaired" in unpaired:
                infiles = "-1 {0:s} -2 {1:s} -s {2:s}".format(paired[lib + "_paired"][0], paired[lib + "_paired"][1], unpaired[lib + "_unpaired"][0])
            elif lib + "_paired" in paired and not lib + "_unpaired" in unpaired:
                infiles = "-1 {0:s} -2 {1:s}".format(paired[lib + "_paired"][0], paired[lib + "_paired"][1])
            elif not lib + "_paired" in paired and lib + "_unpaired" in unpaired:
                infiles = "-s {0:s}".format(unpaired[lib + "_unpaired"][0])
            else:
                logging.error("You found the bug during parsing of SPAdes input files")
                raise EXONtoolsError("Input files command line error")

            libouttmp = os.path.join(tmpdir.path, lib + "_spades")
            makenewdir(libouttmp)

            logpath = os.path.join(log_dir.path, lib + "_spades.log")

            logging.info("Starting SPAdes assembly for '{0:s}' library".format(lib))

            run_executor(executor(
                program="SPAdes",
                params=[params, libouttmp, spades.threads, int(spades.memory), infiles],
                conditions={"positive": [spades.threads], "pathexists": [libouttmp]},
                custom_arg_string=spades.extra + " > /dev/null",
            ))

            # if "--rna" in spades.extra:
            #     outname = "transcripts.fasta"
            # else:
            #     outname = "contigs.fasta"

            if not spades.dryrun:

                spadelogpath = os.path.join(libouttmp, "spades.log")
                shutil.move(spadelogpath, logpath)

                kmerlist = [x for x in os.listdir(libouttmp) if os.path.isdir(os.path.join(libouttmp, x)) and re.search("K\d+", x)]

                for kmer in kmerlist:
                    kmerin = os.path.join(libouttmp, kmer, "simplified_contigs.fasta")
                    kmerout = os.path.join(output.path, lib + "_km" + kmer[1:] + spades.suffix + ".fa")
                    try:
                        shutil.move(kmerin, kmerout)
                        outpaths.append(kmerout)
                    except (IOError, OSError):
                        logging.error("SPAdes assembly output file for the library '{0:s}' is not found".format(lib))
                        raise EXONtoolsError("Non-existing path error")

            logging.debug("SPAdes assembly for the library '{0:s}' is finished: OK".format(lib))

        logging.debug("Successfully finished all assemblies: OK")

        if spades.debug:
            pdb.set_trace()

        # SAVE STATS
        if spades.stats and not spades.dryrun:
            spades_stats(outpaths)

        # DELETE TEMPORARY DIRECTORY
        if not spades.keeptmp:
            tmpdir.delete()


def spades_pars(kmers, coverage, npar, mpirun, cores, parallel):
    """Checking spades input parameters"""

    logging.debug("Testing SPAdes input settings")

    params = "--only-assembler"

    if spades.extra:
        logging.warning("The following arguments will be added to 'SPAdes' command line:")
        logging.warning(spades.extra)

    if kmers:
        if kmers == "program default":
            logging.warning("Using the program default kmer settings")
        else:
            try:
                kmers = [int(x) for x in set(kmers)]
            except ValueError:
                logging.error("K-mer values must be integers")
                raise EXONtoolsError("K-mer format error")
            # if len(kmers) > 1 and "--rna" in spades.extra:
            #     logging.error("You cannot specify multiple k-mer sizes in RNA-Seq mode!")
            #     raise EXONtoolsError("SPAdes parameter error")
            for k in kmers:
                if not isinstance(k, int) or k < 10 or not k % 2:
                    logging.error("All k-mers must be odd integers greater than 10 and less than 128 bp")
                    raise EXONtoolsError("Wrong format of k-mer length")
            kmers.sort()
            logging.warning("SPAdes analysis will use the following k-mers: {0}".format(kmers))
            params = params + " -k {0:s}".format(",".join(map(str, kmers)))
        logging.debug("K-mer settings check: OK")
    else:
        logging.error("No k-mers were provided in -k option")
        raise EXONtoolsError("No k-mer options were provided")

    if cores != 1:
        logging.error("The number of cores for SPAdes run must be set by -T/--threads general option")
        raise EXONtoolsError("Wrong multithread settings")

    if coverage:
        try:
            if coverage == 'program default':
                logging.warning("The current run is not using the coverage cutoff setting")
            elif coverage[0] == 'auto':
                logging.warning("The current run is using the AUTO coverage cutoff setting")
                params = params + " --cov-cutoff auto"
            elif coverage[0] == 'off':
                logging.warning("The current run is not using the coverage cutoff")
                params = params + " --cov-cutoff off"
            elif len(coverage) != 1 or int(coverage[0]) <= 0:
                logging.error("SPAdes coverage settings must be: 'auto', 'off' or a single positive integer")
                raise EXONtoolsError("Wrong settings for coverage parameter")
            else:
                params = params + " --cov-cutoff {0:d} ".format(int(coverage[0]))
                logging.info("The coverage threshold for SPAdes analysis is {0:d}".format(int(coverage[0])))
                logging.debug("Checking kmer coverage parameter: OK")
        except ValueError:
            logging.error("SPAdes assembly analysis accepts only a single positive INTEGER value for the kmer coverage cutoff")
            raise EXONtoolsError("Coverage format error")

    if not isinstance(npar, int) or npar != 5:
        logging.warning("SPAdes assembly does is not using the paired read coverage option (-n/--nreads)")

    if mpirun:
        logging.warning("SPAdes assembly analysis is not using the MPI. The mpirun option will be ignored")

    if parallel:
        logging.warning("SPAdes assembly analysis is not using the MPI. The number of nodes is set to zero")

    logging.debug("SPAdes program settings: OK")
    return params


def spades_stats(paths):
    """Write stats for SPAdes analysis"""

    logging.debug("Collecting statistical data on assemblies")

    statpath = os.path.join(output.path, "STATS_spades_assembly.csv")

    logging.info("All assembly metrics will be saved to '{0:s}'".format(os.path.basename(statpath)))

    header = ["Library", "N", "N>500", "N>1000", "N>5000", "N>10000", "LongestID", "Min", "N75", "N50", "N25", "Max", "L50", "L75", "NG50", "Mean", "E-size", "Total", "GC", "N's", "N's per 100kb"]

    with open(statpath, 'w') as statout:
        csv_writer = csv.writer(statout)
        csv_writer.writerow(header)

        for path in sorted(paths, key=natural_sort):
            x = metrics(path)
            csv_writer.writerow([x.ID.replace(spades.suffix, ""), x.totnum, x.long500, x.long1000, x.long5000, x.long10000, x.longest, x.min, x.N75, x.N50, x.N25, x.max, x.L50, x.L75, x.NG, x.mean, x.esize, x.totlen, x.gc, x.nbases, x.n100])

    logging.debug("Stats data collected and written to a file: OK")
