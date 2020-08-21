# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import pdb
import os
import logging
import shutil
import csv


from mains.EXT_prog import EXTprogram
from mains.EXT_executor import executor
from mains.EXT_parallel import run_executor
from mains.EXT_errors import EXONtoolsError
from mains.EXT_validator import import_module_test
from mains.EXT_IO import getinput, parseinput, output, makenewdir
from progs.assasses import metrics
from utils.sorting import natural_sort


class trinity(EXTprogram):
    """Trinity program assembles sequencing reads into contigs for transcriptome data"""

    name = "trinity"

    def execute_program(self):
        args = self.args
        self.assemble(args.inpath, args.forward, args.reverse, args.unpaired, args.outdir, args.kmers, args.coverage, args.npar, args.mpirun, args.cores, args.parallel)

    def assemble(self, inpath, forward, reverse, unpaired_in, outdir, kmers, coverage, npar, mpirun, cores, parallel):

        if trinity.debug:
            pdb.set_trace()

        import_module_test("numpy")
        logging.debug("The python module 'numpy' is checked: OK")

        if trinity.debug:
            pdb.set_trace()

        executor.setconfig("trinity", "bowtie2", "samtools")

        if trinity.debug:
            pdb.set_trace()

        logging.debug("Exporting paths for Trinity dependencies")
        executor.exportpath("bowtie2", "samtools")
        os.environ['TRINITY_HOME'] = os.path.dirname(executor.program_paths['trinity'])
        logging.debug("Path exports: OK")

        if trinity.debug:
            pdb.set_trace()

        # TESTING INPUT PARAMETER VALUES
        params = trinity_pars(kmers, coverage, npar, mpirun, cores, parallel)

        if trinity.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        trinity.run_dry(getinput, output, makenewdir, executor)
        trinity.set_debug(executor)

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
        log_dir = makenewdir(name="trinity_logs", fullname="LOG")
        logging.info("Find all Trinity log files in {0:s}".format(os.path.basename(log_dir.path)))

        logging.debug("IO settings: OK")

        if trinity.debug:
            pdb.set_trace()

        # RUNNING THE ASSEMBLY ANALYSIS
        logging.debug("RUNNING TRINITY ASSEMBLY ANALYSES")
        outpaths = []
        for lib in librs:

            if lib + "_paired" in paired and lib + "_unpaired" in unpaired:
                newpath = os.path.join(tmpdir.path, os.path.basename(paired[lib + "_paired"][0]))
                shutil.copy2(paired[lib + "_paired"][0], newpath)
                with open(unpaired[lib + "_unpaired"][0], 'r') as oldfile:
                    with open(newpath, 'a') as newfile:
                        for line in oldfile:
                            newfile.write(line)
                infiles = "--left {0:s} --right {1:s}".format(newpath, paired[lib + "_paired"][1])
            elif lib + "_paired" in paired and not lib + "_unpaired" in unpaired:
                infiles = "--left {0:s} --right {1:s}".format(paired[lib + "_paired"][0], paired[lib + "_paired"][1])
            elif not lib + "_paired" in paired and lib + "_unpaired" in unpaired:
                infiles = "--single {0:s}".format(unpaired[lib + "_unpaired"][0])
            else:
                logging.critical("You found the bug during parsing of trinity input files")
                raise EXONtoolsError("Input files command line error")

            libouttmp = os.path.join(tmpdir.path, lib + "_trinity")
            makenewdir(libouttmp)

            logpath = os.path.join(log_dir.path, lib + "_trinity.log")

            logging.info("Starting Trinity assembly of '{0:s}' library".format(lib))
            if trinity.debug:
                pdb.set_trace()

            run_executor(executor(
                program="Trinity",
                params=[int(trinity.memory), params, trinity.threads, infiles, libouttmp],
                conditions={"positive": [trinity.threads], "pathexists": [libouttmp]},
                custom_arg_string=trinity.extra + " > " + logpath + " 2>&1"
            ))
            if not trinity.dryrun:
                trinpath = libouttmp + ".Trinity.fasta"
                outpath = os.path.join(output.path, lib + trinity.suffix + ".fa")
                try:
                    shutil.move(trinpath, outpath)
                    outpaths.append(outpath)
                except (IOError, OSError):
                    logging.error("The Trinity assembly output file for the library '{0:s}' is not found".format(lib))
                    raise EXONtoolsError("Non-existing path error")
            logging.debug("Trinity assembly for the library '{0:s}' is finished: OK".format(lib))

        # SAVE STATS
        if trinity.stats and not trinity.dryrun:
            trinity_stats(outpaths)

        # DELETE TEMPORARY DIRECTORY
        if not trinity.keeptmp:
            tmpdir.delete()


def trinity_pars(kmers, coverage, npar, mpirun, cores, parallel):
    """Checking trinity input parameters"""

    logging.debug("Testing trinity input settings")

    if trinity.extra:
        logging.warning("The following arguments will be added to 'Trinity' command line:")
        logging.warning(trinity.extra)

    params = ""
    if kmers:
        if kmers == "program default":
            logging.info("Using program default kmer settins")
            kmers = None
        else:
            try:
                kmers = [int(x) for x in set(kmers)]
            except ValueError:
                logging.error("K-mer values must be integers")
                raise EXONtoolsError("K-mer format error")
            if len(kmers) != 1:
                logging.error("Trinity assembly analysis accepts only a single value for the required kmer length")
                raise EXONtoolsError
            elif kmers[0] > 32 or kmers[0] < 10:
                logging.error("K-mer length for trinity assembly analysis cannot exceed 32 or to be lower than 10")
                raise EXONtoolsError("Wrong format of k-mer length")
            else:
                params = params + "--KMER_SIZE {0:d} ".format(kmers[0])
        logging.debug("Checking kmer length parameter: OK")

    if cores != 1:
        logging.error("The number of cores for Trinity run must be set by -T/--threads general option")
        raise EXONtoolsError("Wrong multithread settings")

    if coverage:
        try:
            if coverage == 'program default':
                logging.info("Trinity kmer coverage is set to 1 by using the default settings")
            elif len(coverage) != 1 or int(coverage[0]) <= 0:
                logging.error("Trinity assembly analysis accepts only a single positive integer for kmer coverage setting")
                raise EXONtoolsError("Wrong settings for coverage parameter")
            else:
                params = params + "--min_kmer_cov {0:d} ".format(int(coverage[0]))
                logging.info("The coverage threshold for Trinity analysis is {0:d}".format(int(coverage[0])))
                logging.debug("Checking kmer coverage parameter: OK")
        except ValueError:
            logging.error("Trinity assembly analysis accepts only a single positive INTEGER value for the required kmer coverage")
            raise EXONtoolsError("Coverage format error")

    if not isinstance(npar, int) or npar < 2:
        logging.error("The number of required paired reads for contig assembly in Trinity must be above or equal 2 (-n option)")
        raise EXONtoolsError("Paired read per k-mer input parameter format error")
    else:
        params = params + "--min_glue {0:d}".format(npar)
        logging.debug("Minimum number of paired reads for Trinity contig assembly is set: OK")

    if mpirun:
        logging.warning("Trinity assembly analysis is not using the MPI. The mpirun option will be ignored")

    if parallel:
        logging.warning("Trinity assembly analysis is not using the MPI. The number of nodes is set to zero")

    logging.debug("Trinity program settings: OK")
    return params


def trinity_stats(paths):
    """Write stats for Trinity analysis"""

    logging.debug("Collecting statistical data on assemblies")

    statpath = os.path.join(output.path, "STATS_trinity_assembly.csv")

    logging.info("All assembly metrics will be saved to '{0:s}'".format(os.path.basename(statpath)))

    header = ["Library", "N", "N>500", "N>1000", "N>5000", "N>10000", "LongestID", "Min", "N75", "N50", "N25", "Max", "L50", "L75", "NG50", "Mean", "E-size", "Total", "GC", "N's", "N's per 100kb"]

    with open(statpath, 'w') as statout:
        csv_writer = csv.writer(statout)
        csv_writer.writerow(header)

        for path in sorted(paths, key=natural_sort):
            x = metrics(path)
            csv_writer.writerow([x.ID.replace(trinity.suffix, ""), x.totnum, x.long500, x.long1000, x.long5000, x.long10000, x.longest, x.min, x.N75, x.N50, x.N25, x.max, x.L50, x.L75, x.NG, x.mean, x.esize, x.totlen, x.gc, x.nbases, x.n100])

    logging.debug("Stats data collected and written to a file: OK")
