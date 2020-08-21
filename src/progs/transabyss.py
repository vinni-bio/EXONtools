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
import subprocess
import csv

from mains.EXT_prog import EXTprogram
from mains.EXT_executor import executor
from mains.EXT_validator import import_module_test
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, run_executor, set_threads
from mains.EXT_errors import EXONtoolsError
from mains.EXT_IO import getinput, output, makenewdir, parseinput
from progs.assasses import metrics
from utils.sorting import natural_sort


class transabyss(EXTprogram):
    """Tranabyss program assembles sequencing reads into contigs for transcriptome data"""

    name = "transabyss"

    def execute_program(self):
        args = self.args
        self.assemble(args.inpath, args.forward, args.reverse, args.unpaired, args.outdir, args.kmers, args.coverage, args.npar, args.mpirun, args.cores, args.parallel)

    def assemble(self, inpath, forward, reverse, unpaired_in, outdir, kmers, coverage, npar, mpirun, cores, parallel):

        if transabyss.debug:
            pdb.set_trace()

        import_module_test("igraph")
        logging.debug("The python module 'igraph' is checked: OK")

        executor.setconfig("transabyss", "abyss-pe", "blat")

        # EXPORTING PATHS
        logging.debug("Exporting paths for transABySS dependencies")
        executor.exportpath("blat", "abyss-pe")
        logging.debug("Path exports: OK")

        # TESTING INPUT PARAMETER VALUES
        mpi, kmers, coverage = transabyss_pars(kmers, coverage, npar, mpirun, cores, parallel)

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        transabyss.run_dry(getinput, output, makenewdir, executor)
        transabyss.set_debug(executor)

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

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE TMP DIRECTORIES
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        # MAKE LOG DIR
        log_dir = makenewdir(name="transabyss_logs", fullname="LOG")
        logging.info("Find all TransABySS log files in '{0:s}'".format(os.path.basename(log_dir.path)))

        logging.debug("IO settings: OK")

        # RUNNING THE ASSEMBLY ANALYSIS
        TASKS = transabyss_tasks(paired, unpaired, kmers, coverage, npar, mpi, cores, log_dir.path, tmpdir.path)

        if transabyss.debug:
            pdb.set_trace()

        if mpirun:
            jobs = []
            for task in TASKS:
                result = transabyss_run(task[0], task[1], task[2], task[3], task[4], task[5])
                jobs.append(result)
        else:
            processes_requested = set_threads("TransABySS", len(TASKS), transabyss.threads // cores, userdef=False)
            pool = create_pool(processes_requested)
            jobs = hard_worker(run_instance, TASKS, pool)
            close_pool(pool)

        if not jobs:
            logging.error("EXONtools ERROR: No valid results were produced by the ABYSS analysis.")
            raise EXONtoolsError

        logging.debug("All ABySS analyses have been completed: OK")

        if transabyss.debug:
            pdb.set_trace()

        # SAVE STATS
        if transabyss.stats and not transabyss.dryrun:
            transabyss_stats(jobs)

        # DELETE TEMPORARY DIRECTORY
        if not transabyss.keeptmp:
            tmpdir.delete()


def transabyss_pars(kmers, coverage, npar, mpirun, cores, parallel):
    """Checking transabyss input parameters"""

    logging.debug("Testing TransABySS input settings")

    if transabyss.extra:
        logging.warning("The following arguments will be added to 'TransABySS' command line:")
        logging.warning(transabyss.extra)

    if kmers:
        if kmers == "program default":
            logging.info("Using 25 kmer length (default)")
            kmers = [25]
        else:
            try:
                kmers = [int(x) for x in set(kmers)]
            except ValueError:
                logging.error("K-mer values must be integers")
                raise EXONtoolsError("K-mer format error")
            for k in kmers:
                if k < 10:
                    logging.error("All selected kmer lengths must be greater or equal 10 bp")
                    raise EXONtoolsError("Wrong format of k-mer length")
            kmers.sort()
            logging.info("TransABySS analysis will use the following k-mers: {0}".format(kmers))
        logging.debug("K-mer settings check: OK")
    else:
        logging.error("No k-mers were provided in -k option")
        raise EXONtoolsError("No k-mer options were provided")

    if coverage == 'program default':
        coverage = None
        logging.info("All contigs with coverage lower than the squared median will be dismissed")
    elif len(coverage) >= 1 and "".join(coverage).isdigit():
        try:
            coverage = [int(x) for x in set(coverage)]
        except ValueError:
            logging.error("Coverage values must be integers")
            raise EXONtoolsError("Coverage format error")
        for c in coverage:
            if c < 5:
                logging.error("Every coverage value provided in -c option must be positive integer greater or equal 5")
                raise EXONtoolsError("Wrong format for coverage option")
        logging.info("TransABySS analysis will use the following coverage thresholds: {0}".format(coverage))
        logging.debug("Coverage settings check: OK")
    else:
        logging.error("The unknown value type in coverage settings. Please check the option -c and try again")
        raise EXONtoolsError("Wrong format for coverage option")

    if not isinstance(npar, int) or npar < 5:
        logging.error("Number of required paired reads for contig assembly must be above or equal five (-n option)")
        raise EXONtoolsError("Number of paired reads per k-mer must be at least 5 or more")
    else:
        logging.debug("Minimum number of paired reads is set: OK")

    if mpirun:
        p = subprocess.Popen('which ' + mpirun, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        stdoutdata, _ = p.communicate()
        if stdoutdata and os.path.exists(stdoutdata.strip()):
            if parallel < 1:
                logging.error("The number of parallel machines cannot be less than one")
                raise EXONtoolsError("Wrong number of nodes")
            try:
                os.environ['PATH'] = os.path.dirname(stdoutdata) + ":" + os.environ['PATH']
            except OSError:
                logging.error("This is a potential bug. The mpi path cannot be exported")
                raise EXONtoolsError("Wrong path for mpirun")
            mpi = parallel
            logging.debug("MPI is set with {0:d} nodes: OK".format(parallel))
        else:
            logging.error("Provided path for mpirun does not exist. Please verify and try again")
            raise EXONtoolsError("Wrong path for mpirun")
    else:
        if parallel != 0:
            logging.warning("The mpirun path was not provided. The number of nodes is changed to 0.")
        if cores > transabyss.threads:
            logging.error("Single machine analysis uses -T/-threads general EXONtools option to set the number of parallel processes. Its argument cannot be larger than -j option (number of cores per job)")
            raise EXONtoolsError("Wrong -j option settings")
        logging.warning("Using a single machine for the ABySS analysis with the number of cores set by -T option")
        mpi = 0

    logging.debug("TransABySS program settings: OK")
    return mpi, kmers, coverage


def transabyss_tasks(paired, unpaired, kmers, coverage, npar, mpi, cores, logpath, tmppath):
    """Prepare TASK list for TransABySS analysis"""

    logging.debug("Preparing TransABySS input command arguments for the parallel analysis")

    librs = list(set([x.split("_")[0] for x in list(paired.keys()) + list(unpaired.keys())]))

    TASKS = []
    for lib in librs:

        if lib + "_paired" in paired:
            pefiles = " --pe {0:s}".format(" ".join(paired[lib + "_paired"]))
            sefiles = "--se {0:s}".format(" ".join(paired[lib + "_paired"]))
            if lib + "_unpaired" in unpaired:
                sefiles = sefiles + " {0:s}".format(unpaired[lib + "_unpaired"][0])
        else:
            if lib + "_unpaired" in unpaired:
                pefiles = ""
                sefiles = "--se {0:s}".format(unpaired[lib + "_unpaired"][0])
            else:
                logging.critical("Input files for '{}' library are not provided".format(lib))
        infiles = sefiles + pefiles

        libouttmp = os.path.join(tmppath, lib)
        makenewdir(libouttmp)

        if kmers and coverage:
            for k in kmers:
                for c in coverage:
                    kmername = "{0:s}_km{1:d}_cov{2:d}".format(lib, k, c)
                    kmeroutpath = os.path.join(libouttmp, kmername)
                    makenewdir(kmeroutpath)
                    pars = "--kmer {0:d} --cov {1:d} --eros {1:d} --pairs {2:d}".format(k, c, npar)
                    if mpi:
                        TASKS.append([kmeroutpath, logpath, infiles, pars, mpi, cores])
                    else:
                        TASKS.append(worker(transabyss_run, [kmeroutpath, logpath, infiles, pars, mpi, cores]))
        elif kmers and not coverage:
            for k in kmers:
                kmername = "{0:s}_km{1:d}".format(lib, int(k))
                kmeroutpath = os.path.join(libouttmp, kmername)
                makenewdir(kmeroutpath)
                pars = "--kmer {0:d} --pairs {1:d}".format(k, npar)
                if mpi:
                    TASKS.append([kmeroutpath, logpath, infiles, pars, mpi, cores])
                else:
                    TASKS.append(worker(transabyss_run, [kmeroutpath, logpath, infiles, pars, mpi, cores]))
        else:
            logging.error("Cannot create task because kmer length parameter is missing")
            raise EXONtoolsError("This bug occurred duting task creation")

    logging.debug("TransABySS input command arguments are created: OK")
    return TASKS


def transabyss_run(kmeroutpath, logpath, infiles, pars, mpi, cores):
    """Runs one instance of TransABySS analysis"""

    suffix = "-final.fa"

    outpath = os.path.dirname(logpath)
    kmername = os.path.basename(kmeroutpath)
    kmeroutlog = os.path.join(logpath, kmername + '.log')
    lib = kmername.split("_")[0]
    kmer = int(kmername.split("_")[1][2:])
    if len(kmername.split("_")) > 2:
        logging.info("Starting the assembly of '{0:s}' library with {1:d} k-mer length and {2:d} minimum coverage".format(lib, kmer, int(kmername.split("_")[2][3:])))
    else:
        logging.info("Starting the assembly of '{0:s}' library with {1:d} k-mer length".format(lib, kmer))
    makenewdir(kmeroutpath)
    run_executor(executor(
        program="TransABySS",
        params=[infiles, kmername, kmeroutpath, pars, cores, mpi],
        conditions={"positive": [cores], "pathexists": [kmeroutpath]},
        custom_arg_string=transabyss.extra + " > " + kmeroutlog + " 2>&1"
    ))
    if transabyss.dryrun:
        return "DRYRUN"
    else:
        if len(os.listdir(kmeroutpath)) > 0:
            outfile = os.path.join(kmeroutpath, kmername + suffix)
            assembly_file = os.path.join(outpath, kmername + transabyss.suffix + '.fa')
            shutil.copy2(outfile, assembly_file)
            shutil.rmtree(kmeroutpath)
            logging.debug("Finished assembling for '{0:s}' library. All resulted contigs were saved to {1:s}".format(lib, kmername + transabyss.suffix + ".fa"))
            return assembly_file
        else:
            logging.error("No files were found in the TransABySS output folder")
            raise EXONtoolsError("TransABySS output error")


def transabyss_stats(paths):
    """Write stats for ABySS analysis"""

    logging.debug("Collecting statistical data on assemblies")

    statpath = os.path.join(output.path, "STATS_transabyss_assembly.csv")

    logging.info("All assembly metrics will be saved to '{0:s}'".format(os.path.basename(statpath)))

    header = ["Library", "N", "N>500", "N>1000", "N>5000", "N>10000", "LongestID", "Min", "N75", "N50", "N25", "Max", "L50", "L75", "NG50", "Mean", "E-size", "Total", "GC", "N's", "N's per 100kb"]

    with open(statpath, 'w') as statout:
        csv_writer = csv.writer(statout)
        csv_writer.writerow(header)

        for path in sorted(paths, key=natural_sort):
            x = metrics(path)
            csv_writer.writerow([x.ID.replace(transabyss.suffix, ""), x.totnum, x.long500, x.long1000, x.long5000, x.long10000, x.longest, x.min, x.N75, x.N50, x.N25, x.max, x.L50, x.L75, x.NG, x.mean, x.esize, x.totlen, x.gc, x.nbases, x.n100])

    logging.debug("Stats data collected and written to a file: OK")
