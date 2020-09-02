# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import os
import logging
# import csv
import pdb
# import shutil

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import output, getinput, makenewdir, parseinput
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, run_executor, set_threads
from mains.EXT_errors import EXONtoolsError
# from mains.EXT_validator import positive
# from utils.sorting import natural_sort
from utils.seqIO import SeqIO


class cutadapt(EXTprogram):
    """
    Trimming and filtering FASTQ reads, including adapter removal
    """

    name = "cutadapt"

    def execute_program(self):
        args = self.args
        self.clean_reads(args.inpath, args.outdir, args.forward, args.reverse, args.adapters, args.select, args.minlen, args.minqual, args.leading, args.trailing, args.cut5end, args.cut3end, args.rqc, args.phred)

    def clean_reads(self, inpath, outdir, forward, reverse, adapters, select, minlen, minqual, leading, trailing, cut5end, cut3end, rqc, phred):

        if cutadapt.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        cutadapt.run_dry(getinput, output, makenewdir, executor)
        cutadapt.set_debug(executor)

        executor.setconfig("cutadapt")
        if rqc:
            executor.setconfig("fastqc", "multiqc")

        if cutadapt.debug:
            pdb.set_trace()

        # TESTING INPUT PARAMETERS
        if minlen < 0:
            logging.error("Minimal length (-m) must be greater or equal 0")
            raise EXONtoolsError("Wrong input parameter value")
        if minqual < 0:
            logging.error("Minimal quality (-q) must be greater or equal 0")
            raise EXONtoolsError("Wrong input parameter value")
        if leading < 0:
            logging.error("5'-end quality trimming (-q5) must be greater or equal 0")
            raise EXONtoolsError("Wrong input parameter value")
        if trailing < 0:
            logging.error("3'-end quality trimming (-q3) must be greater or equal 0")
            raise EXONtoolsError("Wrong input parameter value")
        if cut5end < 0:
            logging.error("5'-read end base trimming (-c5) must be greater or equal 0")
            raise EXONtoolsError("Wrong input parameter value")
        if cut3end < 0:
            logging.error("3'-read end base trimming (-c3) must be greater or equal 0")
            raise EXONtoolsError("Wrong input parameter value")
        logging.info("Input parameter check: OK")

        if cutadapt.debug:
            pdb.set_trace()

        # READING ADAPTERS
        if adapters:
            logging.info("Uploading adapter sequences")
            adaptfile = SeqIO(adapters, "FASTA")
            adaptfile.readall()
            logging.debug("Uploading adapter sequences: OK")

        if select:
            logging.warning("'ADAPTER SELECT MODE is ON': each library will be analyzed with correspoding adapter")

        if cutadapt.debug:
            pdb.set_trace()

        # GET READ INPUT FILES
        getinput.format(['.fq', '.fastq', '.fastq.gz', '.fq.gz'])
        paired, unpaired = parseinput(inpath, forward, reverse)

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE LOG DIR
        log_dir = makenewdir(name="dependency_logs", fullname="LOG")
        logging.info("Find all dependency program logs in the '{0:s}' folder".format(os.path.basename(log_dir.path)))

        # MAKE TMP DIRECTORY
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        logging.debug("IO settings: OK")

        if cutadapt.debug:
            pdb.set_trace()

        # COMBINING PARAMETERS
        add_params = ""
        if phred:
            add_params = "{0:s} --quality-base={1:s}".format(add_params, phred[1:])
        if leading and not trailing:
            add_params = "{0:s} -q {1:d},0".format(add_params, leading)
        elif trailing and not leading:
            add_params = "{0:s} -q {1:d}".format(add_params, trailing)
        elif trailing and leading:
            add_params = "{0:s} -q {1:d},{2:d}".format(add_params, leading, trailing)
        if cut5end:
            add_params = "{0:s} -u {1:d}".format(add_params, cut5end)
        if cut3end:
            add_params = "{0:s} -u -{1:d}".format(add_params, cut3end)
        if minlen:
            add_params = "{0:s} -m {1:d}".format(add_params, minlen)

        TASKS = []

        for sample in unpaired.keys():

            final_params = "" + add_params

            if adapters and select:
                for inseq in adaptfile.SEQS:
                    if sample.replace("unpaired", "") in adaptfile.SEQS[inseq].name:
                        final_params = "{0:s} -b {1:s}".format(final_params, adaptfile.SEQS[inseq][0].seq)
            if adapters and not select:
                for inseq in adaptfile.SEQS:
                    final_params = "{0:s} -b {1:s}".format(final_params, adaptfile.SEQS[inseq][0].seq)

        outmask = os.path.splitext(os.path.basename(unpaired[sample][0]))[0]
        log_file = os.path.join(log_dir.path, sample + "_clean.log")
        outpath = os.path.join(tmpdir.path, outmask + "_cclean.fq")
        filepars = "{0:s} > {1:s} 2> ".format(unpaired[sample][0], outpath)
        # print(outmask)
        # print(log_file)
        # print(outpath)
        # print(filepars)
        TASKS.append(executor(
            program="cutadapt",
            params=[final_params, filepars, log_file],
            conditions={"positive": [cutadapt.threads]},
            custom_arg_string=cutadapt.extra)
        )

        logging.info("Running the FASTQ read checkup analysis for all files")
        if TASKS:
            processes_requested = set_threads("Read cleaner (cutadapt)", len(TASKS), cutadapt.threads)
            pool = create_pool(processes_requested)
            jobs = hard_worker(run_executor, TASKS, pool)
            close_pool(pool)

        print(jobs)

        # for sample in paired.keys():
        #     # CHOOSING SELECTED ADAPTERS
        #     final_params = "" + add_params
        #     if adapters and select:
        #         for inseq in adaptfile.SEQS:
        #             if sample.replace("paired", "") in inseq.name:
        #                 final_params = "{0:s} -b {1:s} -B {2:s}".format(final_params, inseq.seq, inseq.seq)
        #     if adapters and not select:
        #         for adptr in adaptfile.SEQS:
        #             final_params = "{0:s} -b {1:s} -B {2:s}".format(final_params, inseq.seq, inseq.seq)

        #     if cut5end:
        #         final_params = "{0:s} -U {1:d}".format(final_params, cut5end)

        #     if cut3end:
        #         final_params = final_params + "-U -%d " % cut3end

        #     outmask1 = os.path.splitext(os.path.basename(paired[sample][0]))[0]
        #     outmask2 = os.path.splitext(os.path.basename(paired[sample][1]))[0]
        #     log_file = os.path.join(log_path, sample + "_clean.log")
        #     output1 = os.path.join(tmpdir, outmask1 + "_cclean_paired.fq")
        #     output2 = os.path.join(tmpdir, outmask2 + "_cclean_paired.fq")

        #     filepars = "--output=%s --paired-output=%s %s %s > " % (output1, output2, paired[sample][0], paired[sample][1])

        #     TASKS.append(executor(
        #         program="cutadapt",
        #         params=[final_params, filepars, log_file],
        #         conditions={"positive": [threads]},
        #         custom_arg_string=extra,
        #         doprint=False
        #     )
        #     )
