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
    """This program finds and fixes sequencing errors within fastq files"""

    name = "hammer"

    def execute_program(self):
        args = self.args
        self.correct_errors(args.inpath, args.forward, args.reverse, args.unpaired, args.outdir, args.gzoutput, args.torrent, args.rqc)

    def correct_errors(self, inpath, forward, reverse, unpaired_in, outdir, gzoutput, torrent, rqc):

        if hammer.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        hammer.run_dry(SeqIO, getinput, output, makenewdir, executor)
        hammer.set_debug(executor)

        if rqc:
            executor.setconfig("spades","fastqc", "multiqc")

        if hammer.debug:
            pdb.set_trace()

        logging.debug("Print warnings")
        if torrent:
            logging.warning("Ion Torrent data type is set for SPAdes")

        if gzoutput:
            logging.warning("All output FASTQ files will be compressed with gzip")

        if hammer.debug:
            pdb.set_trace()

        # IO SETTINGS 
        getinput.format(['.fq', '.fastq', '.fastq.gz', '.fq.gz'])
        paired, unpaired = parseinput(inpath, forward, reverse)

        # GET INPUT FILES
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

        # MAKE TMP DIRECTORY
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        # MAKE LOG DIRECTORY
        log_dir = makenewdir(name="spades_logs", fullname="LOG")
        logging.info("All SPAdes log files will be saved in '{0:s} directory'".format(os.path.basename(log_dir.path)))

        logging.debug("Test and IO settings: OK")

        if hammer.debug:
            pdb.set_trace()

        # RUNING BAYESHAMMER
        if torrent:
            logging.info("Correcting sequencing errors with SPAdes IonHammer")
        else:
            logging.info("Correcting sequencing errors with SPAdes BayesHammer")

        ###UNPAIRED READS###
        logging.debug("Setting up SPAdes parameters")
        SPADES_params = '--only-error-correction'
        if torrent:
            SPADES_params = SPADES_params + " --iontorrent"
        if not gzoutput:
            SPADES_params = SPADES_params + " --disable-gzip-output"
        if hammer.extra:
            logging.warning("The following arguments will be added to 'SPAdes' command line:")
            logging.warning(spades.extra)

        logging.debug("SPAdes parameters: OK")

        if hammer.debug:
            pdb.set_trace()
 
        logging.debug("Starting SPAdes analysis")
        for lib in librs:
            logging.info("Correcting '{0:s}' library".format(lib))
            if lib + "_paired" in paired and lib + "_unpaired" in unpaired:
                infiles = "-1 {0:s} -2 {1:s} -s {2:s}".format(paired[lib + "_paired"][0], paired[lib + "_paired"][1], unpaired[lib + "_unpaired"][0])
            elif lib + "_paired" in paired and not lib + "_unpaired" in unpaired:
                infiles = "-1 {0:s} -2 {1:s}".format(paired[lib + "_paired"][0], paired[lib + "_paired"][1])
            elif not lib + "_paired" in paired and lib + "_unpaired" in unpaired:
                infiles = "-s {0:s}".format(unpaired[lib + "_unpaired"][0])
            else:
                logging.error("You found the bug during parsing of SPAdes input files")
                raise EXONtoolsError("Input files command line error")
            
            libouttmp = os.path.join(tmpdir.path, lib + "_hammer")
            makenewdir(libouttmp)
            logpath = os.path.join(log_dir.path, lib + "_spades.log")


            run_executor(executor(
                program="SPAdes",
                params=[SPADES_params, libouttmp, hammer.threads, int(hammer.memory), infiles],
                conditions={"positive": [hammer.threads], "pathexists": [libouttmp]},
                custom_arg_string=hammer.extra + " > /dev/null",
            ))

            if not hammer.dryrun:
                spadelogpath = os.path.join(libouttmp, "spades.log")
                shutil.move(spadelogpath, logpath)
        
            logging.debug("SPAdes read cprrection for the library '{0:s}' is finished: OK".format(lib))

        logging.debug("Successfully finished all corrections: OK")

        if hammer.debug:
            pdb.set_trace()

        # stats_collector = {}
        # TASKS = []
  


def copy_log(from_path, to_path, sample):
    """copies SPADES log file"""
    try:
        shutil.copy2(os.path.join(from_path,"spades.log"),os.path.join(to_path,sample+".log"))
    except IOError:
        # This error should not occur
        logging.error("I can't find spades.log file for '{}' sample".format(sample))
    except Exception as e:
        logging.error(e) 
        raise EXONtoolsError