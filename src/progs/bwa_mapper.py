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
import shutil

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import parseinput, getinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, run_executor, set_threads
from mains.EXT_errors import EXONtoolsError
from mains.EXT_validator import positive
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from progs.asssam import analyze_sam


class bwa_mapper(EXTprogram):
    """This program maps reads to a reference genome or to assembly"""

    name = "bwa"

    def execute_program(self):
        args = self.args
        self.process_contigs(args.inpath, args.forward, args.reverse, args.unpaired, args.reference, args.outdir, args.mismatch, args.gapopen, args.gapextension, args.clipping, args.discordant)

    def process_contigs(self, inpath, forward, reverse, unpaired_in, reference, outdir, mismatch, gapopen, gapextension, clipping, discordant):

        if bwa_mapper.debug:
            pdb.set_trace()

        executor.setconfig("samtools", "bwa")

        # testing parameters
        positive([mismatch, gapopen, gapextension, clipping, discordant], strict=False)

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        bwa_mapper.run_dry(SeqIO, getinput, output, makenewdir, executor)
        bwa_mapper.set_debug(executor)

        # GET READ INPUT FILES
        getinput.format(['.fq', '.fastq'])
        paired, unpaired = parseinput(inpath, forward, reverse)
        if unpaired_in and not inpath:
            unpaired = {os.path.basename(unpaired_in).split("_")[0].split(".")[0] + "_unpaired": getinput(unpaired_in).files}
        elif unpaired_in and inpath:
            logging.error("Input arguments '-i' and '-U' cannot be used simultaneously")
            raise EXONtoolsError("Input path error")
        else:
            pass

        # GET REFERENCE FASTA FILES
        getinput.format(['.fa', '.fasta'])
        FileList = getinput(reference).files
        ref_files = {}
        for x in FileList:
            lib = os.path.basename(x).split("_")[0].split(".")[0]
            if lib + "_paired" in paired or lib + "_unpaired" in unpaired:
                logging.info("'{0:s}' library has matched sequencing read and assembly files".format(lib))
                ref_files[lib] = x
        if not ref_files:
            logging.error("No reference files matched with any read files")
            raise EXONtoolsError("bwa_mapper input error")
        del FileList

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE TMP DIRECTORIES
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        # MAKE LOG DIR
        log_dir = makenewdir(name="dependency_logs", fullname="LOG")
        logging.info("Find all dependency program logs in '{0:s}' folder".format(os.path.basename(log_dir.path)))

        logging.debug("IO settings: OK")

        logging.info("***********************************************************")
        logging.info("Starting READ MAPPING procedure for each reference library")

        if bwa_mapper.debug:
            pdb.set_trace()

        if bwa_mapper.extra:
            logging.warning("The following extra arguments will be added to 'bwa' command line:")
            logging.warning(bwa_mapper.extra)

        bam_final = {}
        for lib in sorted(ref_files.keys(), key=natural_sort):

            liblog = os.path.join(log_dir.path, lib + "_deps.log")
            libdir = makenewdir(os.path.join(tmpdir.path, lib))
            libdbref = makenewdir(os.path.join(libdir.path, lib + "_dbref"))
            dbprefix = os.path.join(libdbref.path, lib)

            logging.info("***********************************************************")
            logging.info("'bwa' program is building the database for '{0:s}' library".format(lib))
            run_executor(executor(
                program="bwa",
                params=['index', " -p " + dbprefix, ref_files[lib]],
                conditions={"pathexists": [ref_files[lib]]},
                custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                quiet=True)
            )

            if not bwa_mapper.dryrun:
                shutil.copy2(ref_files[lib], libdbref.path)
            refpath = os.path.join(libdbref.path, os.path.basename(ref_files[lib]))

            logging.info("'samtools' program is constructing the dictionary for '{0:s}' library".format(lib))
            run_executor(executor(
                program="SAMTOOLS",
                params=["faidx", refpath],
                conditions={},
                custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                quiet=True)
            )

            bam_files = []

            # PAIRED READS
            try:
                bwa_files = "{0:s} {1:s} {2:s}".format(dbprefix, paired[lib + "_paired"][0], paired[lib + "_paired"][1])
                bwa_outpath = os.path.join(libdir.path, lib + "_paired.sam")
                bam_paired = os.path.splitext(bwa_outpath)[0] + '.bam'
                fixed_path = os.path.splitext(bam_paired)[0] + "_fixed.bam"
                paired_sort_path = os.path.join(libdir.path, lib + "_sort_paired.bam")
                bwapars = "{6:s} -t {5:d} -B {0:d} -O {1:d} -E {2:d} -L {3:d} -U {4:d} -Y".format(mismatch, gapopen, gapextension, clipping, discordant, bwa_mapper.threads, bwa_mapper.extra)

                logging.info("'bwa' program is mapping paired reads for '{0:s}' library".format(lib))
                run_executor(executor(
                    program="bwa",
                    params=["mem", bwapars, bwa_files],
                    conditions={"positive": [bwa_mapper.threads]},
                    custom_arg_string='> {0:s} 2>>{1:s}'.format(bwa_outpath, liblog),
                    quiet=True)
                )

                logging.info("'samtools' program is converting SAM to BAM format for '{0:s}' library".format(lib))
                samtoolspars = "--threads {0:d} -b -o {1:s} {2:s}".format(bwa_mapper.threads, bam_paired, bwa_outpath)
                run_executor(executor(
                    program="SAMTOOLS",
                    params=["view", samtoolspars],
                    conditions={"pathexists": [bwa_outpath], "positive": [bwa_mapper.threads]},
                    custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                    quiet=True)
                )

                if not bwa_mapper.dryrun:
                    os.remove(bwa_outpath)

                logging.info("'samtools' program is fixing paired reads in '{0:s}' library".format(lib))
                samtoolspars = "--threads {0:d} -m -O bam {1:s} {2:s}".format(bwa_mapper.threads, bam_paired, fixed_path)
                run_executor(executor(
                    program="SAMTOOLS",
                    params=["fixmate", samtoolspars],
                    conditions={"pathexists": [bam_paired], "positive": [bwa_mapper.threads]},
                    custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                    quiet=True)
                )

                if not bwa_mapper.dryrun:
                    os.remove(bam_paired)

                logging.info("'samtools' program is sorting the bam file with paired reads from '{0:s}' library".format(lib))
                samtoolspars = "--threads {0:d} -O bam -o {1:s} {2:s}".format(bwa_mapper.threads, paired_sort_path, fixed_path)
                run_executor(executor(
                    program="SAMTOOLS",
                    params=["sort", samtoolspars],
                    conditions={"pathexists": [fixed_path], "positive": [bwa_mapper.threads]},
                    custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                    quiet=True)
                )

                if not bwa_mapper.dryrun:
                    os.remove(fixed_path)

                bam_files.append(paired_sort_path)

            except KeyError:                        # skips the whole part if there are no paired reads
                pass

            # UNPAIRED READS
            try:
                bwa_files = "{0:s} {1:s}".format(dbprefix, unpaired[lib + "_unpaired"][0])
                bwa_outpath = os.path.join(libdir.path, lib + "_unpaired.sam")
                bam_unpaired = os.path.splitext(bwa_outpath)[0] + '.bam'
                fixed_path = os.path.splitext(bam_unpaired)[0] + "_fixed.bam"
                unpaired_sort_path = os.path.join(libdir.path, lib + "_sort_unpaired.bam")
                bwapars = "{6:s} -t {5:d} -B {0:d} -O {1:d} -E {2:d} -L {3:d} -U {4:d}".format(mismatch, gapopen, gapextension, clipping, discordant, bwa_mapper.threads, bwa_mapper.extra)

                logging.info("'bwa' program is mapping unpaired reads for '{0:s}' library".format(lib))
                run_executor(executor(
                    program="bwa",
                    params=["mem", bwapars, bwa_files],
                    conditions={"positive": [bwa_mapper.threads]},
                    custom_arg_string='> {0:s} 2>>{1:s}'.format(bwa_outpath, liblog),
                    quiet=True)
                )

                logging.info("'samtools' program is converting SAM to BAM format for '{0:s}' library".format(lib))
                samtoolspars = "--threads {0:d} -b -o {1:s} {2:s}".format(bwa_mapper.threads, bam_unpaired, bwa_outpath)
                run_executor(executor(
                    program="SAMTOOLS",
                    params=["view", samtoolspars],
                    conditions={"pathexists": [bwa_outpath], "positive": [bwa_mapper.threads]},
                    custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                    quiet=True)
                )

                if not bwa_mapper.dryrun:
                    os.remove(bwa_outpath)

                logging.info("'samtools' program is fixing unpaired reads in '{0:s}' library".format(lib))
                samtoolspars = "--threads {0:d} -m -O bam {1:s} {2:s}".format(bwa_mapper.threads, bam_unpaired, fixed_path)
                run_executor(executor(
                    program="SAMTOOLS",
                    params=["fixmate", samtoolspars],
                    conditions={"pathexists": [bam_unpaired], "positive": [bwa_mapper.threads]},
                    custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                    quiet=True)
                )

                if not bwa_mapper.dryrun:
                    os.remove(bam_unpaired)

                logging.info("'samtools' program is sorting the bam file with unpaired reads from '{0:s}' library".format(lib))
                samtoolspars = "--threads {0:d} -O bam -o {1:s} {2:s}".format(bwa_mapper.threads, unpaired_sort_path, fixed_path)
                run_executor(executor(
                    program="SAMTOOLS",
                    params=["sort", samtoolspars],
                    conditions={"pathexists": [fixed_path], "positive": [bwa_mapper.threads]},
                    custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                    quiet=True)
                )

                if not bwa_mapper.dryrun:
                    os.remove(fixed_path)

                bam_files.append(unpaired_sort_path)

            except KeyError:                            # skip the whole part if there are no unpaired files
                pass

            # MERGED FILES
            merged_path = os.path.join(libdir.path, lib + "_merged.bam")
            bampath = os.path.join(output.path, lib + bwa_mapper.suffix + ".bam")
            if not bam_files:
                logging.error("No bams files are produced after mapping")
                raise EXONtoolsError("bwa error")
            elif len(bam_files) == 1:
                if not bwa_mapper.dryrun:
                    shutil.move(bam_files[0], merged_path)
            else:
                logging.info("'samtools' program is merging bam files for '{0:s}' library".format(lib))
                samtoolspars = "--threads {0:d} -O bam -cp {1:s} {2:s}".format(bwa_mapper.threads, merged_path, " ".join(bam_files))
                run_executor(executor(
                    program="SAMTOOLS",
                    params=["merge", samtoolspars],
                    conditions={'pathexists': bam_files, 'positive': [bwa_mapper.threads]},
                    custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                    quiet=True)
                )

                if not bwa_mapper.dryrun:
                    [os.remove(x) for x in bam_files]

            logging.info("'samtools' program is marking duplicated reads in '{0:s}' library".format(lib))
            samtoolspars = "--threads {0:d} -S -O bam {1:s} {2:s}".format(bwa_mapper.threads, merged_path, bampath)
            run_executor(executor(
                program="SAMTOOLS",
                params=["markdup", samtoolspars],
                conditions={'pathexists': [merged_path], 'positive': [bwa_mapper.threads]},
                custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                quiet=True)
            )

            if not bwa_mapper.dryrun:
                os.remove(merged_path)

            bam_final.update({lib: (bampath, refpath)})
            logging.info("'samtools' program is indexing the bam file for '{0:s}' library".format(lib))
            run_executor(executor(
                program="SAMTOOLS",
                params=["index", "-b " + bampath],
                conditions={'pathexists': [bampath]},
                custom_arg_string=' >> {0:s} 2>&1'.format(liblog),
                quiet=True)
            )
            libdir.delete()

        if bwa_mapper.debug:
            pdb.set_trace()

        if bwa_mapper.stats:
            statpath = os.path.join(output.path, "mapping_stats.csv")
            logging.info("Collecting stats for read mapping")
            logging.info("All read mapping metrics will be saved to {0:s}".format(os.path.basename(statpath)))
            logging.debug("Running multiprocessing analysis")
            TASKS = [worker(analyze_sam, [bam_final[x][0], tmpdir.path, 0, 0]) for x in bam_final]
            processes_requested = set_threads("asssam", len(TASKS), bwa_mapper.threads, userdef=False)
            pool = create_pool(processes_requested)
            stat_collector = hard_worker(run_instance, TASKS, pool)
            close_pool(pool)

            logging.debug("All BAM files were successfully analyzed: OK")

            header = ["Library", "Total_reads", "Total_mapped", "Total_paired", "All_R1_map", "All_R2_map", "Concord_map", "Concord_map_unique", "Discord_map", "Discord_map_unique", "R1_single", "R1_single_unique", "R2_single", "R2_single_unique", "Total_unpaired", "Unpaired_map", "Unpaired_map_unique", "Secondary", "Supplementary", "Duplicates", "Filtered", "Insert_length", "Coverage", "", "%_total_map", "%_all_R1_map", "%_all_R2_map", "%_concord_map", "%_concord_map_unique", "%_discord_map", "%_discord_map_unique", "%_R1_single", "%_R1_single_unique", "%_R2_single", "%_R2_single_unique", "%_unpaired_map", "%_unpaired_map_unique", "%_secondary", "%_duplicates"]

            logging.debug("Ready to write stats to the output file")

            if not bwa_mapper.dryrun:
                with open(statpath, 'w') as statout:
                    csv_writer = csv.writer(statout)
                    csv_writer.writerow(header)

                    for x in sorted(stat_collector, key=lambda x: natural_sort(x["library"])):
                        csv_writer.writerow([
                            x["library"],
                            x["total"],
                            x["total_mapped"],
                            x["paired"],
                            x["map_all_r1"],
                            x["map_all_r2"],
                            x["map_concord"],
                            x["map_concord_unique"],
                            x["map_discord"],
                            x["map_discord_unique"],
                            x["map_only_r1"],
                            x["map_only_r1_unique"],
                            x["map_only_r2"],
                            x["map_only_r2_unique"],
                            x["unpaired"],
                            x["map_unpaired"],
                            x["map_unpaired_unique"],
                            x["secondary"],
                            x["supplementary"],
                            x["duplicates"],
                            x["filtered_reads"],
                            x["averinsert"],
                            x["avercoverage"],
                            "",
                            x["%_total_map"],
                            x["%_all_R1_map"],
                            x["%_all_R2_map"],
                            x["%_concord_map"],
                            x["%_concord_map_unique"],
                            x["%_discord_map"],
                            x["%_discord_map_unique"],
                            x["%_R1_single"],
                            x["%_R1_single_unique"],
                            x["%_R2_single"],
                            x["%_R2_single_unique"],
                            x["%_unpaired_map"],
                            x["%_unpaired_map_unique"],
                            x["%_secondary"],
                            x["%_duplicates"]
                        ])

            logging.info("All read mapping metrics were successfully collected")

        if not bwa_mapper.keeptmp:
            tmpdir.delete()
