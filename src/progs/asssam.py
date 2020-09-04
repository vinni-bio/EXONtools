# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.


from __future__ import print_function, division
from collections import Counter
import os
import logging
import csv
import pdb

from mains.EXT_prog import EXTprogram
from mains.EXT_errors import EXONtoolsError
from mains.EXT_IO import getinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, run_executor, set_threads
from mains.EXT_validator import positive
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from utils.phred import phred, phredmode


class asssam(EXTprogram):
    """Estimate read mapping metrics for all sam/bam files within the provided path"""

    name = "asssam"

    def execute_program(self):
        args = self.args
        self.evaluate(args.inpath, args.outdir, args.minlen, args.minqual)

    def evaluate(self, inpath, outdir, minlen, minqual):

        if asssam.debug:
            pdb.set_trace()

        asssam.run_dry(getinput, output, makenewdir, executor, SeqIO)

        positive([minlen], False)

        logging.info("Collecting stats for read mapping success in {0:s}".format(inpath))

        # GET INPUT FILES
        getinput.format(['.sam', '.bam'])
        FileList = getinput(inpath).files

        # GET OUTPUT PATH
        output.nostrict()
        output(outdir)
        outpath = os.path.join(output.path, "mapping_stats.csv")

        # MAKE TMP DIRECTORIES
        check = False
        for x in FileList:
            if x.endswith(".bam"):
                check = True
                executor.setconfig("samtools")
                break
        if check and not asssam.dryrun:
            tmpdir = makenewdir(name="tmp", fullname="temporary")
            tmppath = tmpdir.path
        else:
            tmppath = ""

        logging.info("Read mapping metrics will be saved to:")
        logging.info(outpath)

        if os.path.exists(outpath):
            logging.warning("The file 'mapping_stats.csv' exists and will be overwritten")

        logging.debug("Input files and checks: OK")

        if asssam.debug:
            pdb.set_trace()

        # COLLECT STATS
        logging.info("Starting to collect read mapping stats")

        logging.debug("Running multiprocessing analysis")
        TASKS = [worker(analyze_sam, [x, tmppath, minlen, minqual]) for x in FileList]
        processes_requested = set_threads("asssam", len(TASKS), asssam.threads, userdef=False)
        pool = create_pool(processes_requested)
        stat_collector = hard_worker(run_instance, TASKS, pool)
        close_pool(pool)

        logging.debug("All SAM files were successfully analyzed: OK")

        header = ["Library", "Total_reads", "Total_mapped", "Total_paired", "All_R1_map", "All_R2_map", "Concord_map", "Concord_map_unique", "Discord_map", "Discord_map_unique", "R1_single", "R1_single_unique", "R2_single", "R2_single_unique", "Total_unpaired", "Unpaired_map", "Unpaired_map_unique", "Secondary", "Supplementary", "Duplicates", "Filtered", "Insert_length", "Coverage", "", "%_total_map", "%_all_R1_map", "%_all_R2_map", "%_concord_map", "%_concord_map_unique", "%_discord_map", "%_discord_map_unique", "%_R1_single", "%_R1_single_unique", "%_R2_single", "%_R2_single_unique", "%_unpaired_map", "%_unpaired_map_unique", "%_secondary", "%_duplicates"]

        logging.debug("Ready to write stats to the output file")

        if asssam.debug:
            pdb.set_trace()

        if not asssam.dryrun:
            with open(outpath, 'w') as statout:
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

        if check and not asssam.dryrun:
            tmpdir.delete()


def analyze_sam(inpath, tmppath, minlen, minqual):
    """This function parses sam file and estimates mapping metrics"""

    libname = os.path.basename(inpath).split("_")[0].split(".")[0]
    logging.info("Collecting read mapping stats for '{0:s}' library".format(libname))

    samstats = {}
    filtered = 0
    libflags = []
    coverage = {}
    inserts = []
    supplementary = []
    discordsec = []
    concordsec = []
    forwardsec = []
    reversesec = []
    unpairedsec = []

    if inpath.endswith(".bam"):
        sampath = os.path.join(tmppath, os.path.basename(os.path.splitext(inpath)[0] + ".sam"))
        samviewpars = "--threads 1 -O sam -o {0:s} {1:s}".format(sampath, inpath)
        run_executor(executor(
            program="SAMTOOLS",
            params=["view", samviewpars],
            conditions={"pathexists": [inpath]},
            custom_arg_string=" > /dev/null 2>&1",
            quiet=True)
        )
    else:
        sampath = inpath

    pmode = phredmode(sampath, "SAM")

    samfile = SeqIO(sampath, "SAM")

    count = 0
    for sam in samfile.read():
        if count > 100000:
            unpairedsec = list(set(unpairedsec))
            concordsec = list(set(concordsec))
            discordsec = list(set(discordsec))
            forwardsec = list(set(forwardsec))
            reversesec = list(set(reversesec))
            supplementary = list(set(supplementary))
            count = 0
        if sam.readlen >= minlen and phred(sam.qual, mode=pmode).average() >= minqual:
            if sam.secondary or sam.secs:
                count += 1
                if sam.type == 0 and not sam.flag & 0x4:
                    unpairedsec.append(sam.name)
                elif sam.flag & 0x2:
                    concordsec.append(sam.name)
                elif not sam.flag & 0x4 and not sam.flag & 0x8:
                    discordsec.append(sam.name)
                elif sam.type == 1 and not sam.flag & 0x4:
                    forwardsec.append(sam.name)
                elif sam.type == 2 and not sam.flag & 0x4:
                    reversesec.append(sam.name)
                else:
                    logging.warning("Unknown secondary read flag '{0:d}' in read '{1:0}'".format(sam.flag, sam.name))
            if sam.supplementary or sam.suppls:
                count += 1
                supplementary.append(sam.name)
            if not sam.supplementary and not sam.secondary:
                libflags.append(sam.flag)
                if sam.insert:
                    inserts.append(sam.insert)
                if sam.target != "*" and sam.start != -1:
                    try:
                        coverage[sam.target] += 1
                    except KeyError:
                        coverage[sam.target] = 1
        else:
            if not sam.supplementary and not sam.secondary:
                filtered += 1

    if inpath.endswith(".bam"):
        samfile.delete()

    discordsec = len(set(discordsec))
    concordsec = len(set(concordsec))
    forwardsec = len(set(forwardsec))
    reversesec = len(set(reversesec))
    unpairedsec = len(set(unpairedsec))
    supplementary = len(set(supplementary))
    secondary = discordsec + concordsec + forwardsec + reversesec + unpairedsec

    if not asssam.dryrun and coverage:
        samstats = mapping_stats(libflags)
        samstats["map_concord_unique"] = samstats["map_concord"] - concordsec
        samstats["map_discord_unique"] = samstats["map_discord"] - discordsec
        samstats["library"] = libname
        samstats["secondary"] = secondary
        samstats["supplementary"] = supplementary
        samstats["filtered_reads"] = filtered
        samstats["map_only_r1_unique"] = samstats["map_only_r1"] - forwardsec
        samstats["map_only_r2_unique"] = samstats["map_only_r2"] - reversesec
        samstats["map_unpaired_unique"] = samstats["map_unpaired"] - unpairedsec
        try:
            samstats["averinsert"] = int(round(sum(inserts) / len(inserts), 0))
        except ZeroDivisionError:
            samstats["averinsert"] = 0
        try:
            samstats["%_concord_map_unique"] = round(samstats["map_concord_unique"] / samstats["paired"] * 100, 2)
        except ZeroDivisionError:
            samstats["%_concord_map_unique"] = 0
        try:
            samstats["%_discord_map_unique"] = round(samstats["map_discord_unique"] / samstats["paired"] * 100, 2)
        except ZeroDivisionError:
            samstats["%_discord_map_unique"] = 0
        try:
            samstats["%_unpaired_map_unique"] = round(samstats["map_unpaired_unique"] / samstats["unpaired"] * 100, 2)
        except ZeroDivisionError:
            samstats["%_unpaired_map_unique"] = 0
        try:
            samstats["%_R1_single_unique"] = round(samstats["map_only_r1_unique"] / samstats["paired"] * 100, 2)
        except ZeroDivisionError:
            samstats["%_R1_single_unique"] = 0
        try:
            samstats["%_R2_single_unique"] = round(samstats["map_only_r2_unique"] / samstats["paired"] * 100, 2)
        except ZeroDivisionError:
            samstats["%_R2_single_unique"] = 0
        try:
            samstats["%_secondary"] = round(samstats["secondary"] / samstats["total_mapped"] * 100, 2)
        except ZeroDivisionError:
            samstats["%_secondary"] = 0
        samstats["avercoverage"] = int(round(sum(coverage.values()) / len(coverage), 0))

    logging.debug("Mapping stats for '{0:s}' library: OK".format(libname))
    return samstats


def mapping_stats(flaglist):
    """Collects mapping stats from flag ints"""

    flagids = {
        0x1: "paired",
        0x2: "proper",
        0x4: "unmapped",
        0x8: "unmappedpair",
        0x10: "revcompl",
        0x20: "revcomplpair",
        0x40: "forward",
        0x80: "reverse",
        0x100: "secondary",
        0x200: "filtered",
        0x400: "duplicate",
        0x800: "supplementary"
    }

    statresults = {
        "total": 0,
        "paired": 0,
        "unpaired": 0,
        "fail_concord": 0,
        "map_concord": 0,
        "map_discord": 0,
        "map_all_r1": 0,
        "map_only_r1": 0,
        "map_all_r2": 0,
        "map_only_r2": 0,
        "fail_unpaired": 0,
        "map_unpaired": 0,
        "duplicates": 0,
        "total_mapped": 0
    }

    flagdict = Counter(flaglist)

    for samflag in flagdict:
        result = []
        for flag in flagids:
            if samflag & flag:
                result.append(flagids[flag])
        flagset = set(result)

        if len(set(["supplementary", "secondary"]) - flagset) != 2:
            continue

        statresults["total"] += flagdict[samflag]
        if "duplicate" in flagset:
            statresults["duplicates"] += flagdict[samflag]
        if "paired" in flagset:
            statresults["paired"] += flagdict[samflag]
            if "proper" in flagset:
                statresults["total_mapped"] += flagdict[samflag]
                statresults["map_concord"] += flagdict[samflag]
            else:
                statresults["fail_concord"] += flagdict[samflag]
                if all(elem not in flagset for elem in ["unmapped", "unmappedpair"]):
                    statresults["total_mapped"] += flagdict[samflag]
                    statresults["map_discord"] += flagdict[samflag]
                else:
                    if "unmapped" in flagset:
                        pass
                    elif "forward" in flagset:
                        statresults["total_mapped"] += flagdict[samflag]
                        statresults["map_only_r1"] += flagdict[samflag]
                    elif "reverse" in flagset:
                        statresults["total_mapped"] += flagdict[samflag]
                        statresults["map_only_r2"] += flagdict[samflag]
                    else:
                        logging.error("You found the bug in 'mapping_stats' function")
                        raise EXONtoolsError("Oops.. This error should not happen. See log message")
        else:
            statresults["unpaired"] += flagdict[samflag]
            if "unmapped" not in flagset:
                statresults["total_mapped"] += flagdict[samflag]
                statresults["map_unpaired"] += flagdict[samflag]

    statresults["paired"] = statresults["paired"] // 2
    statresults["map_concord"] = statresults["map_concord"] // 2
    statresults["map_discord"] = statresults["map_discord"] // 2
    statresults["fail_concord"] = statresults["fail_concord"] // 2
    statresults["map_all_r1"] = statresults["map_only_r1"] + statresults["map_concord"] + statresults["map_discord"]
    statresults["map_all_r2"] = statresults["map_only_r2"] + statresults["map_concord"] + statresults["map_discord"]
    statresults["fail_unpaired"] = statresults["unpaired"] - statresults["map_unpaired"]

    try:
        statresults["%_total_map"] = round(statresults["total_mapped"] / statresults["total"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_total_map"] = 0
    try:
        statresults["%_all_R1_map"] = round(statresults["map_all_r1"] / statresults["paired"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_all_R1_map"] = 0
    try:
        statresults["%_all_R2_map"] = round(statresults["map_all_r2"] / statresults["paired"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_all_R2_map"] = 0
    try:
        statresults["%_concord_map"] = round(statresults["map_concord"] / statresults["paired"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_concord_map"] = 0
    try:
        statresults["%_discord_map"] = round(statresults["map_discord"] / statresults["paired"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_discord_map"] = 0
    try:
        statresults["%_R1_single"] = round(statresults["map_only_r1"] / statresults["paired"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_R1_single"] = 0
    try:
        statresults["%_R2_single"] = round(statresults["map_only_r2"] / statresults["paired"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_R2_single"] = 0
    try:
        statresults["%_unpaired_map"] = round(statresults["map_unpaired"] / statresults["unpaired"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_unpaired_map"] = 0
    try:
        statresults["%_duplicates"] = round(statresults["duplicates"] / statresults["total_mapped"] * 100, 2)
    except ZeroDivisionError:
        statresults["%_duplicates"] = 0

    return statresults
