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
from math import log as LN
from functools import reduce as reduce2
from random import choice as choose

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, run_executor, set_threads
from mains.EXT_errors import EXONtoolsError
from utils.sorting import natural_sort, mostcommon
from utils.seqIO import SeqIO
from utils.phred import phredmode
from utils.gccheck import GC, Ncheck


class basecaller(EXTprogram):
    """This program analyzes BAM/SAM mapping results, estimates
    the coverage of each contig and identifies its heterozygous sites.
    It can also trim contig ends by required threshold of read coverage.
    Finally, it produces fasta file with called bases."""

    name = "basecaller"

    def execute_program(self):
        args = self.args
        self.process_contigs(args.inpath, args.reference, args.outdir, args.coverage, args.flanking, args.maxerror, args.minor, args.minmapqual, args.minreadlen, args.mininsert, args.mintarglen, args.unique, args.chimeric, args.duplicates, args.noiupac, args.nbases, args.nogaps)

    def process_contigs(self, inpath, reference, outdir, coverage, flanking, maxerror, minoral, minmapqual, minreadlen, mininsert, mintarglen, unique, chimeric, duplicates, noiupac, nbases, nogaps):

        if basecaller.debug:
            pdb.set_trace()

        executor.setconfig("samtools")

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        basecaller.run_dry(SeqIO, getinput, output, makenewdir, executor, target)
        basecaller.set_debug(executor)

        # CHECK ANALYSIS SETTINGS
        caller_pars(coverage, flanking, maxerror, minoral, minmapqual, minreadlen, mininsert, mintarglen, unique, chimeric, duplicates, noiupac, nbases, nogaps)
        if maxerror < 0.002:
            maxerror = 0.002
            logging.warning("Maximum allowed allignment error was changed to 0.002. This threshold should not be lower to make the  EXONtools algorithm to work properly")

        # GET INPUT PATHS
        getinput.format(['.sam', '.bam'])
        samlist = getinput(inpath).files
        sam_files = {}
        for x in samlist:
            sam_files.update({os.path.basename(x).split("_")[0].split(".")[0]: x})

        getinput.format(['.fa', '.fasta'])
        FileList = getinput(reference).files
        libfiles = {}
        logging.info("The following libraries have matched reference and mapping files")
        for x in FileList:
            lib = os.path.basename(x).split("_")[0].split(".")[0]
            try:
                libfiles.update({lib: (sam_files[lib], x)})
                logging.info("'{0:s}' library --> {1:s} + {2:s}".format(lib,os.path.basename(sam_files[lib]),os.path.basename(x)))
            except KeyError:
                logging.error("The reference file '{0:s}' does not correspond to any BAM/SAM file".format(os.path.basename(x)))
                raise EXONtoolsError("Basecaller input error")
        del FileList
        del samlist
        del sam_files

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE TMP DIRECTORIES
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        logging.debug("IO settings: OK")

        logging.info("**************************************************************")
        logging.info("Starting the BASE CALLING procedure for each reference library")
        if basecaller.debug:
            pdb.set_trace()

        TASKS = []
        for lib in sorted(libfiles.keys(), key=natural_sort):
            liboutfile = os.path.join(output.path, lib + "_called.fa")
            TASKS.append(worker(call_bases, [lib, libfiles[lib], liboutfile, tmpdir.path, coverage, minmapqual, mintarglen, minreadlen, minoral, flanking, maxerror, mininsert, unique, chimeric, duplicates, noiupac, nbases, nogaps]))
        processes_requested = set_threads("BaseCaller", len(TASKS), basecaller.threads, userdef=False)
        pool = create_pool(processes_requested)
        results = hard_worker(run_instance, TASKS, pool)
        close_pool(pool)

        if not basecaller.keeptmp:
            tmpdir.delete()

        if not basecaller.dryrun and not results:
            logging.error("No results were produced during base calling analysis")
            raise EXONtoolsError("Base Caller Fatal Error")

        #### SAVE STATS
        if not basecaller.dryrun and basecaller.stats:
            call_stats = {}
            for result in results:
                call_stats.update(result)
            del results

            statdir = makenewdir(name="STATS", fullname="STATS")
            logging.info("Collecting stats for base calling")
            logging.info("All base calling metrics for each library will be saved to the folder:")
            logging.info(statdir.path)

            for lib in call_stats:
                insertpath = os.path.join(statdir.path,lib+".inserts")
                with open(insertpath,'w') as insertfile:
                    for insert in call_stats[lib]['inserts']:
                        insertfile.write("{0:d}\n".format(insert))

                scaffoldpath = os.path.join(statdir.path,lib+".scaffolds")
                with open(scaffoldpath,'w') as scaffoldfile:
                    for targ in call_stats[lib]['scaffolds']:
                        if call_stats[lib]['scaffolds'][targ]:
                            scaffoldfile.write("{0:s}\t{1:s}\n".format(targ,",".join(call_stats[lib]['scaffolds'][targ])))

                isoformpath = os.path.join(statdir.path,lib+".isoforms")
                with open(isoformpath,'w') as isoformfile:
                    for targ in call_stats[lib]['isoforms']:
                        if call_stats[lib]['isoforms'][targ]:
                            isoformfile.write("{0:s}\t{1:s}\n".format(targ,",".join(call_stats[lib]['isoforms'][targ])))

                chimerpath = os.path.join(statdir.path,lib+".chimeras")
                with open(chimerpath,'w') as chimerfile:
                    for targ in call_stats[lib]['chimeras']:
                        if call_stats[lib]['chimeras'][targ]:
                            chimerfile.write("{0:s}\t{1:s}\n".format(targ,",".join(call_stats[lib]['chimeras'][targ])))

                statpath = os.path.join(statdir.path, lib+"_callstats.csv")
                header = ["No", "Target","Length","BaseCov","Nreads","Npaired","Insert","GC","HTRZG","Nindels","Nbases","Scaffold","Isoform","Chimera"]
                with open(statpath, 'w') as statout:
                    csv_writer = csv.writer(statout)
                    csv_writer.writerow(header)

                    for i,targ in enumerate(sorted(call_stats[lib]["targets"], key=natural_sort)):
                        csv_writer.writerow([
                            i+1,
                            targ,
                            call_stats[lib]["targets"][targ]["length"],
                            call_stats[lib]["targets"][targ]["coverage"],
                            call_stats[lib]["targets"][targ]["nreads"],
                            call_stats[lib]["targets"][targ]["concordance"],
                            call_stats[lib]["targets"][targ]["insertsize"],
                            call_stats[lib]["targets"][targ]["gc"],
                            call_stats[lib]["targets"][targ]["heterozygosity"],
                            call_stats[lib]["targets"][targ]["indels"],
                            call_stats[lib]["targets"][targ]["nbase"],
                            call_stats[lib]["targets"][targ]["scaffold"],
                            call_stats[lib]["targets"][targ]["isoform"],
                            call_stats[lib]["targets"][targ]["chimera"]
                            ])
            logging.info("All read mapping metrics were successfully collected")


def call_bases(lib, inpath, outpath, tmppath, mincov, minmapqual, mintarglen, minreadlen, minoral, flanking, maxerror, mininsert, unique, chimeric, duplicates, noiupac, nbases, nogaps):
    """This function calls bases in the reference fasta file using the aligned reads in the SAM/BAM file,
    estimates heterozygous sites, and trims sequence ends using the coverage threshold """

    logging.info("Calling bases in '{0:s}' library".format(lib))

    libinfo = {"targets": {}, "inserts": [], "scaffolds": {}, "isoforms": {}, "chimeras": {}}

    # INITIALIZE TARGETS
    targetfile = SeqIO(inpath[1])
    targetfile.readall()
    targets = {x.name:x.seq for x in sum(targetfile.SEQS.values(),[])}
    libinfo["total_targets"] = len(targets)
    logging.info("{0:d} reference targets are provided for '{1:s}' library".format(len(targets),lib))
    del targetfile

    if inpath[0].endswith('bam'):
        logging.debug("Converting BAM to SAM file")
        sampath = os.path.join(tmppath, os.path.basename(os.path.splitext(inpath[0])[0] + ".sam"))
        samviewpars = "--threads 1 -O sam -o {0:s} {1:s}".format(sampath, inpath[0])
        run_executor(executor(
            program="SAMTOOLS",
            params=["view", samviewpars],
            conditions={"pathexists": [inpath[0]]},
            custom_arg_string=" > /dev/null 2>&1",
            quiet=True)
        )
    else:
        sampath = inpath

    # TEST PHRED
    logging.debug("Estimating phred score for '{0:s}' library".format(lib))
    if basecaller.dryrun:
        phredadjust = None
    else:
        pmode = phredmode(sampath, filetype="SAM", nchecks=100000)
        if pmode in ["L", "S"]:
            phredadjust = 33
        elif pmode in ["I", "J", "X"]:
            phredadjust = 64
        else:
            logging.error("Unable to identify PHRED score mode for further analysis")
            raise EXONtoolsError("Phred mode prediction error")
        libinfo["phred"] = phredadjust
        logging.debug("'Phred{0:d}' score is set to check the read qualities".format(phredadjust))

    # SETUP TARGET CLASS
    target.setup(
        lib = lib,                          # library name
        flank=flanking,                     # bases from read flanking regions
        mininsert=mininsert,                # minimum number of inserts to allow addition
        mincov=mincov,                      # minimum base coverage for triming
        mintarglen=mintarglen,              # minimum target length to be accepted
        phred=phredadjust,                  # type of phred score
        noiupac=noiupac,                    # use IUPAC
        nbases=nbases,                      # insert N's when low coverage
        nogaps=nogaps                       # avoid inserting N's in unmapped regions
        )

    lasttargname = ""                       # check name and presence of the previous target (switcher)
    libtarget = None                        # empty variable for future target class
    maxalign = 1
    maxmap = 1

    samfile = SeqIO(sampath, "SAM")
    for samread in samfile.read():
        #FILTER READS: mapping success / length / mapping quality / duplicates / uniqueness / chimeras / poor alignemnt
        if not samread.flag & 0x4 and samread.readlen >= minreadlen and samread.mapq >= minmapqual and ((duplicates and samread.duplicate) or not samread.duplicate) and ((unique and not (samread.secondary or samread.secs)) or not unique) and ((chimeric and samread.supplementary) or not samread.supplementary) and (samread.alqual >= -10*LN(maxerror)):

            #### UPDATE MAX VALUES FOR THRESHOLD
            if maxalign and samread.alqual > maxalign:
                maxalign = samread.alqual
            if maxmap and samread.mapq > maxmap:
                maxmap = samread.mapq

            if basecaller.stats and samread.sister and samread.flag & 0x1 and not samread.flag & 0x100 and not samread.flag & 0x800:
                try:
                    libinfo["scaffolds"][samread.target].append(samread.sister)
                except KeyError:
                    libinfo["scaffolds"][samread.target] = [samread.sister]

            if basecaller.stats and samread.secs:
                try:
                    libinfo["isoforms"][samread.target].append(samread.secs[0][0])
                except KeyError:
                    libinfo["isoforms"][samread.target] = [samread.secs[0][0]]

            if basecaller.stats and samread.suppls:
                try:
                    libinfo["chimeras"][samread.target].append(samread.suppls[0][0])
                except KeyError:
                    libinfo["chimeras"][samread.target] = [samread.suppls[0][0]]

            #### FINALIZE THE PREVIOUS TARGET AND START THE NEW TARGET INSTANCE
            if samread.target != lasttargname and lasttargname:     # runs only if target already defined
                libtarget.set_threshold(coefreq=minoral, maxqual=40, maxmap=maxmap, maxalign=maxalign)
                resultedseq = libtarget.finalize()                         # finalize the previous target

                if resultedseq:
                    with open(outpath,'a') as outfile:
                        outfile.write(">{0:s}\n{1:s}\n".format(libtarget.name,libtarget.finalseq))    # SAVE THE OUTPUT
                    if basecaller.stats:
                        libinfo["inserts"] = libinfo["inserts"] + libtarget.inserts
                        libinfo["targets"][libtarget.name]={
                            "coverage":libtarget.basecoverage,
                            "heterozygosity" : libtarget.heterozygosity,
                            "length" : len(libtarget.finalseq),
                            "nreads" : libtarget.reads,
                            "gc" : libtarget.gccont,
                            "concordance" : len(libtarget.inserts),
                            "indels" : libtarget.indels,
                            "insertsize" : libtarget.averinsert,
                            "nbase" : libtarget.nprop,
                            "scaffold" : libtarget.scaffold,
                            "isoform" : libtarget.isoform,
                            "chimera": libtarget.chimera
                        }

                #### CREATE NEW TARGET INSTANCE
                libtarget = target(
                    name=samread.target,                # target name
                    seq=targets[samread.target]         # target sequence
                )
                libtarget.add_read(samread)             # add the first read to the new target instance
                lasttargname = samread.target           # assign target name to the switcher

            ### APPEND NEW SAM READ CLASS TO THE CURRENT TARGET INSTANCE
            elif samread.target == lasttargname:
                libtarget.add_read(samread)

            #### INITIALIZE THE FIRST TARGET (THIS IS LAUNCHED JUST ONCE)
            else:
                lasttargname = samread.target           # assign target name to the switcher
                libtarget = target(
                    name=samread.target,                # target name
                    seq=targets[samread.target]         # target sequence
                )
                libtarget.add_read(samread)             # add the first read to the first target instance

    #### THIS STEP IS REQUIRED TO FINALIZE THE LAST TARGET
    try:
        if not basecaller.dryrun and not libtarget.finalized:
            libtarget.set_threshold(coefreq=minoral, maxqual=40, maxmap=maxmap, maxalign=maxalign)
            resultedseq = libtarget.finalize()                         # finalize the previous target
            if resultedseq:
                with open(outpath,'a') as outfile:
                    outfile.write(">{0:s}\n{1:s}\n".format(libtarget.name,libtarget.finalseq))    # SAVE THE OUTPUT
                if basecaller.stats:
                    libinfo["inserts"] = libinfo["inserts"] + libtarget.inserts
                    libinfo["targets"][libtarget.name]={
                        "coverage":libtarget.basecoverage,
                        "heterozygosity" : libtarget.heterozygosity,
                        "length" : len(libtarget.finalseq),
                        "nreads" : libtarget.reads,
                        "gc" : libtarget.gccont,
                        "concordance" : len(libtarget.inserts),
                        "indels" : libtarget.indels,
                        "insertsize" : libtarget.averinsert,
                        "nbase" : libtarget.nprop,
                        "scaffold" : libtarget.scaffold,
                        "isoform" : libtarget.isoform,
                        "chimera": libtarget.chimera
                    }
    except (NameError,AttributeError):
        logging.error("No mapped reads were found in sam file of '{0:s}' library".format(lib))
        raise EXONtoolsError("Empty SAM/BAM file")

    if not basecaller.dryrun and inpath[0].endswith('bam'):
        logging.debug("Deleting temporary SAM file")
        os.remove(sampath)

    if basecaller.stats:
        for targ in libinfo["scaffolds"]:
            libinfo["scaffolds"][targ] = [x for x in libinfo["scaffolds"][targ] if libinfo["scaffolds"][targ].count(x) > 5]
            libinfo["scaffolds"][targ] = sorted(list(set(libinfo["scaffolds"][targ])),key=natural_sort)
            if targ in libinfo["scaffolds"][targ]:
                libinfo["scaffolds"][targ].remove(targ)

        for targ in libinfo["isoforms"]:
            libinfo["isoforms"][targ] = [x for x in libinfo["isoforms"][targ] if libinfo["isoforms"][targ].count(x) > 10]
            libinfo["isoforms"][targ] = sorted(list(set(libinfo["isoforms"][targ])),key=natural_sort)
            if targ in libinfo["isoforms"][targ]:
                libinfo["isoforms"][targ].remove(targ)

        for targ in libinfo["chimeras"]:
            libinfo["chimeras"][targ] = [x for x in libinfo["chimeras"][targ] if libinfo["chimeras"][targ].count(x) > 20]
            libinfo["chimeras"][targ] = sorted(list(set(libinfo["chimeras"][targ])),key=natural_sort)
            if targ in libinfo["chimeras"][targ]:
                libinfo["chimeras"][targ].remove(targ)

    return {lib:libinfo}


class target():
    """This class takes one target instance and calculates its attributes based on SAM reads.
        First, each SAM read must be added with "add_read" method.
        Attributes will be calculated only after target finalization: use "finalize" method.
    """

    #### IUPAC codes that will be used for annotation of heterozygous sites
    iupac = {"AG":"R", "CT":"Y", "CG":"S", "AT":"W", "GT" :"K", "AC":"M", "CGT":"B","AGT":"D","ACT":"H","ACG":"V"}

    dryrun = False

    library = None                # library name
    flank = None                  # number of bases in flanking regions
    mininsert = None              # minimum insertion occurrences
    mincov = None                 # coverage threshold
    mintarglen = None             # minimum length of the resulted target
    phred = None                  # type of PHRED score used to calculate base qualities
    noiupac = False               # do not use IUPAC
    nbases = False                # insert N's when low coverage
    nogaps = False              # insert dashes in unmapped regions
    threshold = None              # minimum base likelihood threshold

    def __init__(self, name, seq):
        """Target attributes that will be initialized along with the target"""
        self.name = name                  # target name
        self.seq = seq                    # target original seq
        self.finalseq = None              # processed target sequence
        self.heterozygosity = None        # revealed heterozygosity after base calling
        self.nprop = None                 # N content in the resulted target seq
        self.length = len(seq)            # seq length
        self.basecoverage = None          # revealed average base coverage after base calling
        self.gccont = None                # gc content in the resulted target seq
        self.reads = 0                    # number of reads accepted for this target
        self.inserts = []                 # insert distribution for this target
        self.averinsert = 0
        self.indels = {}                  # add indel insertions tuples
        self.bases = [[] for x in range(self.length)]  # collection of called bases for each position on the target
        self.finalized = False            # checks if the target is finalized
        self.isoform = []                 # target names in valid secondary mappings
        self.scaffold = []                # scaffold mappings
        self.chimera = []                 # supplementary mappings

        """
        Base calling threshold calculation.
        We assume that each base must have sufficient sequencing quality and read mapping quality.
        Even if sequence has good quality, its mapping can be wrong.
        Default 0.2 minoral is the calibration coefficient that adjusts the sensitivity of base calling.
        Here we think that well sequenced and mapped allels will have > 20% frequencies
        """

    @classmethod
    def setup(cls, lib, flank=0, mininsert=0, mincov=0, mintarglen=0, phred=33, noiupac=False, nbases=False,nogaps=False):
        cls.library = lib
        cls.flank = flank
        cls.mininsert = mininsert
        cls.mincov = mincov
        cls.mintarglen = mintarglen
        cls.phred = phred
        cls.noiupac = noiupac
        cls.nbases = nbases
        cls.nogaps = nogaps

    @classmethod
    def set_threshold(cls, coefreq, maxqual, maxmap, maxalign):
        ncheck = 4
        allitems = [maxqual,maxmap,maxalign]
        for i,check in enumerate(allitems):
            if check <= 1:
                ncheck-=1
                allitems[i]=1
        cls.threshold = round(coefreq**ncheck * reduce2(lambda x,y: x*y,allitems),2)

    @classmethod
    def setdry(cls):
        cls.dryrun = True

    def add_read(self,samread):
        """Adds new read from SAM file"""

        if self.finalized:
            logging.error("Target '{0:s}' in '{1:s}' library is already finalized and cannot take any more reads".format(self.name, target.library))
            raise EXONtoolsError("Target read SAM error")
        elif target.dryrun:
            pass
        else:
            # correct the problem when read extends further than the reference sequence
            # by the extending reference contig with '-'s
            if samread.end > self.length:
                self.seq = self.seq + "-"*(samread.end - self.length)
                self.bases = self.bases + [[] for li in range(samread.end - self.length)]
                self.length = samread.end

            #### Update base collection with the tupple:
            #### position:(nucleotide, quality , mapping quality, alignment quality, read name)
            weighted_flanks = [j for j in range(samread.start,samread.start+target.flank)]+[j for j in range(samread.end-target.flank,samread.end)]
            for i,x,y in zip(range(samread.start,samread.end),samread.seq,samread.qual):
                if i in weighted_flanks:    # this will put greater weight on quality of marginal bases within each read
                    self.bases[i].append((x,chr(int(ord(y)/1.2)),samread.mapq//2, samread.alqual//2, samread.name))
                else:                       # all other bases are treated normally
                    self.bases[i].append((x,y,samread.mapq, samread.alqual, samread.name))

            #### COUNT READ
            self.reads +=1

            #### ADD INSERTS TO THE TARGET
            if samread.insert:
                self.inserts.append(samread.insert)

            if samread.insertions:
                for pos_insert in samread.insertions:
                    try:
                        self.indels[pos_insert].append(samread.insertions[pos_insert])
                    except KeyError:
                        self.indels[pos_insert] = [samread.insertions[pos_insert]]

            #### ADD SECONDARY MAPPINGS (SCAFFOLDS OR POSSIBLE DUPLICATES)
            if samread.sister:
                if samread.flag & 0x1 and not samread.flag & 0x4 and not samread.flag & 0x8 and not samread.flag & 0x100 and not samread.flag & 0x800:
                    self.scaffold.append(samread.sister)
                else:
                    self.isoform.append(samread.sister)

            if samread.secs:
                self.isoform.append(samread.secs[0][0])

            if samread.suppls:
                self.chimera.append(samread.suppls[0][0])

    def finalize(self):
        """This function finalizes the target instance by calculating all metrics
        Reads cannot be added after finalization
        """
        if target.dryrun:
            self.finalized = True
            return

        seq_start = 0                                    # initial start position before trimming adjusted to Python index
        seq_end = self.length - 1                        # initial end position before trimming adjusted to Python index

        #### TRIM BASES ACCORDING MINIMUM COVERAGE THRESHOLD

        # forward end trimming
        while(seq_start < self.length and len(self.bases[seq_start]) < self.mincov):
            seq_start +=1
        if seq_start == self.length:
            self.finalized = True
            return
        else:
            # reverse end trimming
            while(seq_end > seq_start and len(self.bases[seq_end]) < self.mincov):
                seq_end -=1
            if seq_end < seq_start:
                logging.error("This program bug occurred during contig trimming")
                raise EXONtoolsError("Flank region trimmimg error")

        acceptedrange = [x for x in range(seq_start,seq_end)]           # range of indexes for accepted sequence positions

        # finalize with no output if sequence length is smaller than threshold
        if len(acceptedrange) < target.mintarglen:
            self.finalized = True
            return
        else:
            #### DO THE ACTUAL BASE CALLING
            indel_count = 0
            heterozygosity = 0                                          # counter for heterozygous sites
            coverage = []                                               # list of coverage for each base
            newseq = ""                                                 # new sequence with called bases
            bases = self.bases[seq_start:seq_end]                       # take subset of base collection
            breaks = []                                                 # bases with zero coverage
            basecount = 0
            self.bases = None                                           # empty memory bucket
            for i,base in enumerate(self.seq[seq_start:seq_end]):       # do calling for each selected base position
                base_list = bases[i]                                    # create the list of selected bases
                baselen = len(base_list)                                # number of calls
                coverage.append(baselen)                                # save coverage
                basecall = ""                                           # empty variable for base call
                test = {"A":[0,0,0,0],"C":[0,0,0,0],"G":[0,0,0,0],"T":[0,0,0,0], "-":[0,0,0,0]}     # testing framework
                for x in base_list:                                     # for each selected base
                    try:
                        test[x[0]][0]+= 1                               # count test values
                        test[x[0]][1]+= ord(x[1])-target.phred          # sum test qualities for each base
                        test[x[0]][2]+= x[2]+1                          # sum read alignment qualities
                        test[x[0]][3]+= x[3]+1                          # sum read mapping
                    except KeyError:                                    # do not count gaps or Ns
                        pass
                for nuc in sorted(test.keys()):                         # calculate the likelihood for each possible nucleotide
                    stats = test[nuc]                                   # checks if the nucleotide passes threshold
                    try:
                        #### MAKE COMPARISON
                        #### proportion * weighted quality * weighted alignment quality * weighted mapping > threshold
                        if (stats[0]/baselen) * (stats[1]/baselen) * (stats[2]/baselen) * (stats[3]/baselen) > target.threshold:
                            basecall += nuc                         # add nucleotide to base call variable
                    except ZeroDivisionError:                       # skip nucleotide if one of the measures equals 0
                        pass
                    except TypeError:
                        print([(stats[0]/baselen) * (stats[1]/baselen) * (stats[2]/baselen) * (stats[3]/baselen),target.threshold])
                        raise EXONtoolsError("Threshold error")

                #### ADD INSERTIONS IF THEY ARE PRESENT
                try:
                    if len(self.indels[i+seq_start]) >= target.mininsert and len(self.indels[i+seq_start])/baselen >= 0.5:
                        indelseq_list = [indelseq[0] for indelseq in self.indels[i+seq_start]]
                        selected_indel = mostcommon(indelseq_list)
                        basecount+=len(selected_indel)
                        newseq += selected_indel
                        indel_count+=1
                except (KeyError, ZeroDivisionError):               # do nothing if there are any problems
                    pass

                #### COLLECT EMPTY COVERAGES
                basecount+=1
                if baselen ==0:
                    breaks.append(basecount)

                #### CALL AMBIGUOUS SITES
                if target.nbases and baselen < target.mincov:
                    newseq += 'N'
                elif "-" in basecall and test["-"][0]/baselen > 0.5:
                    indel_count+=1
                else:
                    basecall=basecall.replace("-","")
                    if len(basecall) == 1:
                        newseq += basecall
                    elif len(basecall) == 2 and not target.noiupac:
                        newseq += target.iupac[basecall]
                        heterozygosity +=1
                    elif len(basecall) == 2 and target.noiupac:
                        if base in basecall:
                            newseq +=base
                        else:
                            newseq +=choose(basecall)
                        heterozygosity +=1
                    elif not basecall or 2 < len(basecall) < 5:
                        newseq += 'N'
                    else:
                        logging.error("FATAL: This error should never appear. Report program bug")
                        raise EXONtoolsError("FATAL ERROR IN BASE CALLING")

            #### FILTER SHORT TARGETS
            if len(newseq.replace("-","").replace("N","")) < target.mintarglen:
                self.finalized = True
                return

            #### ESTIMATE GC CONTENT FOR THE TARGET
            self.gccont = GC(newseq, iupac=True)
            self.nprop = round(Ncheck(newseq)/len(newseq.replace("-","")),2)

            if self.inserts:
                self.averinsert = int(round(sum(self.inserts)/len(self.inserts),2))
            self.indels = indel_count
            self.heterozygosity = round(heterozygosity / len(newseq.replace("-","")),3)
            self.basecoverage = int(round(sum(coverage) / len(newseq),2))               # estimate average coverage

            #### TRIM TARGETS BY COVERAGE GAPS AND CHOOSE THE LONGEST PIECE
            if target.nogaps and breaks:
                breaks = [0]+breaks
                breaks.append(None)
                parts = [newseq[breaks[i]:breaks[i+1]] for i in range(len(breaks)-1)]
                parts = [x[:-1] if x.endswith("N") else x for x in parts]
                newseq = sorted(parts,key=len,reverse=True)[0]
                if len(newseq.replace("-","").replace("N","")) < target.mintarglen:
                    self.finalized = True
                    return

            #### CHOOSE THE MOST COMMON SCAFFOLD
            self.scaffold = [x for x in self.scaffold if x != self.name]
            if self.scaffold:
                check = mostcommon(self.scaffold)
                if self.scaffold.count(check) >= 5:
                    self.scaffold = check
                else:
                    self.scaffold = "NA"
            else:
                self.scaffold = "NA"

            #### CHOOSE THE MOST COMMON ISOFORM
            self.isoform = [x for x in self.isoform if x != self.name]
            if self.isoform:
                check = mostcommon(self.isoform)
                if self.isoform.count(check) > 10:
                    self.isoform = check
                else:
                    self.isoform = "NA"
            else:
                self.isoform = "NA"

            #### CHOOSE THE MOST COMMON CHIMERA
            self.chimera = [x for x in self.chimera if x != self.name]
            if self.chimera:
                check = mostcommon(self.chimera)
                if self.chimera.count(check) > 20:
                    self.chimera = check
                else:
                    self.chimera = "NA"
            else:
                self.chimera = "NA"

            #### REMOVE FLANKING Ns
            newseq = newseq.strip("N")

            self.finalseq = newseq                                                      # save new sequence
            self.finalized = True                                                       # turn switcher on

            return newseq                                                               # returns new sequence


def caller_pars(coverage, flanking, maxerror, minoral, minmapqual, minreadlen, mininsert, mintarglen, unique, chimeric, duplicates, noiupac, nbases, nogaps):
    """Testing mapper parameters"""

    logging.debug("Testing base caller parameters")
    if basecaller.extra:
        logging.warning("Sorry.. The 'Basecaller' program does not accept any extra arguments in '-E' option")
        raise EXONtoolsError("Argument value error")
    if coverage < 0:
        logging.error("Minumum coverage cannot be negative")
        raise EXONtoolsError("Argument value error")
    if flanking < 0:
        logging.error("Flanking region cannot be negative")
        raise EXONtoolsError("Argument value error")
    if maxerror < 0 or maxerror > 1:
        logging.error("Error rate for read mapping alignment should be set within the range between 0 and 1")
        raise EXONtoolsError("Argument value error")
    if minoral < 0 or minoral > 1:
        logging.error("Minimum read proportion for minor allels should be between 0 and 1")
        raise EXONtoolsError("Argument value error")
    if minmapqual < 0 or minmapqual > 60:
        logging.error("Minimum mapping quality should be set within the range between 0 and 60")
        raise EXONtoolsError("Argument value error")
    if minreadlen < 0:
        logging.error("Minumum read length cannot be negative")
        raise EXONtoolsError("Argument value error")
    if mininsert < 0:
        logging.error("Minumum number of insert occurences cannot be negative")
        raise EXONtoolsError("Argument value error")
    if mintarglen < 0:
        logging.error("Minumum number of insert occurences cannot be negative")
        raise EXONtoolsError("Argument value error")
    if unique:
        logging.warning("Only uniquely mapped reads will be used in base calling procedure")
    if chimeric:
        logging.warning("Chimeric read alignments will be used in base calling procedure")
    if duplicates:
        logging.warning("Read duplicates will be used in base calling procedure")
    if noiupac:
        logging.warning("IUPAC codes will NOT be used for heterozygous sites. Using the reference value instead")
    if nbases:
        logging.warning("All bases with coverage below the threshold will be changed to N's")
    if nogaps:
        logging.warning("Avoiding coverage gaps. Only the longest fragment will be produced in each contig")
