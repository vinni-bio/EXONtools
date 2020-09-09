#! /usr/bin/env python
# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2020.
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

# PYTHON BUILT-IN MODULES
from __future__ import print_function, division
import pdb
import logging
import sys
import os
import platform
import argparse
import signal
from multiprocessing import cpu_count
from time import time
from datetime import timedelta


# EXONtools version and copyright info
VERSION = "v0.3b"
COPYRIGHT = "\nCopyright 2020, Kirill Vinnikov. All rights reserved.\nRead LISENCE.txt for BSD 3-Clause Clear License statement."


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """SETTING THE OUTPUT CLASS FOR ARGPARSE HELP MENUS"""
    pass


def main(argv):
    """Parses command line args, builds an argparse.ArgumentParser, and runs the chosen EXONtools command.
    Otherwise, prints usage."""

    options = argparse.ArgumentParser(
        prog="EXONtools",
        description="""

       ***************************************************************
       *                       EXONtools v0.3b                       *
       *            a complete pipeline for exon capture             *
       *         sequencing data analysis on non-model organisms     *
       ***************************************************************
        """,
        usage='%(prog)s.py [general options] COMMAND [command options]',
        formatter_class=CustomFormatter,
        epilog="Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues\nor send me an email to <vinni(at)hawaii.edu> with your suggestions"
    )
    options._optionals.title = "General EXONtools options"
    options.add_argument('-c', '--citations', action=EXT_citat.CITATION)
    options.add_argument('-D', '--debug', action='store_true', help="Run the debugging mode with pdb (advanced users only)", default=False, dest="debugmode")
    options.add_argument('-E', '--extra', action="store", type=str, default="", help="Auxillary command line arguments for a dependency program. See details below (advanced users only)", dest="extra", metavar="<str>")
    options.add_argument('-K', '--keeptmp', action='store_true', help="Do not delete the directory with temporary output files", default=False, dest="keeptmp")
    options.add_argument('-L', '--log', action='store_true', help="Additionally save all console output to EXONtools.log file in the current directory", default=False, dest="logmode")
    options.add_argument('-l', '--license', action=EXT_LCNS.LICENSE)
    options.add_argument('-M', '--memory', type=float, action="store", help="RAM usage limit for dependency programs (Gb)", default=4, dest="memory", metavar="<int>")
    options.add_argument('-Q', '--quiet', action='store_true', help="Hide console output except any critical errors", default=False, dest="quietmode")
    options.add_argument('-R', '--dryrun', action='store_true', help="Run the pipeline command without implementing any real actions", default=False, dest="rundry")
    options.add_argument('-S', '--stats', action='store_true', help="Write a csv file with some summary statistics for the current pipeline step", default=False, dest="stats")
    options.add_argument("-T", "--threads", type=int, action="store", help="Number of threads to run in parallel mode", default=1, dest="threads", metavar="<int>")
    options.add_argument('-v', '--version', help="Show EXONtools current version", action='version', version='%(prog)s' + VERSION + COPYRIGHT)
    options.add_argument("-W", "--warnigns", help="Show warnings only in the log output", action="store_true", default=False, dest="warnings")
    commands = options.add_subparsers(dest='action', title="Available EXONtools commands", description=EXT_help.description, prog="EXONtools.py [options]", metavar="")

    # STEPS A1 and D1: split FASTQ files with multiplexed data
    demultiplex_reads = commands.add_parser("demultiplex_reads", description=demultiplex_reads_help.description, usage="EXONtools.py [general options] demultiplex_reads [command options]", formatter_class=CustomFormatter, epilog=demultiplex_reads_help.epilog)
    demultiplex_reads._optionals.title = "'demultiplex_reads' command options"
    IOoptions = demultiplex_reads.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-R1", "--forward", action="store", help="Input path to a single FASTQ file with forward or unpaired reads", type=str, required=True, metavar="<path>", dest="forward")
    IOoptions.add_argument("-R2", "--reverse", action="store", help="Input path to a single FASTQ file with reverse reads", type=str, required=False, metavar="<path>", dest="reverse")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory", type=str, default="./DEMULTIPLEX", metavar="<path>", dest="outdir")
    Dsettings = demultiplex_reads.add_argument_group(title="Demultiplexing settings")
    Dsettings.add_argument("-b", "--barcode", action="store", required=False, help="Path to the text file with library indexes. Also launches demultiplexing using barcodes", type=str, dest="barcode", metavar="<path>")
    Dsettings.add_argument("--inseq", action="store_true", default=False, help="Enables searching barcodes within sequences", dest="inseqsearch")
    Dsettings.add_argument("--start", action="store", required=False, help="Start position of a barcode in each read sequence (with --inseq option)", type=int, dest="start", default=1, metavar="<int>")
    Dsettings.add_argument("--trim", action="store", required=False, help="Number of bases to trim additionally from each read after barcode removal (with --inseq option)", type=int, dest="trim", default=0, metavar="<int>")
    Dsettings.add_argument("-p","--pattern", action="store", required=False, help="Custom grep pattern to split reads by their name identifiers (type default to use '^@([^_:]+)[_:].*$' pattern", type=str, dest="pattern", metavar="<str>")
    Dsettings.add_argument("-t","--tolerance", help="Rename read names using library identifiers", action="store", default=0, dest="tolerance", type=int, choices=[0, 1, 2], metavar="<int>")
    miscopts = demultiplex_reads.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--rqc", help="Performs read quality check on splitted files", action="store_true", dest="rqc")
    miscopts.add_argument("--gzip", help="Compress output files", action="store_true", default=False, dest="gzoutput")
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="", metavar="<str>")
    miscopts.add_argument("--program", default="demultiplexer", action="store", choices=['demultiplexer'], help="Indicates which program to use read demultiplexing. Choices are: ['demultiplexer']", dest="program", metavar="<program name>")
    demultiplex_reads.set_defaults(command=demultiplex_reads_com.command)
    
    # STEPS A2 and D2: verify fastq format, check correspondence of paired reads, check read quality, check Illumina filters, rename reads using file name pattern, convert to FASTA format or compress files
    format_reads = commands.add_parser("format_reads", description=format_reads_help.description, usage="EXONtools.py [general options] format_reads [command options]", formatter_class=CustomFormatter, epilog=format_reads_help.epilog)
    format_reads._optionals.title = "'format_reads' command options"
    IOoptions = format_reads.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", help="Path to the directory containing read fastq files (alternative option to R1/R2/U)", type=str, required=False, metavar="<path>", dest="inpath")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory", type=str, default="./PRECLEAN", metavar="<path>", dest="outdir")
    IOoptions.add_argument("-R1", "--forward", action="store", required=False, help="Path to a single fastq file containing forward (_R1) paired reads", type=str, metavar="<path>", dest="forward")
    IOoptions.add_argument("-R2", "--reverse", action="store", required=False, help="Path to a single fastq file containing reverse (_R2) paired reads", type=str, metavar="<path>", dest="reverse")
    IOoptions.add_argument("-U", "--unpaired", action="store", required=False, help="Path to a single fastq file containing unpaired reads", type=str, metavar="<path>", dest="unpaired")
    Fsettings = format_reads.add_argument_group(title="Read format settings")
    Fsettings.add_argument("--type", action="store", help="Read name pattern. Choices are: 'Illumina', 'Torrent', EXONtools' and 'CUSTOM'", default="Illumina", choices = ["illumina","torrent", "exontools","custom"], type=str.lower, dest="pattern", metavar="<str>")
    Fsettings.add_argument("--custom", action="store", help="Custom grep pattern to capture complete read name identifier", default="", type=str, dest="customgrep", metavar="<str>")
    Fsettings.add_argument("--rename", help="Rename read names using library identifiers", action="store_true", default=False, dest="rename")
    Fsettings.add_argument("--fasta", help="Convert FASTQ to FASTA", action="store_true", default=False, dest="fq2fa")
    miscopts = format_reads.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--rqc", help="Performs read quality check on formatted files", action="store_true", dest="rqc")
    miscopts.add_argument("--gzip", help="Compress output files", action="store_true", default=False, dest="gzoutput")
    miscopts.add_argument("--skip", help="Skip FASTQ format check", action="store_true", dest="skipcheck")
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="", metavar="<str>")
    miscopts.add_argument("--program", default="readformatter", action="store", choices=['readformatter'], help="Indicates which program to use for read formatting. Choices are: ['readformatter']", dest="program", metavar="<program name>")
    format_reads.set_defaults(command=format_reads_com.command)

    # STEPS A3 and D3: find and fix sequencing errors within fastq files
    correct_reads = commands.add_parser("correct_reads", description=correct_reads_help.description, usage="EXONtools.py [general options] correct_reads [command options]", formatter_class=CustomFormatter, epilog=correct_reads_help.epilog)
    correct_reads._optionals.title = "'correct_reads' command options"
    IOoptions = correct_reads.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", help="Path to the directory containing read fastq files (alternative option to R1/R2/U)", type=str, required=False, metavar="<path>", dest="inpath")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory", type=str, default="./CORRECTED", metavar="<path>", dest="outdir")
    IOoptions.add_argument("-R1", "--forward", action="store", required=False, help="Path to a single fastq file containing forward (_R1) paired reads", type=str, metavar="<path>", dest="forward")
    IOoptions.add_argument("-R2", "--reverse", action="store", required=False, help="Path to a single fastq file containing reverse (_R2) paired reads", type=str, metavar="<path>", dest="reverse")
    IOoptions.add_argument("-U", "--unpaired", action="store", required=False, help="Path to a single fastq file containing unpaired reads", type=str, metavar="<path>", dest="unpaired")
    miscopts = correct_reads.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--torrent", help="Required for IonTorrent data", action="store_true", default=False, dest="torrent")
    miscopts.add_argument("--rqc", help="Performs read quality check on corrected files", action="store_true", dest="rqc")
    miscopts.add_argument("--gzip", help="Compress output files", action="store_true", default=False, dest="gzoutput")
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="", metavar="<str>")
    miscopts.add_argument("--program", default="hammer", action="store", choices=['hammer'], help="Indicates which program to use for read formatting. Choices are: ['readformatter']", dest="program", metavar="<program name>")
    correct_reads.set_defaults(command=correct_reads_com.command)

    # STEP A4 and D4: find and remove reads with PCR duplicates
    deduplicate_reads = commands.add_parser("deduplicate_reads", description=deduplicate_reads_help.description, usage="EXONtools.py [general options] deduplicate_reads [command options]", formatter_class=CustomFormatter, epilog=deduplicate_reads_help.epilog)
    deduplicate_reads._optionals.title = "'deduplicate_reads' command options"
    IOoptions = deduplicate_reads.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", help="Path to the directory containing read fastq files (alternative option to R1/R2/U)", type=str, required=False, metavar="<path>", dest="inpath")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory", type=str, default="./DEDUPLICATED", metavar="<path>", dest="outdir")
    IOoptions.add_argument("-R1", "--forward", action="store", required=False, help="Path to a single fastq file containing forward (_R1) paired reads", type=str, metavar="<path>", dest="forward")
    IOoptions.add_argument("-R2", "--reverse", action="store", required=False, help="Path to a single fastq file containing reverse (_R2) paired reads", type=str, metavar="<path>", dest="reverse")
    IOoptions.add_argument("-U", "--unpaired", action="store", required=False, help="Path to a single fastq file containing unpaired reads", type=str, metavar="<path>", dest="unpaired")
    miscopts = deduplicate_reads.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--mask", help="Bases to mask from start and end positions in R1 and R2 reads when searching for duplicates", action="store", type=int, nargs=4, default=[0,0,0,0], metavar="<int>", dest="skip")
    miscopts.add_argument("--phred", help="Phred score type in read quality line. Choices are: 33 and 64.", action="store", type=int, default=33, choices=[33,64], metavar="<int>", dest="phred")
    miscopts.add_argument("--rqc", help="Performs read quality check on deduplicated files", action="store_true", dest="rqc")
    miscopts.add_argument("--gzip", help="Compress output files", action="store_true", dest="gzoutput")
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="", metavar="<str>")
    miscopts.add_argument("--program", default="deduplicator", action="store", choices=['deduplicator'], help="Indicates which program to use for read formatting. Choices are: [deduplicatorNEW', 'deduplicatorOLD']", dest="program", metavar="<program name>")
    deduplicate_reads.set_defaults(command=deduplicate_reads_com.command)


    # STEPS B1 and E1: reconstruct de novo assemblies
    assemble_reads = commands.add_parser("assemble_reads", description=assemble_reads_help.description, usage="EXONtools.py [general options] assemble_reads [command options]", formatter_class=CustomFormatter, epilog=assemble_reads_help.epilog)
    assemble_reads._optionals.title = "'assemble_reads' command options"
    IOoptions = assemble_reads.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", help="Path to the directory with fastq files", type=str, required=False, dest="inpath", metavar="<path>")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory for storing resulted assemblies", type=str, default="RAW_ASSEMBLIES", dest="outdir", metavar="<path>")
    IOoptions.add_argument("-R1", "--forward", action="store", required=False, help="Path to a single fastq file containing forward (_R1) paired reads", type=str, metavar="<path>", dest="forward")
    IOoptions.add_argument("-R2", "--reverse", action="store", required=False, help="Path to a single fastq file containing reverse (_R2) paired reads", type=str, metavar="<path>", dest="reverse")
    IOoptions.add_argument("-U", "--unpaired", action="store", required=False, help="Path to a single fastq file containing unpaired reads", type=str, metavar="<path>", dest="unpaired")
    ASSoptions = assemble_reads.add_argument_group(title="Assembly settings")
    ASSoptions.add_argument("-k", "--kmers", action="store", help="List of kmers. Number of resulted assemblies will be to equal #kmers * #coverages", nargs='+', type=str, default="program default", dest="kmers", metavar="<int>")
    ASSoptions.add_argument("-c", "--coverages", action="store", help="List of required kmer coverages (default: sqrt[median] for ABySS or [1] for Trinity)", nargs='+', type=str, default='program default', dest="coverage", metavar="<int>")
    ASSoptions.add_argument("-n", "--nreads", action="store", help="Minimum number of paired reads required to reconstruct a contig", type=int, default=5, dest="npar", metavar="<int>")
    MPIoptions = assemble_reads.add_argument_group(title="Open MPI settings for ABySS and TransABySS assembly runs")
    MPIoptions.add_argument("-j", "--jobs", action="store", help="Number of cores to use per each node (mpirun) or each job (single)", default=1, type=int, dest="cores", metavar="<int>")
    MPIoptions.add_argument("-p", "--parallel", action="store", help="Total number of processes across all machines (nodes) to run with the mpi (works only together with --mpirun option)", default=0, type=int, dest="parallel", metavar="<int>")
    MPIoptions.add_argument("--mpirun", action="store", help="Path to Open-MPI mpirun. EXONtools will use a single machine if nothing is provided", required=False, type=str, dest="mpirun", metavar="<path>")
    miscopts = assemble_reads.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--suffix", action="store", help="The ending that will be automatically added to all output file names", type=str, dest="suffix", default="_rawcontigs", metavar="<str>")
    miscopts.add_argument("--program", default="abyss", action="store", choices=['abyss', 'spades', 'transabyss', 'trinity'], help="Indicates which assembler to use. Choices are: ['abyss','spades','transabyss','trinity']", dest="program", type=str.lower, metavar="<program name>")
    assemble_reads.set_defaults(command=assemble_reads_com.command)

    # STEPS B2 and E2: making the consensus assembly
    consensus_assembly = commands.add_parser("consensus_assembly", description=consensus_assembly_help.description, usage="EXONtools.py [general options] consensus_assembly [command options]", formatter_class=CustomFormatter, epilog=consensus_assembly_help.epilog)
    consensus_assembly._optionals.title = "'consensus_assembly' command options"
    IOoptions = consensus_assembly.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", help="Path to the directory with assembly files or to a single assembly file", type=str, required=True, dest="inpath", metavar="<path>")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory for storing consensus assemblies", type=str, default="CONSENSUS", dest="outdir", metavar="<path>")
    IOoptions.add_argument("--parse", action="store_true", help="Parse input files to make a separate assembly for each library", default=False, dest="parsing")
    FILTERoptions = consensus_assembly.add_argument_group(title="Filtering options")
    FILTERoptions.add_argument("-c", "--cluster", help="Contig identity level in CDHIT clusterization analysis", type=float, action="store", default=0.99, dest="cluster", metavar="<float>")
    FILTERoptions.add_argument("-l", "--length", help="Minimum contig length allowed in the output", type=int, action="store", default=200, dest="minlen", metavar="<int>")
    ASSoptions = consensus_assembly.add_argument_group(title="Assembly options")
    ASSoptions.add_argument("-a", "--overlap", help="Allowed minimum overlap proportion 'BLAT alignment/contigLength overlap' to assemble two contigs", type=float, default=0.95, action="store", dest="overlap", metavar="<float>")
    ASSoptions.add_argument("-d", "--decrement", help="Decrement in similarity (-s) per each cycle in the BLAT clustering procedure", type=float, default=0.01, action="store", dest="decrement", metavar="<float>")
    ASSoptions.add_argument("-r", "--repeat", help="Number of times to repeat decrement steps in the BLAT clustering procedure", type=int, action="store", default=2, dest="repeats", metavar="<int>")
    ASSoptions.add_argument("-s", "--similarity", help="Starting similarity threshold for accepting the results of pairwise BLAT alignment", type=float, action="store", default=0.99, dest="similarity", metavar="<float>")
    miscopts = consensus_assembly.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="_consensus", metavar="<str>")
    miscopts.add_argument("--program", default="assemblator", action="store", choices=['assemblator'], help="Indicates which program to use for the reassembling. Choices are: ['assemblator']", dest="program", metavar="<program name>")
    consensus_assembly.set_defaults(command=consensus_assembly_com.command)

    # STEPS B3 and E3: mapping reads to reference
    map_reads = commands.add_parser("map_reads", description=map_reads_help.description, usage="EXONtools.py [general options] map_reads [command options]", formatter_class=CustomFormatter, epilog=map_reads_help.epilog)
    map_reads._optionals.title = "'map_reads' command options"
    IOoptions = map_reads.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", help="Path to the directory containing read fastq files (alternative option to R1/R2/U)", type=str, required=False, metavar="<path>", dest="inpath")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory", type=str, default="MAPPEDREADS", metavar="<path>", dest="outdir")
    IOoptions.add_argument("-R1", "--forward", action="store", required=False, help="Path to a single fastq file containing forward (_R1) paired reads", type=str, metavar="<path>", dest="forward")
    IOoptions.add_argument("-R2", "--reverse", action="store", required=False, help="Path to a single fastq file containing reverse (_R2) paired reads", type=str, metavar="<path>", dest="reverse")
    IOoptions.add_argument("-U", "--unpaired", action="store", required=False, help="Path to a single fastq file containing unpaired reads", type=str, metavar="<path>", dest="unpaired")
    IOoptions.add_argument("-r", "--reference", action="store", required=True, help="Path to the directory or to a single fasta file containing reference (target) sequences", type=str, metavar="<path>", dest="reference")
    SCOptions = map_reads.add_argument_group(title="Alignment score penalties")
    SCOptions.add_argument("--mismatch", help="Mismatch penalty", type=int, action="store", default=9, dest="mismatch", metavar="<int>")
    SCOptions.add_argument("--gapopen", help="Gap open penalty", type=int, action="store", default=16, dest="gapopen", metavar="<int>")
    SCOptions.add_argument("--gapextension", help="Gap extension penalty", type=int, action="store", default=1, dest="gapextension", metavar="<int>")
    SCOptions.add_argument("--clipping", help="End clipping penalty", type=int, action="store", default=5, dest="clipping", metavar="<int>")
    SCOptions.add_argument("--discordance", help="Discordant paired read penalty", type=int, action="store", default=17, dest="discordant", metavar="<int>")
    miscopts = map_reads.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="_mapped", metavar="<str>")
    miscopts.add_argument("--program", default="bwa", action="store", choices=['bwa', 'bowtie2'], help="Indicates which program to use. Choices are: ['bwa','bowtie2']", dest="program", metavar="<program name>")
    map_reads.set_defaults(command=map_reads_com.command)

    # STEPS B4 and E4: calling bases and trimming contigs based on read coverage
    call_bases = commands.add_parser("call_bases", description=call_bases_help.description, usage="EXONtools.py [general options] call_bases [command options]", formatter_class=CustomFormatter, epilog=call_bases_help.epilog)
    call_bases._optionals.title = "'call_bases' command options"
    IOoptions = call_bases.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", required=True, help="Path to the directory with sorted BAM/SAM files or to a single sorted BAM/SAM file containing mapped reads", type=str, metavar="<path>", dest="inpath")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory", type=str, default="CALLEDBASES", metavar="<path>", dest="outdir")
    IOoptions.add_argument("-r", "--reference", action="store", required=True, help="Path to the directory or to a single fasta file containing reference (target) sequences previously used for mapping reads", type=str, metavar="<path>", dest="reference")
    BaseCall = call_bases.add_argument_group(title="Base calling settings")
    BaseCall.add_argument("-m", "--minor", help="Minimum proportion of reads required to accept the minor allele", type=float, action="store", default=0.2, dest="minor", metavar="<int>")
    BaseCall.add_argument("-c", "--coverage", help="Minimum number of reads required to retain bases in flanking regions", type=int, action="store", default=10, dest="coverage", metavar="<int>")
    BaseCall.add_argument("-f", "--flanking", help="Number of bases from each read flanking region to be judged much more strictly than all other bases", type=int, action="store", default=5, dest="flanking", metavar="<int>")
    BaseCall.add_argument("--mininsert", help="Minimum number of inserts per base required to accept it as a true insert if its occurrence exceeds 50%%", type=int, action="store", default=5, dest="mininsert", metavar="<int>")
    BaseCall.add_argument("--chimeric", action="store_true", help="Apply supplementary (chimeric) reads for base calling", default=False, dest="chimeric")
    BaseCall.add_argument("--duplicates", action="store_true", help="Allow read duplicates in base calling", default=False, dest="duplicates")
    BaseCall.add_argument("--unique", action="store_true", help="Use only uniquely mapped reads for base calling", default=False, dest="unique")
    FilterOpts = call_bases.add_argument_group(title="Read and target filtering settings")
    FilterOpts.add_argument("--minreadlen", help="Minimum read length allowed for base calling", type=int, action="store", default=50, dest="minreadlen", metavar="<int>")
    FilterOpts.add_argument("--minmapqual", help="Minimum read mapping quality allowed for base calling", type=int, action="store", default=0, dest="minmapqual", metavar="<int>")
    FilterOpts.add_argument("--maxerror", help="Maximum read alignment error allowed for base calling (NM/READLENGTH)", type=float, action="store", default=0.05, dest="maxerror", metavar="<float>")
    FilterOpts.add_argument("--mintarglen", help="Minimum sequence length in the final output", type=int, action="store", default=200, dest="mintarglen", metavar="<int>")
    SeqAnnot = call_bases.add_argument_group(title="Sequence annotation settings")
    SeqAnnot.add_argument("--nbases", action="store_true", help="Convert all bases that are below the coverage threshold (in option -c) to Ns", default=False, dest="nbases")
    SeqAnnot.add_argument("--noiupac", action="store_true", help="Do not use the IUPAC codes to designate all heterozygous sites", default=False, dest="noiupac")
    SeqAnnot.add_argument("--nogaps", action="store_true", help="Do not apply gaps instead of the bases with zero read coverage (see the details below)", default=False, dest="nogaps")
    miscopts = call_bases.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="_called", metavar="<str>")
    miscopts.add_argument("--program", default="basecaller", action="store", choices=['basecaller'], help="Indicates which program to use. Choices are: ['basecaller']", dest="program", metavar="<program name>")
    call_bases.set_defaults(command=call_bases_com.command)

    # STEPS B5 and E5: assembly annotation based on provided reference sequence annotation or database
    annotate_contigs = commands.add_parser("annotate_contigs", description=annotate_contigs_help.description, usage="EXONtools.py [general options] annotate_contigs [command options]", formatter_class=CustomFormatter, epilog=annotate_contigs_help.epilog)
    annotate_contigs._optionals.title = "'annotate_contigs' command options"
    IOoptions = annotate_contigs.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", help="Path to the directory with assembly files or to a single assembly file in fasta format", type=str, required=True, metavar="<path>", dest="inpath")
    IOoptions.add_argument("--filter", action="store", help="Optional path to the directory with reference files or to a single file in fasta format which will be used to filter out the query contigs (e.g., mtDNA genome). File extensions are used to recognize sequence types: '.pep'  - for protein data and '.fasta' for nucleotide data", type=str, metavar="<path>", default=None, dest="filtering")
    IOoptions.add_argument("--scaffolds", action="store", help="Optional path to the directory with scaffold files or to a single scaffold file which was produced by the EXONtools 'call_bases' command", type=str, metavar="<path>", default=None, dest="scaffolds")
    IOoptions.add_argument("--isoforms", action="store", help="Optional path to the directory with isoform files or to a single isoform file which was produced by the EXONtools 'call_bases' command", type=str, metavar="<path>", default=None, dest="isoforms")
    IOoptions.add_argument("--chimeras", action="store", help="Optional path to the directory with chimeras files or to a single chimeras file which was produced by the EXONtools 'call_bases' command", type=str, metavar="<path>", default=None, dest="chimeras")
    IOoptions.add_argument("-o", "--out", action="store", required=False, help="Path to the output directory for storing resulted files with annotations", type=str, default="ANNOTATION", metavar="<path>", dest="outdir")
    IOoptions.add_argument("-r", "--reference", action="store", required=True, help="Path to the file with reference annotation", type=str, metavar="<path>", dest="reference")
    DBoptions = annotate_contigs.add_argument_group(title="Sequence format settings")
    DBoptions.add_argument("-q", "--qtype", action="store", help="Type of the query file. Choices are: ['nucl', 'prot']", choices=["nucl", "prot"], required=True, type=str, dest="query", metavar="<str>")
    DBoptions.add_argument("-t", "--ttype", action="store", help="Type of the target reference file. Choices are: ['nucl', 'prot']", choices=["nucl", "prot"], required=True, type=str, dest="target", metavar="<str>")
    DBoptions.add_argument("-g", "--gencode", action="store", help="NCBI genetic code ID", choices=[1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31], default=1, type=int, metavar="<int>", dest="gencode")
    DBoptions.add_argument("-d", "--database", action="store", choices=["swissprot", "ensembl", "exontools", "custom"], help="Type of reference database. Choices are: ['swissprot','ensembl','exontools','custom']", required=True, dest="database", metavar='<str>')
    DBoptions.add_argument("--grep", action="store", help="Custom grep pattern to parse reference IDs. Sequence id (1st parentheses) and gene ID (2nd parentheses). Example: 'DBID_(seqID)_(geneID)_ANYCHARS'", required=False, dest="custom", metavar='<grep>')
    PROCESSoptions = annotate_contigs.add_argument_group(title="Query contig preprocessing settings")
    PROCESSoptions.add_argument("-l", "--length", action="store", help="Minimum length for contigs to be accepted for annotation analysis", default=150, type=int, metavar="<int>", dest="minlen")
    PROCESSoptions.add_argument("-c", "--cluster", action="store", help="Sequence identity threshold for initial clusterization. Default 1.0 value will cluster only identical sequences", default=1.00, type=float, metavar="<float>", dest="cluster")
    PROCESSoptions.add_argument("--orf", action="store_true", help="Each query contig will be searched for long ORFs", default=False, dest="orf")
    PROCESSoptions.add_argument("--translate", action="store_true", help="Each query contig will be searched for long ORFs and translated to AA sequence", default=False, dest="translating")
    ALIGNoptions = annotate_contigs.add_argument_group(title="Alignment filtering settings")
    ALIGNoptions.add_argument("-a", "--alignment", action="store", help="Minimum alignment length required for contig annotation", default=20, type=int, metavar="<int>", dest="overlap")
    ALIGNoptions.add_argument("-e", "--evalue", action="store", help="Minimum e-value for filtering blast matches", default=1e-10, type=float, metavar="<float>", dest="evalue")
    ALIGNoptions.add_argument("-s", "--similarity", action="store", help="Minimum similarity between query and target sequences (proportion of match/mismatched bases in the blast alignment)", default=0.80, type=float, metavar="<float>", dest="similarity")
    miscopts = annotate_contigs.add_argument_group(title="Miscellaneous command options")
    # miscopts.add_argument("--sort", action="store", help="Select field to sort annotation gff file", choices=["seqid","score", "gene"], default=1, type=int, metavar="<int>", dest="sortout")
    miscopts.add_argument("--norename", action="store_true", help="Keep the original names of annotated contigs", default=False, dest="norename")
    miscopts.add_argument("--oneway", action="store_true", help="Convert all negative strand contigs to positive strand", default=False, dest="oneway")
    miscopts.add_argument("--organism", action="store", help="Species name to be provided in the annotation (organism field)",required = False, type=str, metavar="<str>", dest="organism")
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="_annotated", metavar="<str>")
    miscopts.add_argument("--program", required=False, default="annotator", action="store", choices=['annotator'], help="Indicates which program to use. Choices are: ['annotator']", dest="program", metavar="<program name>")
    annotate_contigs.set_defaults(command=annotate_contigs_com.command)

    # STEP B6a: assembly evaluation (contiguity)
    evaluate_assembly = commands.add_parser("evaluate_assembly", description=evaluate_assembly_help.description, usage="EXONtools.py [general options] evaluate_assembly [command options]", formatter_class=CustomFormatter, epilog=evaluate_assembly_help.epilog)
    evaluate_assembly._optionals.title = "'evaluate_assembly' command options"
    evaluate_assembly.add_argument("-i", "--in", action="store", help="Path to the directory with assembly files or to a single assembly file", type=str, required=True, dest="inpath", metavar="<path>")
    evaluate_assembly.add_argument("-o", "--out", action="store", help="Path to the output directory for storing stat file.", type=str, default=os.path.realpath(os.path.curdir), dest="outdir", metavar="<path>")
    evaluate_assembly.add_argument("-m", "--minlen", action="store", help="Minimum contig length for evaluation", type=int, default=500, dest="minlen", metavar="<int>")
    evaluate_assembly.add_argument("--program", default="assass", action="store", choices=['assass'], help="EXONtools program for assembly quality assessment", dest="program", type=str.lower, metavar="<program name>")
    evaluate_assembly.set_defaults(command=evaluate_assembly_com.command)

    # STEP B6b: assembly evaluation (read mapping)
    evaluate_mapping = commands.add_parser("evaluate_mapping", description=evaluate_mapping_help.description, usage="EXONtools.py [general options] evaluate_mapping [command options]", formatter_class=CustomFormatter, epilog=evaluate_mapping_help.epilog)
    evaluate_mapping._optionals.title = "'evaluate_mapping' command options"
    evaluate_mapping.add_argument("-i", "--in", action="store", help="Path to the directory with BAM/SAM files or to a single BAM/SAM file", type=str, required=True, dest="inpath", metavar="<path>")
    evaluate_mapping.add_argument("-o", "--out", action="store", help="Path to the output directory for storing stat file", type=str, default=os.path.realpath(os.path.curdir), dest="outdir", metavar="<path>")
    evaluate_mapping.add_argument("-l", "--minlen", action="store", help="Minimum read length (CIGAR) to include SAM record", type=int, default=0, dest="minlen", metavar="<int>")
    evaluate_mapping.add_argument("-q", "--minqual", action="store", help="Minimum read average base quality to include SAM record", type=int, default=0, dest="minqual", metavar="<int>")
    evaluate_mapping.add_argument("--program", default="asssam", action="store", choices=['asssam'], help="EXONtools program for mapping quality assessment", dest="program", type=str.lower, metavar="<program name>")
    evaluate_mapping.set_defaults(command=evaluate_mapping_com.command)

    # STEP C1: search exon boundaries
    search_exons = commands.add_parser("search_exons", description=search_exons_help.description, usage="EXONtools.py [general options] search_exons [command options]", formatter_class=CustomFormatter, epilog=search_exons_help.epilog)
    search_exons._optionals.title = "'search_exons' command options"
    IOoptions = search_exons.add_argument_group(title="Input/Output settings")
    IOoptions.add_argument("-i", "--in", action="store", help="Path to the directory or to a single file with annotated assembly", type=str, required=True, dest="inpath", metavar="<path>")
    IOoptions.add_argument("-g", "--gff", action="store", required=True, help="Path to the directory or to a single gff file with assembly annotation", type=str, metavar="<path>", dest="gff")
    IOoptions.add_argument("-o", "--out", action="store", help="Path to the output directory for storing results", type=str, default="EXONS", dest="outdir", metavar="<path>")
    IOoptions.add_argument("-r", "--reference", action="store", required=True, help="Path to the directory or to a single fasta file with reference genomes", type=str, metavar="<path>", dest="reference")
    ALIGNoptions = search_exons.add_argument_group(title="Filtering settings")
    ALIGNoptions.add_argument("-a", "--alignment", action="store", help="Minimum blast alignment length required for exon prediction", default=20, type=int, metavar="<int>", dest="overlap")
    ALIGNoptions.add_argument("-e", "--evalue", action="store", help="Minimum e-value for filtering blast matches", default=1e-6, type=float, metavar="<float>", dest="evalue")
    ALIGNoptions.add_argument("-l", "--lag", action="store", help="Minimum distance between start positions of two exons (e.g., minumum exon length in predictions)", default=60, type=int, metavar="<int>", dest="lag")
    ALIGNoptions.add_argument("-s", "--similarity", action="store", help="Minimum alignment similarity between query contig and reference region", default=0.70, type=float, metavar="<float>", dest="similarity")
    miscopts = search_exons.add_argument_group(title="Miscellaneous command options")
    miscopts.add_argument("--chimeric", action="store_true", help="Parse chimeric contigs", dest="chimeric", default=False)
    miscopts.add_argument("--isomeric", action="store_true", help="Filter contig isomers leaving only the longest sequence", dest="isomeric", default=False)
    miscopts.add_argument("--nagenes", action="store_true", help="Remove annotations with the unknown gene names (e.g., gene=NA)", dest="naexclude", default=False)
    miscopts.add_argument("--orf", action="store_true", help="Include only coding exon sequences in the output fasta file (reference gff file must have CDS annotations)", dest="orfcheck", default=False)
    miscopts.add_argument("--unique", action="store_true", help="Choose the longest contig from each annotation for exon prediction", dest="unique", default=False)
    miscopts.add_argument("--suffix", action="store", help="Ending that will be automatically added to all output file names", type=str, dest="suffix", default="_exons", metavar="<str>")
    miscopts.add_argument("--program", required=False, default="exosearcher", action="store", choices=['exosearcher'], help="Indicates which program to use. Choices are: ['exosearcher']", dest="program", metavar="<program name>")
    search_exons.set_defaults(command=search_exons_com.command)

###############################################################################################

    # PARSE COMMAND ARGUMENTS
    args, unknown = options.parse_known_args()

    # WELCOMING SCREEN
    if args.quietmode and not args.debugmode:
        pass
    else:
        print("\n-----------------------------------\nWelcome to the EXONtools pipeline!\n-----------------------------------\n")

        if args.rundry:
            print("IMPORTANT: You selected the DRY RUN mode for the current run of the EXONtools. No actual analyses will be performed\n")

        if args.debugmode:
            print("IMPORTANT: RUNNING THE DEBUGGING MODE\n")
        print('The EXONtools command:')
        print(os.path.realpath(sys.argv[0]), " ".join(sys.argv[1:]))

        print("\nYOUR SYSTEM INFORMATION")
        print("EXONtools version: {}".format(VERSION[1:]))
        print("Python version: {}".format(platform.python_version()))
        print("OS version: {}".format(platform.platform()))
        print("Number of cores: {}".format(cpu_count()))
        print("Memory: {0:.1f} GB\n".format(float(psutil.virtual_memory()[0]) / 1000000000))

    # SETTING UP THE EXONTOOLS LOGGER
    if args.warnings and not args.debugmode:
        logger = setup_logger(level="WARNING", logmode=args.logmode, quietmode=args.quietmode)
        logger.warning("The output verbosity is set to show warnings only")
        loglevel = "WARNING"
    elif args.debugmode:
        logger = setup_logger(level="DEBUG", logmode=args.logmode)
        logger.warning("Launching the EXONtools debugging mode")
        loglevel = "DEBUG"
        pdb.set_trace()
    else:
        logger = setup_logger(logmode=args.logmode, quietmode=args.quietmode)
        loglevel = "INFO"

    # START PIPELINE
    logger.info("Starting the EXONtools pipeline")
    start_time = time()
    logger.debug("Setting up the pipeline start time: OK")

    # WRITE JSON FILES
    log2json(level=loglevel, logmode=args.logmode, quietmode=args.quietmode)
    logger.debug("Writing log config json files: OK")

    args2json(args)
    logger.debug("Writing json file with command parameters: OK")

    new_logger_settings()

    # Check processors
    if args.threads > cpu_count() or args.threads <= 0:
        logging.error("Please choose the correct number of threads. Maximum allowed = {max_thr:d}\n".format(max_thr=cpu_count()))
        raise EXONtoolsError("Please correct the number of threads (-t option). See the last log message")
    logging.debug("Check that the number of threads is real: OK")
    # Check OS
    if sys.platform.startswith("win"):
        raise EXONtoolsError("Sorry, but the EXONtools pipeline currently doesn't work on Windows systems. You may try to use Linux instead.")
    # Check Python
    if not testpy():
        raise EXONtoolsError("Sorry, but the EXONtools pipeline currently doesn't work with any python version lower than 2.7 or 3.5")

    logging.debug("All system compatibility checks: OK")

    if args.logmode:
        logpath = os.path.realpath(os.path.join(os.path.curdir, "EXONtools.log"))
        logging.warning("The logging file mode is turned on")
        logging.warning("All logs will be written to {}".format(logpath))
    if args.quietmode:
        logging.debug("The quiet mode is automatically turned off during debugging")
    if args.keeptmp:
        logging.warning("The temporary directory will not be deleted after completion")
    if unknown:
        logging.warning("Ignoring the following unknown args: {un_arg:s}".format(un_arg=', '.join(unknown)))

    logging.debug("The initial pipeline commands are:\n {}".format(args.__dict__))
    signal.signal(signal.SIGTSTP, signal.SIG_IGN)
    logging.debug("Setting up signals: OK")

    # GET the EXONtools command object and RUN it
    try:
        cmd = args.command
    except AttributeError:
        logging.error("The EXONtools command was not provided in the arguments. This error should never occur.")
        raise EXONtoolsError
    logging.debug("Setting up the EXONtools command: OK")
    logging.debug("Ready to execute the EXONtools command '{}'".format(args.action))
    if loglevel == "DEBUG":
        pdb.set_trace()
    cmd(args).execute_command()

    # FINAL WORDS
    logging.debug("EXONtools command '{command:s}' has succesfully finished the analysis".format(command=args.action))
    if not args.quietmode:
        print("\nThe '{prog:s}' program has finished analysis in {time:s}".format(prog=args.program,time=str(timedelta(seconds=(time() - start_time)))))
    logging.debug("Reached the end of the EXONtools pipeline")
    if not args.quietmode:
        print("\nThank you for using the EXONtools pipeline!\nPlease cite us in your publications.")
###############################################################################################


if __name__ == "__main__":
    # Import/Install Extra Python MODULES
    from mains.EXT_validator import testpy, import_module_test

    if len(sys.argv) == 1:
        sys.exit("Welcome to the EXONtools pipeline!\nPlease use -h/--help option to read the help menu.")
    elif sys.argv[1].lower() == "force_install":
        import_module_test("psutil", force_install=True)
        import_module_test("python-igraph", force_install=True)
        import_module_test("configparser", force_install=True)
        import_module_test("importlib", force_install=True)
        import_module_test("numpy", force_install=True)
        import_module_test("decimal", force_install=True)
        sys.exit()
    elif sys.argv[1].lower() == "test":
        CMD = "find " + os.path.join(os.path.dirname(os.path.abspath(__file__)), "utest") + " -name 'test*' -type f -print -exec echo \; -exec python '{}' -v \; -exec echo \;"
        os.system(CMD)
        sys.exit()
    elif sys.argv[1].lower() == "setdefaultconfig":
        import_module_test("configparser")
        from mains.EXT_executor import DefaultConfig
        DefaultConfig()
        print("\nThe default configuration file is written to ../EXONtools/src/dependencies.ini\n")
        sys.exit()
    else:
        import_module_test("psutil")
        import_module_test("configparser")
        import_module_test("importlib")
        import_module_test("numpy")
        import_module_test("decimal")
        import psutil

        # EXONtools MODULES
        from mains.EXT_logger import setup_logger, log2json, args2json, new_logger_settings
        from mains.EXT_errors import EXONtoolsError
        from chelp import *
        from coms import *

    # Start EXONtools pipeline
    main(sys.argv)
