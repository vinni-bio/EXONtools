# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

description = """
                        ###########################
                             DEMULTIPLEX_READS
                        ###########################

This command splits a single FASTQ file with multiplexed sequencing data. The file
will be subdivided into separate child files. Each child file will include only
those sequences that correspond to a single library identified by their indexes/barcodes
or by a user-specified read name pattern. Unrecognized reads will be
saved separately to 'Undefined_reads.fq' file.
"""

epilog = """
Provide path to a single input file (-R1) with FASTQ reads. If paired sequences are availalble
in the same order, please add in another input file (-R2). Currently, the following 
file extensions are supported for demultiplexing: ['fastq','fq','fastq.gz','fq.gz'].

The program can apply a user-provided regular expression to split reads by
their name identifiers. Within your regular expression pattern, please indicate which
part should be used to compare reads by puting it within parentheses.
For example, in the following grep pattern '^@([^_:]+)[_:].*$':
the parentheses will capture a library name in read sequences (default).
Call this option by using '-p' or '--pattern' flags.

Alternatively, the program can demultiplex libraries based on their indexes
provided in a separate text file. This will only work if the barcode
option (-b) is activated. A text file with library indexes must include
the first column with library names and the second column with corresponging
barcode sequence. If dual indexes were used in your analysis, please add them 
a third column. All columns must be space or tab separated. 
Use '--inseq' option to search provided barcodes within your read sequences.

You can also allow barcodes to have one or two random changes. Use '-t'/'--tolerance'
option to set the number of possible changes. By default this is 0,
assuming the exact match.

Use '--gzip' flag if you want to compress your output files.
IMPORTANT: File compression will greatly increase the overall running time.
But as a trade-off, it will save your disk space
"""
