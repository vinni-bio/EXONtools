# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """
                        ###########################
                             DEMULTIPLEX_READS
                        ###########################

This command splits FASTQ files with multiplexed sequencing data. Each FASTQ file
with multiplexed reads will be subdivided into separate child files. Each child file
will include only those sequences that correspond to a single library identified by
their indexes or by a user-specified read name pattern. Unrecognized reads will be
saved separately.
"""

epilog = """
Provide path to the directory containing FASTQ files (-i), to a single file (-i)
or to files with paired and unpaired reads (-R1, -R2, -U). Currently, the following
file extensions are supported for formatting: ['fastq','fq','fastq.gz','fq.gz'].

The program can apply user-provided 'regular expression' (grep) pattern to split reads by
their name identifier. Within your regular expression pattern, please indicate which
part should be used to compare reads by puting it within parentheses.
Provide a second pair of parentheses if you also would like to split FASTQ files
by read direction. No more than two pairs of parentheses are allowed.
For example, in the following grep pattern '^@([^_:]+)[_:].*$' (default):
the parentheses will capture a library name in read sequences.
Call this option by using '-p'/'--pattern' flag.

Alternatively, the program can map all libraries based on their indexes
provided in a separate text file. This will only work if the barcode
option (-b) is activated. A text file with library indexes must include
the first column with library names and the second column with corresponging
barcode sequence. If dual indexes were used in your analysis, please add a third column.
All columns must be space or tab separated. Use '--inseq' flag to search barcodes
within read sequences.

You can also allow barcodes to have up to two random changes. Use '-t'/'--tolerance'
flag to set the number of possible changes. By default this is zero,
assuming the exact matching.
"""
