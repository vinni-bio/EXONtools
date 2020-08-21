# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """

                    ###########################
                        CONSENSUS_ASSEMBLY
                    ###########################

This command takes multiple raw assemblies of each library and reconstructs
their consensus assemblies by clustering and reassembling overlapping contigs
"""

epilog = """
Provide path to the folder with assemblies in fasta files
or to a single assembly fasta file (-i). Currently, the following
file extensions are supported: ['.fasta','.fa'].

When using '--parse' option, please verify that all file names include
the corresponding library identifiers before the first underscore
(e.g., libraryname_assembly_kmer32.fasta). They will be used for
automatic concatenation of asseblies related to the same library.

Consensus_assembly command allows to adjust two different criteria
for contig reassembling:
1) similarity (i.e., proportion of non-matched bases to the length of alignment);
2) overlap threshold (i.e., proportion of aligned bases to the length of overlap).

You can set the number of BLAT clustering cycles and the level of similarity decrement
with each cycle. Please note that the program will always add one additional BLAT step
at the end of the analysis that will use the exactly same similarity value
as was used on the previous step. For example, if you use repeat = 2
with starting similarity = 0.99 and decrement = 0.01, then you will end up with one
BLAT cycle having similarity = 0.99 and two BLAT cycles with similarity = 0.98.

Use '-S/--stats' general option to create 'consensus_assembly_stats.csv' file containing
the information about number of contigs resulted on each step of the analysis.

Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues
or send me an email to <vinni(at)hawaii.edu> with your suggestions
"""
