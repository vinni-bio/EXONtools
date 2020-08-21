# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """

                    #########################
                            MAP_READS
                    #########################

This command maps reads to a reference genome or to an assembly"""

epilog = """
Provide path to the folder (-i) or to single paired and/or unpaired
FASTQ files (-R1, -R2, -U) containing sequencing reads
that have to be mapped on the reference contigs (-r). Currently,
the following file extensions are supported for FASTQ files
in the mapping procedure: ['.fastq','.fq', '.fastq.gz', '.fq.gz'].

Provide path to the folder with FASTA files containing reference (target)
sequences (-r) or path that leads directly to a single FASTA file (-r).
Currently, the following extensions are supported for FASTA files with
references: ['.fasta','.fa'].

Please verify that all file names (both FASTQ and FASTA files) include
library names before the first underscore (e.g., libraryname_paired_R1.fq),
which will be used for automatic identification of which file belongs
to which library. Paired fastq files must also include '_R1' and '_R2'
identifiers in their names, otherwise they will be treated
as unpaired files.

ALIGNMENT SCORE PENALTIES:
--mismatch (default: 9)
--gapopen (default: 16)
--gapextension (default: 1)
--clipping (default: 5)  BWA only
--discordance (default: 17) BWA only

Use '-S/--stats' general option to create 'mapping_stats.csv' file containing all
information about depth of coverage and many other useful mapping stats.
The same output can be retrieved later by using the EXONtools command 'evaluate_mapping'.

Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues
or shoot me an email to <vinni(at)hawaii.edu> with your suggestions
"""
