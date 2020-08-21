# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """
                        ##########################
                               FORMAT_READS
                        ##########################

This command is used to rename all sequencing reads according
to their file names, to convert FASTQ to FASTA format,
to remove all reads that didn't pass Illumina filters
and provides some other useful tools.

"""

epilog="""
Provide path to the folder containing FASTQ files (-i) or to single
paired and/or unpaired FASTQ files (-R1, -R2, -U).
Currently, the following file extensions are supported for formatting:
['fastq','fq','fastq.gz','fq.gz'].

All names of files with paired reads must include '_R1' for forward
and '_R2' for reverse reads. If these identifiers are not provided,
the files will be treated as unpaired.

When optional '--rename' flag is added to the command, all reads will be
renamed according to their library identifiers provided in their file names
before the first underscore ('_'). All previous read names will be saved
to file 'LIBNAME_names.dat'

Use '--gzip' flag if you want to compress your output files.
IMPORTANT: File compression will greatly increase
the overall time for running the command.
But as a trade-off, it will significantly save your disk space

Use '-S' or '--stats' general command option to create a'preclean_stats.csv' file
with the list of library names, read numbers and barcodes

'--fastqc' flag runs the FastQC analysis on precleaned fastq files
"""
