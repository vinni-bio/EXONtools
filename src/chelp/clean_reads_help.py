# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.


description = """

                    #########################
                            CLEAN_READS
                    #########################

This command deletes reads with low average quality and reads, which size is
lower than the required minimum. It also trims adapter residuals
and bases with low quality"""

epilog = """
Provide path to a directory or to a file with FASTQ files (-i). Optionally,
you can select single paired FASTQ files (-R1, -R2). If selecting directory,
it will be parsed to identify paired files automatically.
Currently, the following file extensions are supported for cleanup analysis:
['fastq','fq','fastq.gz','fq.gz'].

All names of files with paired reads must include '_R1' for forward
and '_R2' for reverse reads. If these identifiers are not provided,
all files will be treated as unpaired.

To set an optional adapter removal, use '--adapters' flag and provide
a path to the file containing adapter sequences in FASTA format.

Read trimmimg will be performed before size selection and quality filtering.

If chosen, adapter removal will be performed before all other cleaning procedures.

Use '-S' or '--stats' as general EXONtools option to create a 'cleanup_stats.csv'
file with a list of library names and their quality stats.

'--rqc' flag starts read quality check analysis on cleaned FASTQ files.

'--program' flag allows you to choose among 'trimmomatic' (default), 'bbduk'
and 'cutadapt' dependencies.

"""
