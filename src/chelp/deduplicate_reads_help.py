# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

description = """
                        ###########################
                             DEDUPLICATE_READS
                        ###########################

This command finds and removes reads with PCR duplicates in fastq files leaving
just one read with the highest sequencing quality or largest size.
PCR duplicates are recognized as identical sequences
"""

epilog = """
This command only removes read duplicates within each single file
or within paired files but not across all files in the input directory.
As a result each final read sequence within a file will be unique only within
that particular file/paired file, but can be also present in other files.
To ensure that all sequences are uniqe across the entire dataset,
you should merge all files into a single file using 'cat' command.

Provide path to a directory or to a file with FASTQ files (-i). Optionally,
you can select single paired FASTQ files (-R1, -R2) or a single unpaired files (-U).
If selecting directory, it will be parsed to identify paired files automatically.
Currently, the following file extensions are supported for cleanup analysis:
['fastq','fq','fastq.gz','fq.gz'].

All names of files with paired reads must include '_R1' for forward
and '_R2' for reverse reads. If these identifiers are not provided,
the files will be treated as unpaired.

You can mask/skip some bases at the beginning and/or at the end of forward (R1)
and/or reverse (R2) reads while searching duplicates. Please use '--mask' option,
providing the number of bases for masking in the following way:
'--mask Nbases_in_R1_start Nbases_in_R1_end Nbases_in_R2_start Nbases_in_R2_start'
Example: '--mask 10 10 10 10'

If you use --mask option, provide all four values even if your input does not include
any files with reverse reads (R2) (use 0s for non-existing files).
Example for a single-end fastq file: '--mask 13 20 0 0'
This pattern will skip the first 13 and the last 20 bases in each read,
and only remaining region will be used for searching duplicates.

Use '-S' or '--stats' general command option to create 'read_stats.csv' file
with the list of library names, read numbers and barcodes

'--rqc' option enables the read quality analysis using 'fastqc' and 'multiqc'


"""