# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

description = """
                        ########################
                             CORRECT_READS
                        ########################

This command finds and fixes sequencing errors within fastq files
"""

epilog = """
Provide path to a directory or to a file with FASTQ files (-i). Optionally,
you can select single paired FASTQ files (-R1, -R2) or single unpaired files (-U).
If selecting directory, it will be parsed to identify paired files automatically.
Currently, the following file extensions are supported for cleanup analysis:
['fastq','fq','fastq.gz','fq.gz'].

All names of files with paired reads must include '_R1' for forward
and '_R2' for reverse reads. If these identifiers are not provided,
the files will be treated as unpaired.

Use '-S' or '--stats' general command option to create 'read_stats.csv' file
with the list of library names, read numbers and barcodes

'--rqc' option enables the read quality analysis using 'fastqc' and 'multiqc'

Please also cite Bayeshammer when correcting Illumina data:
'Nikolenko, S.I., Korobeynikov, A.I. and Alekseyev, M.A., 2013. BayesHammer:
Bayesian clustering for error correction in single-cell sequencing. BMC Genomics 14, S7'

Please also cite Ionhammer when Ion Torrent data:
'Ershov, V., Tarasov, A., Lapidus, A. and Korobeynikov, A., 2019. IonHammer: 
Homopolymer-space hamming clustering for IonTorrent read error correction. 
J. Comput. Biol., 26(2), 124-127.'

"""