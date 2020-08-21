# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """
                        #########################
                               SEARCH_EXONS
                        #########################

This command finds exon boundaries within each annotated contig based on
provided reference mappings

"""

epilog = """
Provide path to the FOLDER containing fasta files OR to a single fasta FILE
with annotated contigs produced on the Step B5 (-i). Currently, the following
file extensions are supported for query and reference files: ['.fasta','.fa', '.fna'].

Provide path to the FOLDER containing gff3 files OR to a single gff3 FILE
with annotation information produced on the Step B5 (-g). Currently, the following
file extensions are supported for query and reference files: ['.gff','.gff3', '.gtf'].

Provide path to the FOLDER containing fasta files OR to a single fasta FILE
with reference genomes that will be used to map contigs (-r). Genomes are
not required to represent closely related species to your organism but should be
at least within the same class or order. Currently, the following
file extensions are supported for query and reference files: ['.fasta','.fa', '.fna'].

Please verify that all assembly files with query contigs have names that include
library names before the first underscore (e.g., LIBRARYNAME_assembly.fasta),
which will be used as assembly names. IMPORTANT: these library identifiers
will be used to associate annotation files with corresponding contig libraries.

You can adjust the following options to improve the quality of exon prediction:
1) minimum exon length (-l)
2) e-value for filtering blast results (-e)
3) minimum blast alignment length (-a)
4) blast alignment similarity threshold (-s)

Use '-S/--stats' general option to create 'search_exons_stats.csv' file
containing annotation stats for each library.

Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues
or shoot me an email to <vinni(at)hawaii.edu> with your suggestions

"""
