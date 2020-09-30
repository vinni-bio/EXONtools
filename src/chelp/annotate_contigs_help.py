# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """
                        #########################
                            ANNOTATE_CONTIGS
                        #########################

This command takes any assembly (or multiple assemblies) and annotates selected contigs
based on provided reference annotation

"""

epilog = """
Provide path to the FOLDER containing fasta files OR to a single fasta FILE
with query contigs that you wish to annotate (-i). Provide path to the reference
fasta FILE (-r). Currently, the following file extensions are supported
for query and reference files: ['.fasta','.fa','.pep','.prot'].

Please verify that all assembly files with query contigs have names that include
library names before the first underscore (e.g., LIBRARYNAME_assembly.fasta),
which will be used as assembly names. IMPORTANT: these library identifiers
will be used to associate optional files (i.e., scaffolds, isoforms, chimeras)
with corresponding contig libraries.

Define the type of query (-q) and target (-t). You may also choose to search
for the best ORFs in each query contigs and then use them for annotation:
'--orf' will trim your input assemblies.
If your target type is 'prot', than all ORFs will be translated to peptide sequences.

The following alignemnt algorithms will be automatically selected for annotation:
blastn:  both query and target are 'nucl' sequences (+ORF)
blastp:  both query and targer are 'prot' sequences (query ORFs are translated)
blastx:  if query has 'nucl' type sequences and target has 'prot' type sequences
tblastn: if query has 'prot' type sequences and target has 'nucl' type sequences

You can change the following options to improve the quality of annotation:
1) minimum length of accepted contigs for annotation (-l)
2) identity threshold for clustering input assemblies and targets (-c)
3) e-value for filtering blast results (-e)
4) minimum alignment length for contig annotaiton (-m)
5) blast alignment similarity threshold (-s)

Don't use any other parentheses in your custom grep pattern except those that
mean to capture sequence and gene identifiers.

Use '--norename' option to avoid renaming annotated sequences in fasta file

Use '-S/--stats' general option to create 'annotation_stats.csv' file
containing annotation stats for each library.

Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues
or shoot me an email to <vinni(at)hawaii.edu> with your suggestions

"""
