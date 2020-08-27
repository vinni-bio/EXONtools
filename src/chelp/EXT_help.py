# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.


description = """
    align_stacks            - align annotated contigs within each target cluster (Step E7)
    annotate_contigs        - annotate contigs based on a reference database or annotated genome (Steps B5, E5)
    assemble_reads          - reconstruct de novo assemblies (Steps B1, E1)
    call_bases              - trim contigs and call bases (Steps B4, E4)
    call_snps               - call and filter SNPs from all alignments (Step E10)
    clean_reads             - filter and trim low quality reads, remove adapter residuals (Steps A5, D5)
    clean_stacks            - filter and trim raw aligments (Step E8)
    consensus_assembly      - merge assemblies and reduce redundancy of the resulted consensus assembly (Steps B2, E2)
    correct_reads           - correct sequencing errors in reads (Steps 3A, 3D)
    decontaminate_reads     - filter contaminated reads based on provided contaminant genomes (Steps A8, D8)
    deduplicate_reads       - remove read duplicates resulted from PCR amplification, i.e. PCR bias (Steps A4, D4)
    demultiplex_reads       - group reads based on barcode identifiers or read name patterns (Steps A1, D1)
    design_baits            - desig oligonucleotide sequences for exon capture baits (Step C3)
    evaluate_assembly       - assess contingency quality of assemblies (Step B6a)
    evaluate_mapping        - assess read mapping quality of assemblies (Step B6b)
    evaluate_reads          - assess sequencing read quality (Stages A & D)
    filter_reads            - filter reads that do not pass user-defined threshold for allowed fraction of
                              repeated bases or 'N' bases, i.e. low-complexity reads (Steps A6, D6)
    format_reads            - verify fastq format, check correspondence of paired reads,
                              check read quality, check Illumina filters, rename reads using
                              file name pattern, convert to FASTA format or compress files (Steps A2, D2).
    map_reads               - map reads to assembly (Steps B3, E3)
    map_exons               - check presence/occurences of exons in the reference genome/pseudoreference (Step C2)
    merge_reads             - merge paired reads if they have an overlap (Steps A7, D7)
    search_exons            - predict exon boundaries within contigs (Step C1)
    split_stacks            - split stacks by exon/intron boundaries(Step E9)
    stack_exons             - cluster exons using their annotations (Step E6)

For getting the detailed information about each command, provide the help menu option after command name:
./EXONtools.py COMMAND --help

Please check that the names of your sequence files follow the EXONtools rules of library name format:
unpaired files: LIBNAME_anyletters.fq or LIBNAME.fasta (e.g., 'ST-1_unpaired.fq')
paired files: LIBNAME_R1.fq or LIBNAME_R2.fq (e.g., 'ST-1_R1_paired.fq', 'ST-1_R2_paired.fq')
'R1' and 'R2' annotations are mandatory requirments for the files with paired reads.
Otherwise they will be considered as unpaired. A single library can only have one file
with unpaired reads and two files with paired reads (R1 & R2).

Step IDs correspond to the EXONtools pipeline (check the tutorial shown on the EXONtools github page).

You can use only one command from the list shown above for each run of the EXONtools pipeline.
To run several commands sequentially, create a bash script file as shown in 'ToyData' examples.

For some commands, you may wish to select a specific dependency program from the list
by using the '--program program_name' flag in command options.

On your own risk, you can specify the additional parameters for the majority of dependency programs
using the '--extra' flag in general EXONtools options (for safety precautions leave
a space before the opening quote, otherwise the EXONtools will raise an error):
EXONtools.py --extra " -argument" COMMAND [command options]

Supplementary EXONtools commands (for advanced users only):
EXONtools.py force_install        - force pip to install all required external python modules
EXONtools.py SetDefaultConfig     - write/overwrite default configuration settings in dependencies.ini
EXONtools.py test                 - run all unit tests

Run ./EXONtools.py -c/--citations to print the complete list of references
including all dependency programs used in the EXONtools pipeline"""
