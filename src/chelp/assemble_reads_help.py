# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """

                        ##########################
                              ASSEMBLE_READS
                        ##########################

This command assembles sequencing reads stored in FASTQ files for each library.
"""

epilog = """
Provide path to the folder with fastq files containing sequencing reads
that need to be assembled (-i) or to single paired and/or unpaired
FASTQ files (-R1, -R2, -U). Currently, the following file extensions
are supported in the assembling procedure ['.fastq','.fq'].

Please verify that all file names include library names before the first
underscore (e.g., libraryname_paired_R1.fq), which will be used for
automatic identification of files belonging to each library. Paired files
must include '_R1' and '_R2' identifiers in their names, ,otherwise they
will be treated as unpaired files.

ABySS and Trans-ABySS assemblers usually need multiple k-mer lengths to
reconstruct multiple assemblies per each library. However, Trinity
assembler usually performs very well (but very slow) with its default k-mer
length. Therefore, use the '--kmers' command option to provide the list
of kmer lengths for ABySS programs only.

IMPORTANT!!!
ABySS and TransABySS programs can run assembly analysis in parallel mode (MPI).
If you have the Open-MPI installed in your computer (cluster),
please provide the MPI path using '--mpirun' option. Then you can set
the number number of required nodes (-p) and the number of cores per each node (-j).
Please note that the EXONtools will be using a single machine with a number of
threads (-T/--threads) divided by the number of cores assigned to each job (-j)
if the --mpirun path is not provided.

Some MPIs require the additional environment variables set before
the command can be run. Here is an example on how I set those on my HPC:
export OMP_NUM_THREADS=1
export I_MPI_FABRICS=tmi
export I_MPI_PMI_LIBRARY=./lib64/libpmi.so

This is an example of the EXONtools command with MPI run on 10 machines with 20 cores in each:
./EXONtools.py assemble_reads -i ./INPUT -o ./OUTPUT --mpirun ./openMPI/bin/1.10.2/bin/mpirun
-k 22 32 42 52 62 72 82 92 102 -n 10 -p 20 -j 20 --program abyss

Use '-S/--stats' general option to create 'assembly_stats.csv' file containing all
contiguity information for each assembly. The same output can be retrieved later
by using the EXONtools command 'evaluate_assembly'.

Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues
or shoot me an email to <vinni(at)hawaii.edu> with your suggestions
"""
