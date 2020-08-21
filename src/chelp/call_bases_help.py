description = """

                    #########################
                            CALL_BASES
                    #########################

This command analyzes BAM/SAM mapping results, estimates the coverage
of each contig in the reference and identifies its heterozygous sites.
It can also trim contig ends by required threshold of read coverage.
Finally, it produces fasta file with called bases"""

epilog = """
Provide path to a folder with BAM/SAM files (-i) or to a single BAM/SAM
file containing reads that have been mapped on the target sequences (-r).
Currently, the following file extensions are supported for FASTQ files
in the mapping procedure: ['.bam','.sam'].

If you are NOT using BAM/SAM files produced by 'map_reads' command in EXONtools,
you must sort them manually in samtools before running the command 'call_bases'.

Provide path to the folder with FASTA files containing reference (target)
sequences (-r) or path that leads directly to a single FASTA file (-r).
These reference files should be exactly the same as were used for mapping.
Currently, the following extensions are supported for FASTA files with
reference sequences: ['.fasta','.fa'].

IMPORTANT: Please verify that all file names (both BAM/SAM and FASTA files) include
library names before the first underscore (e.g., libraryname_paired_R1.fq),
which will be used for automatic identification of which file belongs
to which library. It means that you must rename the reference file if you are
using 'call_bases' command for running the reference-based assembly.

IMPORTANT: internal bases in reference contigs having zero read coverage
(i.e., coverage gaps) will be annotated as N's by default. However, you can select
the option '--nogaps' to output only the largest segment of the contig.
Ideally, all contigs produced by the EXONtools pipeline should not produce
any coverage gaps. But they might appear if you are using 'map_reads' command
for creating a reference-based assembly.

You may choose to specify the coverage limit that will be used to trim
flanking regions in each contig if those will be below the threshold (-c option).

While doing the base calling procedure, you may filter reads which will be treated
in the base calling analysis (leave the default settings if you are not sure
what all these parameters are used for)
    - maximum alignment error rate (default: 0.05) = NM/TRIMMED_READ_LENGTH
    - minimum mapping quality (default: 0) - corresponds to uniqueness of mapping
    - minimum read length (default: 50 bp)
    - number of bases from each read flanking region that will be judged
        more strictly than other read bases (default: 5, use 0 to turn off this option)
    - minimum number of insertion occurrences across all mapped reads to
        begin the checkup if it is true, i.e. exceeds 50% (default: 5)

Select '--noiupac' option if you would like to retain the original
(reference or random if not present in reference) nucleotide bases in all heterozygote sites

Use '-S/--stats' general option to create 'basecalling_stats.csv' file
containing all information about contig processing results.

I will be happy to receive your feedback on how this pipeline step can be improved.

Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues
or shoot me an email to <vinni(at)hawaii.edu> with your suggestions
"""
