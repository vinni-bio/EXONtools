# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.


from __future__ import print_function
import argparse
import sys


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """SETTING THE OUTPUT CLASS FOR ARGPARSE HELP MENUS"""
    pass


class CITATION(argparse.Action):
    """SETTING THE ACTION CLASS FOR CITATION OPTION"""

    citations = """

Please provide this reference in your publications if you were using the EXONtools pipeline:

'Vinnikov, K.A., 2018. EXONtools: a complete pipeline for exon capture sequencing data
analysis of non-model organisms. Version 0.2b. URL: https://github.com/vinni-bio/EXONtools'

The paper describing the EXONtools pipeline has been submitted to Molecular Ecology Resources
journal in August, 2018. So, fingers crossed. Please follow up the EXONtools github page
to receive the most recent updates.

-----------------------------------------------------------------------------
LIST OF REFERENCES FOR ALL DEPENDENCY PROGRAMS USED IN THE EXONTOOLS PIPELINE
-----------------------------------------------------------------------------

ABySS
'Simpson, J.T., Wong, K., Jackman, S.D., Schein, J.E., Jones, S.J. and Birol, I., 2009.
ABySS: a parallel assembler for short read sequence data. Genome Res. 19, 1117-1123'

BayesHammer
'Nikolenko, S.I., Korobeynikov, A.I. and Alekseyev, M.A., 2013. BayesHammer:
Bayesian clustering for error correction in single-cell sequencing. BMC Genomics 14, S7'

BBmap
'Bushnell, B., 2018. BBMap: Short read aligner for DNA and RNA-seq data.
URL: https://sourceforge.net/projects/bbmap/ (accessed 19 July, 2018)'

BLAST
'Altschul, S.F., Gish, W., Miller, W., Myers, E.W. and Lipman, D.J., 1990.
Basic local alignment search tool. J. Mol. Biol. 215, 403-410'

BLAT
'Kent, W.J., 2002. BLAT - the BLAST-like alignment tool. Genome Res. 12, 656-664'

Bowtie2
'Langmead, B. and Salzberg, S.L., 2012. Fast gapped-read alignment with Bowtie 2.
Nat. Methods 9, 357-359.'

bwa
'Li, H. and Durbin, R., 2009. Fast and accurate short read alignment with
Burrows-Wheeler transform. Bioinformatics 25, 1754-1760'

CAP3
'Huang, X. and Madan, A., 1999. CAP3: A DNA sequence assembly program.
Genome Res. 9, 868-877'

cd-hit
'Li, W. and Godzik, A., 2006. Cd-hit: a fast program for clustering and comparing
large sets of protein or nucleotide sequences. Bioinformatics 22, 1658-1659'

Cutadapt
'Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput
sequencing reads. EMBnet J. 17, 10-12'

FastQC
'Andrews, S., 2010. FastQC: a quality control tool for high throughput sequence data.'

MAFFT
'Katoh, K. and Standley, D.M., 2013. MAFFT multiple sequence alignment software version 7:
improvements in performance and usability. Mol. Biol. Evol. 30, 772-780'

MultiQC
'Ewels, P., Magnusson, M., Lundin, S. and Käller, M., 2016. MultiQC: summarize analysis
results for multiple tools and samples in a single report. Bioinformatics 32, 3047-3048.'

PEAR
'Zhang, J., Kobert, K., Flouri, T. and Stamatakis, A., 2013. PEAR: a fast and accurate
Illumina Paired-End reAd mergeR. Bioinformatics 30, 614-620'

SAMTOOLS
'Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G.,
Abecasis, G. and Durbin, R., 2009. The sequence alignment/map format and SAMtools.
Bioinformatics 25, 2078-2079'

SPAdes
'Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S.,
Lesin, V.M., Nikolenko, S.I., Pham, S., Prjibelski, A.D. and Pyshkin, A.V., 2012.
SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing.
J. Comput. Biol. 19, 455-477'

Trans-ABySS
Robertson, G., Schein, J., Chiu, R., Corbett, R., Field, M., Jackman, S.D., Mungall, K.,
Lee, S., Okada, H.M., Qian, J.Q. and Griffith, M., 2010. De novo assembly and analysis
of RNA-seq data. Nat. Methods 7, 909-912.'

Trimmomatic
'Bolger, A.M., Lohse, M. and Usadel, B., 2014. Trimmomatic: a flexible trimmer for
Illumina sequence data. Bioinformatics 30, 2114-2120'

Trinity
'Grabherr, M.G., Haas, B.J., Yassour, M., Levin, J.Z., Thompson, D.A., Amit, I.,
Adiconis, X., Fan, L., Raychowdhury, R., Zeng, Q. and Chen, Z., 2011.
Full-length transcriptome assembly from RNA-Seq data without a reference genome.
Nat. Biotechnol. 29, 644-652'

"""

    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help="Show the list of references for all dependency programs currently used in the EXONtools and exit"):
        super(CITATION, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        formatter = parser._get_formatter()
        formatter.add_text(CITATION.citations)
        parser._print_message(formatter.format_help(), sys.stdout)
        parser.exit()
