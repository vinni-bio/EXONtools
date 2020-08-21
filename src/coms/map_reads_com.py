# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.bwa_mapper import bwa_mapper
from progs.bwt2_mapper import bwt2_mapper


class command(EXTcommand):
    """This command maps reads to a reference genome or to assembly,
estimates the coverage of each contig and identifies heterozygous sites.
It can also trim contig ends by required threshold of read coverage.
Finally, it produces fasta file with called bases from read mapping
results"""

    supported_programs = [bwa_mapper, bwt2_mapper]
    default_program = bwa_mapper
    command_name = "map_reads"
