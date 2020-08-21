# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.basecaller import basecaller


class command(EXTcommand):
    """This command analyzes BAM/SAM mapping results, estimates
the coverage of each contig and identifies its heterozygous sites.
It can also trim contig ends by required threshold of read coverage.
Finally, it produces fasta file with called bases."""

    supported_programs = [basecaller]
    default_program = basecaller
    command_name = "call_bases"
