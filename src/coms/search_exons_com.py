# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.exosearch import exosearcher


class command(EXTcommand):
    """This command takes annotated assembly (or multiple assemblies)
    and predicts exon boundaries within each contig based on provided reference"""

    supported_programs = [exosearcher]
    default_program = exosearcher
    command_name = "search_exons"
