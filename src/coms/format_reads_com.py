# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.formatter import seqformatter


class command(EXTcommand):
    """This command takes raw reads and verifies their names and structure"""

    supported_programs = [seqformatter]
    default_program = seqformatter
    command_name = "format_reads"
