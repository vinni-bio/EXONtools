# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.demultiplexer import demultiplexer


class command(EXTcommand):
    """This command splits file with reads by library names or indexes"""

    supported_programs = [demultiplexer]
    default_program = demultiplexer
    command_name = "demultiplex_reads"
