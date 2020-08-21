# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.abyss import abyss
from progs.spades import spades
from progs.transabyss import transabyss
from progs.trinity import trinity


class command(EXTcommand):
    """This class assembles sequencing reads with different programs"""

    supported_programs = [abyss, spades, transabyss, trinity]
    default_program = abyss
    command_name = "assemble_reads"
