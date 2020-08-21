# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.assasses import assass


class command(EXTcommand):
    """This class estimates various metrics of assembly quality"""

    supported_programs = [assass]
    default_program = assass
    command_name = "evaluate_assembly"
