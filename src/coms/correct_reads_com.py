# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.hammer import hammer


class command(EXTcommand):
    """This command finds and fixes sequencing errors within fastq files"""

    supported_programs = [hammer]
    default_program = hammer
    command_name = "correct_reads"
