# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.bbduk_cleaner import bbduk
from progs.cutadapt_cleaner import cutadapt
from progs.trimmomatic_cleaner import trimmomatic


class command(EXTcommand):
    """This command deletes reads with low average quality and reads,
    which size is lower than the required minimum. It also trims
    adapter residuals and bases with low quality"""

    supported_programs = [bbduk, cutadapt, trimmomatic]
    default_program = trimmomatic
    command_name = "clean_reads"
