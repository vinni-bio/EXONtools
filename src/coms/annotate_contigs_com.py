# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from mains.EXT_com import EXTcommand
from progs.annotator import annotator


class command(EXTcommand):
    """This command takes assembly (or multiple assemblies)
    and annotates its contigs based on provided reference"""

    supported_programs = [annotator]
    default_program = annotator
    command_name = "annotate_contigs"
