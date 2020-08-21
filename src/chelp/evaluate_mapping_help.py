# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """

                    ############################
                          EVALUATE_MAPPING
                    ############################

This command evaluates various metrics for read mapping quality"""

epilog = """
Provide path to the folder with SAM/BAM files or path to a single SAM/BAM file (-i).
Currently, the following file extensions are supported ['.sam','.bam'].

The output file 'mapping_stats.csv' will be saved to the current directory by default.
Alternatively, you can use '-o/--out" option to define the output path that you prefer.

Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues
or shoot me an email to <vinni(at)hawaii.edu> with your suggestions
"""
