# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

description = """

                    ############################
                          EVALUATE_ASSEMBLY
                    ############################

This command evaluates various metrics for contigency quality of assemblies"""

epilog="""
Provide path to the folder with fasta files containing assemblies or path
to a single assembly file (-i). Currently, the following file extensions
are supported for assembly files ['.fasta','.fa'].

The output file 'assembly_stats.csv' will be saved to the current directory by default.
Alternatively, you can use '-o/--out" option to define the output path that you prefer.

Please report all bugs here: https://github.com/vinni-bio/EXONtools/issues
or shoot me an email to <vinni(at)hawaii.edu> with your suggestions
"""
