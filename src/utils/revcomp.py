# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.


from __future__ import print_function, division
import logging
from mains.EXT_errors import EXONtoolsError


def DNArevcomp(seq):
    """Returns the reversed complementary DNA sequence"""

    complements = {'A': 'T', 'T': 'A', 'U': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M', 'W': 'W', 'S': 'S', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}

    if not isinstance(seq, str):
        logging.error("Sequence must be in string format")
        raise EXONtoolsError("DNArevcomp sequence type error")

    try:
        return "".join([complements[base] for base in reversed(seq.upper())])
    except KeyError:
        logging.error("Only the following bases can be processed: 'ABCDHKMNRSUVWY'")
        logging.error("Please check your sequence for any unacceptable characters")
        raise EXONtoolsError("DNA reverse complement function error")
