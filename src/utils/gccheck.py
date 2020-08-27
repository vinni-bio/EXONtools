# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from __future__ import print_function, division
import logging
import re

from mains.EXT_errors import EXONtoolsError


def GC(seq, iupac=False):
    """Checks GC content for the nucleotide sequence skipping all ambiguity codes"""

    if isinstance(seq, str):
        AT = len([x.start() for x in re.finditer('A|T', seq.upper())])
        CG = len([x.start() for x in re.finditer('C|G', seq.upper())])
        total = AT + CG
        if iupac:
            S = len([x.start() for x in re.finditer('R|Y|K|M', seq.upper())])
            RYKM = len([x.start() for x in re.finditer('R|Y|K|M', seq.upper())])
            DH = len([x.start() for x in re.finditer('D|H', seq.upper())])
            BV = len([x.start() for x in re.finditer('B|V', seq.upper())])
            W = len([x.start() for x in re.finditer('B|V', seq.upper())])
            total = total + S + RYKM + DH + BV + W
            CG = CG + S + RYKM * 0.5 + DH * 0.33 + BV * 0.66
        try:
            return round(CG / total, 2)
        except ZeroDivisionError:
            logging.error("No standard nucleotide codes are found in the provided sequence")
            raise EXONtoolsError("Sequence format error")
    else:
        logging.error("Sequence must have a string type")
        raise EXONtoolsError


def Ncheck(seq):
    """Checks N bases in the nucleotide sequence"""

    if isinstance(seq, str):
        return len([x.start() for x in re.finditer('N|n', seq)])
    else:
        logging.error("Sequence must have a string type")
        raise EXONtoolsError
