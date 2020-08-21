# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.


from __future__ import print_function, division
import logging
import os
import json
from random import choice as choose

from mains.EXT_errors import EXONtoolsError
from utils.revcomp import DNArevcomp


def loadcodes(gcode):
    """Load the dictionary with the appropriate genetic codes"""

    coding_schemes = {
        1: 'Standard',
        2: 'VertMt',
        3: 'YeastMt',
        4: 'MoldProMt',
        5: 'InvertMt',
        6: 'CilDasHex',
        9: 'EchFlatMt',
        10: 'Euplotid',
        11: 'BactPlant',
        12: 'YeastNuc',
        13: 'AscidianMt',
        14: 'AltFlatMt',
        15: 'Blepharisma',
        16: 'ChlorMt',
        21: 'TremMt',
        22: 'ScenMt',
        23: 'ThrausMt',
        24: 'PteroMt',
        25: 'SR1Gracil',
        26: 'Pachysolen',
        27: 'Karyorelict',
        28: 'Condylostoma',
        29: 'Mesodinium',
        30: 'Peritrich',
        31: 'Blastocrithidia'
    }

    if gcode not in coding_schemes:
        logging.error("Genetic code ID error. Please choose from {0:s}".format(",".join(map(str, list(coding_schemes.keys())))))
        raise EXONtoolsError("Wrong genetic code ID number")

    jsonpath = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "jsons", "NCBI_GC.json")

    with open(jsonpath, 'r') as gencodefile:
        gene_codes = json.load(gencodefile)["NCBI_GC_4.2"]
    for code in gene_codes.keys():
        for stopcode in gene_codes[code]["stop_codons"]:
            gene_codes[code]["table"][stopcode] = "*"

    return gene_codes[coding_schemes[gcode]]


def clean_iupac(seq):
    """Remove IUPAC encoded ambiguities from heterozygotes"""

    iupac_rev = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "GT": "K", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}

    return "".join([choose(iupac_rev[x]) if x in iupac_rev else x for x in seq.upper()])


def translate(seq, frame, code, finalstop=False):
    """Translate DNA sequence into peptide sequence"""

    frames = [-3, -2, -1, 1, 2, 3]

    if not isinstance(seq, str):
        logging.error("Sequence must be in string format")
        raise EXONtoolsError("Translation sequence type error")

    if frame not in frames:
        logging.error("Please choose one from these {0:s} reading frames".format(",".join(map(str, frames))))
        raise EXONtoolsError("Frame format error")

    if len(seq) < 3:
        return

    AASEQ = ""
    gencode = loadcodes(code)
    codetab = gencode['table']
    corseq = clean_iupac(seq)
    if frame < 0:
        corseq = DNArevcomp(corseq)
    i = abs(frame) - 1
    j = i + 3

    if corseq[i:j] in gencode["start_codons"]:
        AASEQ = AASEQ + "M"
        i += 3
        j += 3

    while(j <= len(seq)):
        try:
            AASEQ = AASEQ + codetab[corseq[i:j]]
            i += 3
            j += 3
        except KeyError:
            AASEQ = AASEQ + "X"
            i += 3
            j += 3
    try:
        if not finalstop and AASEQ[-1] == "*":
            AASEQ = AASEQ[:-1]
    except IndexError:
        pass

    if AASEQ:
        return AASEQ


def findorfs(seq, code, minlen, diff=0.5):
    """Find ORFs within provided sequence"""

    frames = [-3, -2, -1, 1, 2, 3]
    codons = {-3: [], -2: [], -1: [], 1: [], 2: [], 3: []}
    longorfs = {-3: None, -2: None, -1: None, 1: None, 2: None, 3: None}

    gencode = loadcodes(code)

    if not isinstance(seq, str):
        logging.error("Sequence must be in string format")
        raise EXONtoolsError("Translation sequence type error")

    if len(seq) < 3 or len(seq) < minlen:
        return

    forseq = clean_iupac(seq)
    revseq = DNArevcomp(forseq)

    for f in frames:
        i = abs(f) - 1
        j = i + 3
        if f > 0:
            lenforseq = len(forseq)
            while(j <= lenforseq):
                codons[f].append(triplet(forseq[i:j], i + 1, f))
                i += 3
                j += 3
        else:
            lenrevseq = len(revseq)
            while(j <= lenrevseq):
                codons[f].append(triplet(revseq[i:j], lenrevseq - i, f))
                i += 3
                j += 3

    for l in longorfs.keys():
        longest = orfseq()
        newseq = orfseq()
        check = True
        for tripl in codons[l]:
            if tripl.seq in gencode["start_codons"] and check:
                newseq = orfseq()
                newseq.append(tripl)
                check = False
            elif tripl.seq in gencode["stop_codons"]:
                newseq.append(tripl)
                if newseq.length > longest.length and newseq.length >= minlen:
                    longest = newseq
                newseq = orfseq()
                check = True
            else:
                newseq.append(tripl)
        if newseq.length > longest.length and newseq.length >= minlen:
            longest = newseq
        longorfs[l] = longest

    maxlen = max(map(lambda x: x.length, longorfs.values()))

    [longorfs.pop(x) for x in list(longorfs.keys()) if longorfs[x].length <= maxlen * diff]

    return list(longorfs.values())


class triplet(object):
    """Stores triplet sequence and coordinates of its nucleotides"""

    def __init__(self, seq, pos, frame):
        self.seq = seq
        self.pos = pos
        self.frame = frame


class orfseq(object):
    """Stores orf sequence and its coordinates"""

    def __init__(self):
        self.seq = ""
        self.start = None
        self.end = None
        self.frame = None
        self.length = 0

    def append(self, tripl):
        """Appends triplet to ORF sequence"""

        self.seq += tripl.seq

        if self.start is None:
            if tripl.frame > 0:
                self.start = tripl.pos
                self.length = 3
                self.frame = tripl.frame
                self.end = tripl.pos + 2
            else:
                self.start = tripl.pos - 2
                self.length = 3
                self.frame = tripl.frame
                self.end = tripl.pos
        else:
            if self.frame > 0:
                self.length += 3
                self.end += 3
            else:
                self.length += 3
                self.start -= 3
