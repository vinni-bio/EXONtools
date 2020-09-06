# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from __future__ import print_function, division
import os
import re
import logging
import gzip
from math import log as LN

from mains.EXT_errors import EXONtoolsError
from utils.sorting import natural_sort
from utils.mapq import MAPQ


class SeqIO(object):
    """Class for sequence file processing and storage"""

    extfasta = [".fasta", ".fas", ".fa", ".fna", ".prot", ".pep"]
    extfastq = [".fastq", ".fq"]
    extsam = [".sam"]
    custom = []
    dryrun = False
    exttype = None

    def __init__(self, path, fileformat="FASTA", strict=True):
        formdict = {"fasta": SeqIO.extfasta, "fastq": SeqIO.extfastq, "sam": SeqIO.extsam, "custom": SeqIO.custom}
        try:
            if SeqIO.exttype:
                self.type = SeqIO.exttype.lower()
            else:
                self.type = fileformat.lower()
            self.extensions = formdict[fileformat.lower()]
        except KeyError:
            logging.error("SeqIO error. Invalid file format")
            raise EXONtoolsError("SeqIO error")

        if os.path.exists(path) and os.path.isfile(path):
            zcheck = os.path.splitext(path)
            if (zcheck[1] in self.extensions) or (zcheck[1] == ".gz" and os.path.splitext(zcheck[0])[1] in self.extensions):
                self.path = os.path.realpath(path)
                self.SEQS = None
                self.strict = strict
                self.total = 0
            else:
                logging.error("Your file extension is not supported. Please use: '{0:s}'".format("', '".join(self.extensions)))
                raise EXONtoolsError("SeqIO error")
        elif SeqIO.dryrun:
            self.path = ""
        else:
            logging.error("Provided file path does not exist")
            logging.error(path)
            raise EXONtoolsError("SeqIO error")

    @classmethod
    def setdry(cls):
        cls.dryrun = True

    def delete(self):

        logging.debug("Deleting the file '{0:s}'".format(self.path))
        try:
            if not SeqIO.dryrun:
                os.remove(self.path)
            logging.debug("The file '{0:s}' is deleted: OK".format(os.path.basename(self.path)))
        except OSError:
            logging.error("EXONtools failed to delete '{0:s}'".format(os.path.basename(self.path)))
            raise EXONtoolsError("Cannot delete SeqIO file instance")
        except TypeError:
            logging.error("Cannot delete '{0:s}' file because it has not been assigned yet".format(os.path.basename(self.path)))
            raise EXONtoolsError("Cannot delete file from non-existing path")

    def totalcount(self, lines=False):
        """count records"""

        if not SeqIO.dryrun and os.path.getsize(self.path):

            try:
                if self.path.endswith(".gz"):
                    infile = gzip.open(self.path, 'rt')
                else:
                    infile = open(self.path, 'r')
            except IOError:
                logging.error("Cannot open input file {0:s}. Please verify its format.".format(os.base.name(self.path)))
                raise EXONtoolsError("SeqIO error")

            if lines:
                self.total = sum(bl.count('\n') for bl in readblock(infile))
                infile.close()

            else:
                line = infile.readline().strip()
                while(line and line == ""):
                    line = infile.readline().strip()

                if self.type == "sam" and line:
                    while(line):
                        if not line.startswith("@"):            # skip all headers if present
                            self.total += 1
                        line = infile.readline().strip()
                    infile.close()

                elif self.type == "fasta" and line:
                    if line.startswith(">"):
                        while(line):
                            if line.startswith(">"):
                                self.total += 1
                                line = infile.readline().strip()
                            else:
                                line = infile.readline().strip()
                    else:
                        logging.error("The first line in FASTA file should start with '>' sign")
                        raise EXONtoolsError("SeqIO error")
                    infile.close()

                elif self.type == "fastq" and line:
                    if line.startswith("@"):
                        while(line):
                            if line.startswith("@"):
                                self.total += 1
                                infile.readline()
                                infile.readline()
                                line = infile.readline().strip()
                            else:
                                line = infile.readline().strip()
                    else:
                        logging.error("The first line in FASTQ file should start with '@' sign")
                        raise EXONtoolsError("SeqIO error")
                    infile.close()

                elif not line:
                    pass

                else:
                    logging.error("The format for SeqIO operation is not defined. This is probably a program bug")
                    raise EXONtoolsError("SeqIO error")

    def read(self):
        """Read file records"""

        if not SeqIO.dryrun and os.path.getsize(self.path):
            self.total = 0
            try:
                if self.path.endswith(".gz"):
                    infile = gzip.open(self.path, 'rt')
                else:
                    infile = open(self.path, 'r')
            except IOError:
                logging.error("Cannot open input file {0:s}. Please verify its format.".format(os.base.name(self.path)))
                raise EXONtoolsError("SeqIO error")

            line = infile.readline().strip()
            while(line and line == ""):
                line = infile.readline().strip()

            if self.type == "sam" and line:
                while(line):
                    if not line.startswith("@"):            # skip all headers if present
                        self.total += 1
                        yield samline(line)
                    line = infile.readline().strip()
                infile.close()

            elif self.type == "fasta" and line:
                if line.startswith(">"):
                    while(line):
                        if line.startswith(">"):
                            seqname = line[1:]
                            # seqname = line[1:].replace(" ", "_").split()[0]
                            seqline = ""
                            line = infile.readline().strip()
                            while(line):
                                seqline = seqline + line
                                line = infile.readline().strip()
                                if line is None or line.startswith(">"):
                                    break
                            self.total += 1
                            yield fastaseq(seqname, seqline.upper())
                        else:
                            line = infile.readline().strip()
                else:
                    logging.error("The first line in FASTA file should start with '>' sign")
                    raise EXONtoolsError("SeqIO error")
                infile.close()

            elif self.type == "fastq" and line:
                if line.startswith("@"):
                    while(line):
                        if line.startswith("@"):
                            paired = 0
                            identifierline = line[1:]
                            nameline = identifierline.split()
                            seqname = nameline[0]
                            barcode = None
                            filtered = False
                            if re.search("/1|/2", seqname[-2:]):
                                paired = int(seqname[-1])
                                seqname = seqname[:-2]
                            elif re.search("_R1|_R2", seqname):
                                paired = int(re.search("_R1|_R2", seqname).group()[-1])
                            elif re.search("\d:\w:\d:.*", nameline[1]):
                                xlist = nameline[1].split(":")
                                paired = int(xlist[0])
                                barcode = xlist[3]
                                if xlist[1] == "Y":
                                    filtered = True
                            else:
                                pass
                            if len(nameline) > 1:
                                extra = "_".join(nameline[1:])
                            else:
                                extra = ""
                            seqline = infile.readline().strip().upper()
                            infoline = infile.readline().strip()[1:]
                            qualline = infile.readline().strip()
                            line = infile.readline().strip()
                            self.total += 1
                            yield fastqseq(identifierline, seqname, seqline, qualline, infoline, paired, extra, filtered, barcode)
                        else:
                            line = infile.readline().strip()
                else:
                    logging.error("The first line in FASTQ file should start with '@' sign")
                    raise EXONtoolsError("SeqIO error")
                infile.close()
            elif not line:
                pass
            else:
                logging.error("The format for SeqIO operation is not defined. This is probably a program bug")
                raise EXONtoolsError("SeqIO error")

    def readall(self):
        self.seqnames = []
        self.SEQS = {}
        if self.type == "sam":
            self.strict = False
        if not SeqIO.dryrun and os.path.getsize(self.path):
            for seq in self.read():
                if self.strict and (set([seq.name]) & set(self.seqnames)):
                    logging.error("Two similar sequence names are found in one file: ")
                    logging.error("{0:s} in {1:s}".format(seq.name, os.path.basename(self.path)))
                    raise EXONtoolsError("SeqIO error")
                else:
                    self.seqnames.append(seq.name)
                    try:
                        self.SEQS[seq.name].append(seq)
                    except KeyError:
                        self.SEQS[seq.name] = [seq]
            if self.strict and len(self.SEQS) == 0:
                logging.error("No sequences were identified. Please verify that your file has the right format")
                raise EXONtoolsError("SeqIO error")

    def callread(self, seqname):
        if self.SEQS:
            try:
                return self.SEQS[seqname]
            except KeyError:
                pass

    @staticmethod
    def writefasta(seqs, filepath, sortseq=True, overwrite=True):
        if not SeqIO.dryrun and isinstance(seqs, dict):
            if overwrite:
                mode = 'w'
            else:
                mode = 'a'
            with open(filepath, mode) as outfile:
                try:
                    if sortseq:
                        for seq in sorted(seqs.keys(), key=natural_sort):
                            outfile.write(">{0:s}\n{1:s}\n".format(seq, seqs[seq]))
                    else:
                        for seq in seqs.keys():
                            outfile.write(">{0:s}\n{1:s}\n".format(seq, seqs[seq]))
                except (TypeError, KeyError):
                    logging.error("You must provide dictionary with sequence names as keys and sequences as values")
                    raise EXONtoolsError("SeqIO write fasta method error")
        else:
            logging.error("You must provide dictionary with sequence names as keys and sequences as values")
            raise EXONtoolsError("SeqIO write fasta method error")

    @classmethod
    def setformat(cls, extensions=None, exttype="FASTA"):
        """Setting the file extension requirements"""

        cls.exttype = exttype
        if isinstance(extensions, (list,)):
            for ext in extensions:
                if not isinstance(ext, (str,)):
                    logging.error("Unknown format is provided for file extension")
                    raise EXONtoolsError("SeqIO error")
            cls.custom = extensions
        elif extensions == "default":
            cls.custom = [".fasta", ".fas", ".fa", ".prot", ".fastq", ".fq", ".sam"]
        elif isinstance(extensions, (str,)):
            cls.custom = [extensions]
        else:
            logging.error("Unknown format is provided for SeqIO file extension")
            raise EXONtoolsError("SeqIO error")

    def nostrict(self):
        self.strict = False


class fastaseq(object):
    """Class for each read within FASTA file"""

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq.upper()


class fastqseq(object):
    """Class for each read within FASTQ file"""

    def __init__(self, identifier, name, seq, qual, info, pair, extra, filtered=False, barcode=None):
        self.identifier = identifier
        self.name = name
        self.seq = seq.upper()
        self.qual = qual
        self.info = info
        self.pair = pair
        self.extra = extra
        self.filtered = filtered
        self.barcode = barcode


class samline(object):
    """Class for parsing each read line within SAM file"""

    def __init__(self, line):

        samlist = line.split('\t')
        samtags = samlist[11:]
        XAZ = [x[5:].split(';') for x in samtags if x.startswith("XA:Z:")]
        SAZ = [x[5:].split(';') for x in samtags if x.startswith("SA:Z:")]
        AS = [int(x[5:]) for x in samtags if x.startswith("AS:i:")]
        XS = [int(x[5:]) for x in samtags if x.startswith("XS:i:")]
        NM = [int(x[5:]) for x in samtags if x.startswith("NM:i:")]

        """Initial class attributes and procedures"""
        self.name = samlist[0]                              # read name
        self.target = samlist[2]                            # target name
        self.seq = samlist[9].upper()                       # read nucleotide sequence
        self.qual = samlist[10]                             # read quality
        self.start = int(samlist[3]) - 1                    # adjusted mapping start position on the target
        self.mapq = int(samlist[4])                         # read mapping quality
        self.cigar = samlist[5].replace("*", "")            # CIGAR sequence
        self.flag = int(samlist[1])                         # SAM flag
        self.insert = 0                                     # insert count initially equals 0
        self.readlen = len(self.seq)                        # intial read length
        self.insertions = {}                                # dictionary with insertion seqs and their coordinates
        self.duplicate = False                              # if read is a PCR duplicate
        self.secondary = False                              # if read is a secondary mapping
        self.supplementary = False                          # if read is chimeric
        self.suppls = []                                    # list of supplementary mappings
        self.secs = []                                      # list of secondary mappings
        self.AS = None                                      # alignment score
        self.XS = None                                      # suboptimal alignment score
        self.NM = None                                      # number of mismatches
        self.type = None                                    # read type [0,1,2]
        self.reallen = None                                 # adjusted length after cigar processing
        self.sister = None                                  # secondary target

        # PARSE CIGAR STRING
        if self.cigar:
            slen = ""                                           # section length
            check_start = True                                  # switcher to show the start position
            curpos = 0                                          # current position on the read
            for i, s in enumerate(self.cigar):                  # for each CIGAR block
                if s.isdigit():                                 # add number of bases to slen if digit
                    slen = slen + s
                else:
                    if s == "S" and len(self.cigar) == i + 1:   # if skip (S) at the end
                        self.seq = self.seq[:-int(slen)]        # trim end of nucleotide sequence
                        self.qual = self.qual[:-int(slen)]      # trim end of quality sequence
                    elif s == "S" and check_start:              # if skip (S) at the beginning
                        self.seq = self.seq[int(slen):]         # trim the start of nucleotide sequence
                        self.qual = self.qual[int(slen):]       # trim the start of quality sequence
                        check_start = False                     # turn off the switcher
                    elif s == "I":                              # if insertion (I)
                        self.seq = self.seq[:curpos] + self.seq[curpos + int(slen):]     # delete insertion in nuc seq
                        self.qual = self.qual[:curpos] + self.qual[curpos + int(slen):]  # delete insertion in qual seq
                        self.insertions.update({self.start + curpos: (self.seq[curpos + 1:curpos + 1 + int(slen)], self.qual[curpos + 1:curpos + 1 + int(slen)])})
                    elif s == "D":                              # if deletion (D)
                        self.seq = self.seq[:curpos] + "-" * int(slen) + self.seq[curpos:]    # add gaps '-' to nuc seq
                        self.qual = self.qual[:curpos] + "$" * int(slen) + self.qual[curpos:]  # add '$' to qual seq
                        curpos += int(slen)                     # adjust the current position
                    else:
                        curpos += int(slen)                     # adjust the current position
                    slen = ""                                   # reset section length

        self.reallen = len(self.seq)                            # adjusted length of the processed seq
        self.end = self.start + self.reallen                    # mapping end position on the target

        # Identify type of the read
        if self.flag & 0x40:
            self.type = 1
        elif self.flag & 0x80:
            self.type = 2
        else:
            self.type = 0

        # Check if duplicate
        if self.flag & 0x400:
            self.duplicate = True

        # Check if secondary
        if self.flag & 0x100:
            self.secondary = True
        elif XAZ:
            self.secs = [tuple(x.split(',')) for x in XAZ[0]]
        else:
            pass

        # Check if supplementary
        if self.flag & 0x800:
            self.supplementary = True
        elif SAZ:
            self.suppls = [tuple(x.split(',')) for x in SAZ[0]]
        else:
            pass

        # check insert size
        insertcheck = int(samlist[8])
        if insertcheck < 0:
            self.insert = abs(insertcheck)

        nextmap = samlist[6].replace("*", "").replace("=", "").strip()
        if nextmap:
            self.sister = nextmap

        if AS:
            self.AS = abs(AS[0])
        if NM:
            self.NM = abs(NM[0])
        if XS:
            self.XS = abs(XS[0])

        if self.cigar and self.mapq == 255:
            if self.flag & 0x2:
                self.mapq = 42
            else:
                self.mapq = MAPQ(AS=self.AS, XS=self.XS, seqlen=self.reallen)

        if self.cigar and self.AS <= 0:
            if abs(self.AS) < self.reallen:
                self.AS = self.reallen - abs(self.AS)
            else:
                self.AS = 0

        if self.XS and self.XS < 0:
            self.XS = self.reallen - abs(self.XS)

        if self.NM:
            self.alqual = int(-10 * LN(self.NM / self.reallen))
        elif self.reallen:
            self.alqual = int(-10 * LN(0.5 / self.reallen))
        else:
            logging.error("Unknown alignment quality for '{0:s}' read with {1:d} flag".format(self.name, self.flag))
            raise EXONtoolsError("Read addition error")
        if self.alqual < 0:
            self.alqual = 0


def readblock(infile, size=65536):
    """Read file block"""

    while True:
        b = infile.read(size)
        if not b:
            break
        yield b
