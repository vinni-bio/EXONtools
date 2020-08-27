# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.


from __future__ import print_function, division
import os
import logging

from mains.EXT_errors import EXONtoolsError


class BLAST6(object):
    """Class for parsing blast output files and storage their results"""

    dryrun = False

    def __init__(self, path, fileformat=".outfmt6"):

        if os.path.exists(path) and os.path.isfile(path):
            if path.endswith(fileformat):
                self.path = os.path.realpath(path)
            else:
                logging.error("Your BLAST OUTPUT file extension is not supported. Please use: '.outfmt6'")
                raise EXONtoolsError("SeqIO error")
        elif BLAST6.dryrun:
            self.path = ""
        else:
            logging.error("Provided file path does not exist")
            raise EXONtoolsError("BLAST6 error")

    @classmethod
    def setdry(cls):
        cls.dryrun = True

    def read(self):
        if not BLAST6.dryrun and os.path.getsize(self.path):
            with open(self.path, 'r') as infile:
                for line in infile:
                    yield blastline(line.strip())

    def delete(self):
        logging.debug("Deleting the file '{0:s}'".format(self.path))
        try:
            if not BLAST6.dryrun:
                os.remove(self.path)
            logging.debug("The file '{0:s}' is deleted: OK".format(os.path.basename(self.path)))
        except OSError:
            logging.error("EXONtools failed to delete '{0:s}'".format(os.path.basename(self.path)))
            raise EXONtoolsError("Cannot delete BLAST6 file instance")
        except TypeError:
            logging.error("Cannot delete '{0:s}' file because it has not been assigned yet".format(os.path.basename(self.path)))
            raise EXONtoolsError("Cannot delete file from non-existing path")


class blastline(object):
    """Represents one line in blast output (Table format)"""

    source = "EXONtools"
    blastprog = None
    dryrun = False
    totalnum = 0

    def __init__(self, line):
        blastline.totalnum += 1
        self.id = "line" + str(blastline.totalnum)
        matchlist = line.split("\t")
        self.query = matchlist[0]
        self.target = matchlist[1]
        self.identity = round(float(matchlist[2]) / 100, 5)
        self.length = int(matchlist[3])
        self.mismatch = int(matchlist[4])
        self.ngaps = int(matchlist[5])
        self.qstart = int(matchlist[6])
        self.qend = int(matchlist[7])
        self.tstart = int(matchlist[8])
        self.tend = int(matchlist[9])
        self.evalue = float(matchlist[10])
        self.bit = float(matchlist[11])
        self.scaffolds = []
        self.isoforms = []
        self.chimeric = False

    @classmethod
    def setblast(cls, prog):
        cls.blastprog = prog

    @classmethod
    def setsource(cls, name):
        cls.source = name

    def missasembled(self):
        self.chimeric = True

    def add_iso(self, name):
        self.isoforms.append(name)

    def add_scaf(self, name):
        self.scaffolds.append(name)
