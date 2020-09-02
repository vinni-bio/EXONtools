# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from __future__ import print_function
import logging
import os
import shutil
import re

from mains.EXT_errors import EXONtoolsError
from utils.sorting import natural_sort


class output(object):
    """The class that administrates the output directory"""

    path = None
    dryrun = False
    strict = True

    def __init__(self, outpath):
        output.path = output.test_outpath(outpath)

    def __repr__(self):
        return output.path

    @classmethod
    def setdry(cls):
        cls.dryrun = True

    @classmethod
    def nostrict(cls):
        cls.strict = False

    @classmethod
    def delete(cls):
        logging.info("Deleting the output directory {0:s}".format(output.path))
        try:
            if not output.dryrun:
                shutil.rmtree(output.path)
            logging.debug("The output directory was deleted: OK")
        except OSError as e:
            if e.errno == 2:
                pass
            else:
                logging.error("Something wrong happened when EXONtools was trying to delete the output directory")
                raise EXONtoolsError("Directory cannot be deleted")
        except TypeError:
            logging.error("Cannot delete the output directory because it has not been assigned yet")
            raise EXONtoolsError("Cannot delete directory from non-existing path")

    @staticmethod
    def test_outpath(outpath):
        logging.debug("Testing outpath parameters")
        outpath_list = [x for x in os.path.split(outpath) if x]
        if os.path.exists(outpath) and os.path.isdir(outpath):
            outdir = os.path.abspath(outpath)
            if len(os.listdir(outdir)) > 0 and output.strict:
                logging.error("The output directory exists and is not empty. Please provide another path or clean the directory\n")
                raise EXONtoolsError("The output directory exists and is not empty.")
            else:
                logging.debug("The existing output directory will be used: OK")
                return outdir
        elif len(outpath_list) == 1 and outpath_list[0][0].isalpha():
            outdir = os.path.join(os.path.realpath(os.curdir), outpath_list[0])
        elif len(outpath_list) > 1:
            outdir = os.path.realpath(outpath)
        else:
            logging.error("Something wrong happened when EXONtools was trying to create the output directory. Try another path")
            raise EXONtoolsError("Cannot create the output directory")
        try:
            logging.info("Creating the output directory {0:s}".format(outdir))
            if not output.dryrun:
                os.mkdir(outdir)
            logging.debug("The output directory is created: OK")
            return outdir
        except OSError as e:
            if e.errno == 17:
                logging.error("The output directory exists and is not empty. Please provide another path.")
                raise EXONtoolsError
            else:
                logging.error("Something wrong happened when EXONtools was trying to create the output directory. Try to provide another path for output directory")
                raise EXONtoolsError("Cannot create the output directory")


class getinput(object):
    """This class administrates all input paths (-i option)"""

    extensions = ['.fa', '.fasta', '.fna', '.fq', '.fastq', '.gz']
    critical = True
    ignore_empty_files = True
    parentdir = None

    def __init__(self, inpath):
        FileList = getinput.test_inpath(inpath)
        if getinput.ignore_empty_files:
            FileList2 = [x for x in FileList if os.path.getsize(x)]
            if len(FileList2) != len(FileList):
                logging.warning("Some of your files were empty and were dismissed from the analysis.")
                if len(FileList2) == 0:
                    if getinput.critical:
                        logging.error("There are no input files to process after removing all empty files.")
                        raise EXONtoolsError("Input path error. See the log message")
                    else:
                        logging.warning("The input file '{0:s}' is empty!".format(os.path.basename(FileList[0])))
        self.files = FileList2

    def find_paired(self):
        """Parsing the paired files"""

        logging.debug("Parsing the paired files")

        paired = {}
        unpaired = {}

        forward_reads = [x for x in self.files if "_R2" not in x]
        reverse_reads = [x for x in self.files if "_R2" in x]

        if len(forward_reads) + len(reverse_reads) == 0:
            if getinput.critical:
                logging.error("No input files include forward and reverse read identifiers")
                logging.error("Please check that file names include '_R1' or '_R2' for paired reads correspondingly.")
                raise EXONtoolsError("Input path error. See the log message")
            else:
                logging.warning("No input files include forward and reverse read identifiers")
                return paired, unpaired

        for infileR1 in forward_reads:
            R1_libname = re.split("_|\.", os.path.basename(infileR1))[0]
            check = False
            for infileR2 in reverse_reads:
                R2_libname = re.split("_|\.", os.path.basename(infileR2))[0]
                if "_R1" in os.path.basename(infileR1) and R1_libname == R2_libname:
                    paired[R1_libname + "_paired"] = [infileR1, infileR2]
                    check = True
                    break
            if check:
                reverse_reads.remove(infileR2)
            else:
                R1_libname = R1_libname + "_unpaired"
                if R1_libname not in unpaired:
                    unpaired[R1_libname] = [infileR1]
                elif getinput.critical:
                    logging.error("Several files with unpaired reads are found for the same '{0:s}' library".format(R1_libname.split("_")[0]))
                    logging.error("Please correct your input file names according to the rules of the EXONtools pipeline")
                    raise EXONtoolsError("Input path error. See the log message")
                else:
                    while R1_libname in unpaired:
                        R1_libname = R1_libname + "_"
                    unpaired[R1_libname] = [infileR1]

        if reverse_reads:
            for infileR2 in reverse_reads:
                R2_libname = re.split("_|\.", os.path.basename(infileR2))[0] + "_unpaired"
                if R2_libname not in unpaired:
                    unpaired[R2_libname] = [infileR2]
                elif getinput.critical:
                    logging.error("Several files with unpaired reads are found for the same '{0:s}' library".format(R2_libname.split("_")[0]))
                    logging.error("Please correct your input file names according to the rules of the EXONtools pipeline")
                    raise EXONtoolsError("Input path error. See the log message")
                else:
                    while R2_libname in unpaired:
                        R2_libname = R2_libname + "_"
                    unpaired[R2_libname] = [infileR2]

        if unpaired:
            logging.info("The folowing files do not have their matching pair:")
            for infile in sorted(unpaired.keys()):
                logging.info(" --> " + os.path.basename(unpaired[infile][0]))

        logging.debug("Paired files were successfully parsed: OK")
        return paired, unpaired

    @classmethod
    def format(cls, extensions=None):
        """Setting the file extension requirements"""

        if isinstance(extensions, (list,)):
            cls.extensions = extensions
        elif extensions == "default":
            cls.extensions = ['.fa', '.fasta', '.fq', '.fastq', '.gz']
        elif isinstance(extensions, (str,)):
            cls.extensions = [extensions]
        else:
            logging.error("Unknown format is provided for input file extension")
            raise EXONtoolsError("Input path error. See the log message")

    @staticmethod
    def test_inpath(inpath):
        """"Testing the input path parameters"""

        logging.debug("Testing the input path parameters")

        inpath = os.path.realpath(inpath)

        if os.path.isfile(inpath):
            if getinput.check_ext(inpath):
                logging.info("Provided input data --> {0:s}".format(os.path.basename(inpath)))
                getinput.parentdir = os.path.dirname(inpath)
                return [inpath]
            else:
                if getinput.critical:
                    logging.error("Please check your input file format")
                    logging.error(" Only files with the following extensions are allowed in the current EXONtools command:\n{0:s}".format(" ".join(getinput.extensions)))
                    raise EXONtoolsError("Input path error. See the log message")
                else:
                    logging.warning("Your input file format is not supported. Returninig empty list.")
                    return []
        elif os.path.isdir(inpath):
            FileList = [os.path.join(inpath, x) for x in os.listdir(inpath) if any(list(map(lambda p: x.endswith(p), getinput.extensions)))]
            if len(FileList) == 0:
                if getinput.critical:
                    logging.error("Directory '{0:s}' does not contain supported files. Please correct the input path or check that input files have have one of the following extensions: {1:s}".format(inpath, " ".join(getinput.extensions)))
                    raise EXONtoolsError("Input path error. See the log message")
                else:
                    logging.warning("Directory '{0:s}' does not contain supported files.".format(inpath))
                    return FileList
            FileList.sort(key=natural_sort)
            #logging.debug("{0:d} files are provided for the analysis".format(len(FileList)))
            logging.info("The following files are found in")
            logging.info("{0:s}:".format(inpath))
            for infile in FileList:
                logging.info(" --> {0:s}".format(os.path.basename(infile)))
            getinput.parentdir = inpath
            return FileList

        else:
            if getinput.critical:
                logging.error("Directory or file '{0:s}' does not exist.\nPlease check your input path".format(inpath))
                raise EXONtoolsError("Path '{0:s}' does not exist".format(inpath))
            else:
                logging.warning("Directory or file '{0:s}' does not exist".format(inpath))

    @staticmethod
    def check_ext(infile):
        for ext in getinput.extensions:
            if os.path.basename(infile).endswith(ext):
                return True

    @classmethod
    def setdry(cls):
        getinput.critical = False

    @classmethod
    def undry(cls):
        getinput.critical = True


def parseinput(inpath, forward, reverse):
    """ This function takes input files from multiple option arguments and combines them into two dictionaries with paired and unpaired files """

    logging.debug("Parsing input files")

    if (inpath and forward) or (inpath and reverse):
        logging.error("Options -i and -R1/2 are mutually exclusive and cannot be used together")
        raise EXONtoolsError("Input path error. See the log message")
    elif inpath:
        inpath_valid = getinput(inpath)
        paired, unpaired = inpath_valid.find_paired()
    elif forward and reverse:
        if os.path.isfile(forward) and os.path.isfile(reverse):
            forward_valid = getinput(forward)
            reverse_valid = getinput(reverse)
            forward_valid.files = forward_valid.files + reverse_valid.files
            paired, unpaired = forward_valid.find_paired()
            if unpaired:
                logging.error("Paired files must have compatible name identifiers. Verify provided file paths in -R1 and -R2 options")
                raise EXONtoolsError("Input path error. See the log message")
        else:
            logging.error("Paths in -R1 and -R2 options must lead to real files. Verify your input settings")
            raise EXONtoolsError("File does not exist")
    elif forward:
        if os.path.isfile(forward):
            forward_valid = getinput(forward)
            paired, unpaired = forward_valid.find_paired()
        else:
            logging.error("File does not exist. Paths in -R1 option must lead to real file")
            raise EXONtoolsError("Input path error. See the log message")
    elif reverse:
        logging.error("Reverse reads (-R2 option) can only be used together with -R1")
        raise EXONtoolsError("Input path error. See the log message")
    else:
        logging.error("No data was provided for analysis. Verify your input settings")
        raise EXONtoolsError("Input path error. See the log message")

    logging.debug("Inpit file were successfully parsed: OK")
    return paired, unpaired


class makenewdir(object):
    """This class creates and contols the new directory within the output directory or by provided path"""

    dryrun = False

    def __init__(self, name, fullname=""):

        if os.path.expanduser(os.path.dirname(name)):
            self.path = name
            self.name = os.path.basename(name)
            self.fullname = fullname
            logging.debug("Creating the {0:s} directory {1:s}".format(self.fullname, self.path))

        elif output.path:
            self.name = name
            self.fullname = fullname
            self.path = os.path.join(output.path, self.name)
            logging.debug("Creating the {0:s} directory {1:s}".format(self.fullname, self.path))
        else:
            logging.error("The directory can be initialized only after the output directory and name assignation.")
            raise EXONtoolsError("The output directory must be initialized first")

        if not makenewdir.dryrun:
            try:
                os.mkdir(self.path)
                logging.debug("The {0:s} directory was created: OK".format(self.fullname))
            except OSError as e:
                if e.errno == 17 and not makenewdir.dryrun:
                    shutil.rmtree(self.path)
                    os.mkdir(self.path)
                else:
                    logging.error("Some error occured during creation of the {0:s} directory".format(self.fullname))
                    raise EXONtoolsError("Cannot create the directory")

    def delete(self):
        logging.debug("Deleting the {0:s} directory {1:s}".format(self.fullname, self.path))
        try:
            if not makenewdir.dryrun:
                shutil.rmtree(self.path)
            logging.debug("The {0:s} directory was deleted: OK".format(self.fullname))
        except OSError as e:
            if e.errno == 2:
                pass
            else:
                logging.error("Something wrong happened when EXONtools was trying to delete the {0:s} directory".format(self.fullname))
                raise EXONtoolsError("Cannot delete the directory")
        except TypeError:
            logging.error("Cannot delete the {0:s} directory because it has not been assigned yet".format(self.fullname))
            raise EXONtoolsError("Cannot delete directory from non-existing path")

    @classmethod
    def setdry(cls):
        makenewdir.dryrun = True
