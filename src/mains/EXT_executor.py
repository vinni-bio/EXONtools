# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2019
# Copyright 2019 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.


from __future__ import print_function
from configparser import RawConfigParser, ConfigParser
import logging
import subprocess
import os
import sys
import shutil

from mains.EXT_errors import EXONtoolsError
from mains import EXT_validator


def DefaultConfig():
    """Writes the default config file"""

    config = RawConfigParser(allow_no_value=True)
    cfg_path = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "dependencies.ini")

    config.add_section("PATHS")
    config.set("PATHS", "########## PROVIDE LOCAL PATHS FOR THE DEPENDENCIES HERE ##########")
    config.set("PATHS", "#### OR CHANGE TO A BLANK FIELD IF YOU DON'T NEED THE PROGRAM #####")
    config.set("PATHS", "ABYSS", "default")
    config.set("PATHS", "BBTOOLS", "default")
    config.set("PATHS", "BLAST", "default")
    config.set("PATHS", "BLAT", "default")
    config.set("PATHS", "BOWTIE2", "default")
    config.set("PATHS", "BWA", "default")
    config.set("PATHS", "CAP3", "default")
    config.set("PATHS", "CDHIT", "default")
    config.set("PATHS", "CUTADAPT", "default")
    config.set("PATHS", "FASTQC", "default")
    config.set("PATHS", "MULTIQC", "default")
    config.set("PATHS", "SAMTOOLS", "default")
    config.set("PATHS", "SPADES", "default")
    config.set("PATHS", "TRIMMOMATIC", "default")
    config.set("PATHS", "TRINITY", "default")
    config.set("PATHS", "TRANSABYSS", "default")

    config.add_section("FORMATS")
    config.set("FORMATS", "########## PLEASE DON'T CHANGE ANYTHING HERE ##########")
    config.set("FORMATS", "ABYSS", "file")
    config.set("FORMATS", "BBTOOLS", "directory")
    config.set("FORMATS", "BLAST", "directory")
    config.set("FORMATS", "BLAT", "file")
    config.set("FORMATS", "BOWTIE2", "directory")
    config.set("FORMATS", "BWA", "file")
    config.set("FORMATS", "CAP3", "file")
    config.set("FORMATS", "CDHIT", "directory")
    config.set("FORMATS", "CUTADAPT", "file")
    config.set("FORMATS", "FASTQC", "file")
    config.set("FORMATS", "MULTIQC", "file")
    config.set("FORMATS", "SAMTOOLS", "file")
    config.set("FORMATS", "SPADES", "file")
    config.set("FORMATS", "TRIMMOMATIC", "file")
    config.set("FORMATS", "TRINITY", "file")
    config.set("FORMATS", "TRANSABYSS", "file")

    config.add_section("PROGRAMS")
    config.set("PROGRAMS", "########## PLEASE DON'T CHANGE ANYTHING HERE ##########")
    config.set("PROGRAMS", "ABYSS", ["abyss-pe"])
    config.set("PROGRAMS", "BBTOOLS", ["bbduk.sh", "bbmerge.sh"])
    config.set("PROGRAMS", "BLAST", ["makeblastdb", "blastn", "blastx", "blastp", "tblastn", "tblastx"])
    config.set("PROGRAMS", "BLAT", ["blat"])
    config.set("PROGRAMS", "BOWTIE2", ["bowtie2", "bowtie2-build"])
    config.set("PROGRAMS", "BWA", ["bwa"])
    config.set("PROGRAMS", "CAP3", ["cap3"])
    config.set("PROGRAMS", "CDHIT", ["cd-hit-est"])
    config.set("PROGRAMS", "CUTADAPT", ["cutadapt"])
    config.set("PROGRAMS", "FASTQC", ["fastqc"])
    config.set("PROGRAMS", "MULTIQC", ["multiqc"])
    config.set("PROGRAMS", "SAMTOOLS", ["samtools"])
    config.set("PROGRAMS", "SPADES", ["spades.py"])
    config.set("PROGRAMS", "TRIMMOMATIC", ["trimmomatic"])
    config.set("PROGRAMS", "TRINITY", ["trinity"])
    config.set("PROGRAMS", "TRANSABYSS", ["transabyss"])

    config.add_section("VERSIONS")
    config.set("VERSIONS", "########## OPTIONAL DEPENDENCY VERSION CHECK ##########")
    config.set("VERSIONS", "ABYSS", "2")
    config.set("VERSIONS", "BBTOOLS", "38.00")
    config.set("VERSIONS", "BLAST", "2.7")
    config.set("VERSIONS", "BLAT", "NA")
    config.set("VERSIONS", "BOWTIE2", "2.2")
    config.set("VERSIONS", "BWA", "0.7")
    config.set("VERSIONS", "CAP3", "NA")
    config.set("VERSIONS", "CDHIT", "4.2")
    config.set("VERSIONS", "CUTADAPT", "2.8")
    config.set("VERSIONS", "FASTQC", "0.11")
    config.set("VERSIONS", "MULTIQC", "1.8")
    config.set("VERSIONS", "SAMTOOLS", "1.7")
    config.set("VERSIONS", "SPADES", "3.0")
    config.set("VERSIONS", "TRIMMOMATIC", "0.36")
    config.set("VERSIONS", "TRINITY", "2.6")
    config.set("VERSIONS", "TRANSABYSS", "2")

    with open(cfg_path, 'w') as cfg_file:
        config.write(cfg_file)


def getlist(option, sep=',', chars=' []"\''):
    """Return a list from a ConfigParser option. By default,
       split on a comma and strip whitespaces."""
    return [chunk.strip(chars) for chunk in option.split(sep)]


class executor(object):
    """A class to execute the dependency programs"""

    DEFAULT_CONFIG_FILEPATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dependencies.ini")
    dryrun = False
    debug = False
    config_set = False
    commandTemplates = {}
    prognames = {}
    program_paths = {}
    versions = {}

    def __init__(self, program, params, conditions=None, custom_arg_string="", quiet=False):
        """Initalizes Executor object.  Reads EXONtools.cfg and loads configuration settings, builds the command
        string, and configures stdIO pipes """

        if conditions is None:
            conditions = {}
        if not executor.config_set:
            executor.loadConfigs()
            executor.initalize_commands()

        self.program = program
        self.extra_args = custom_arg_string
        self.conditions = conditions
        self.command = "{0:s} {1:s}".format(executor.commandTemplates[program.upper()].format(*params), custom_arg_string)
        self.quiet = quiet

    def run_program(self):
        if not executor.program_paths[self.program.lower()]:
            logging.error("The path for the program {0:s} is not provided. Please verify your config file".format(self.program))
            raise EXONtoolsError("Dependency path error")
        executor.validate_conditions(self.conditions)
        if executor.debug:
            output = None
        else:
            output = open(os.devnull, 'wb')
        if not self.quiet:
            logging.info("Running '{0:s}' program".format(self.program))
        if executor.dryrun:
            logging.info("Implementing the dry run mode. No output will be provided")
            logging.info(self.command)
        else:
            logging.debug(self.command)
            try:
                subprocess.check_call(self.command, shell=True, stdout=output, stderr=output)
                if not executor.debug:
                    output.close()
            except subprocess.CalledProcessError:
                logging.error("Cannot start the '{0:s}' program".format(self.program))
                raise EXONtoolsError("Program execution error")

    @classmethod
    def loadConfigs(cls):
        """Parsing the configuration file with dependency paths"""

        config = ConfigParser()
        if not os.path.isfile(executor.DEFAULT_CONFIG_FILEPATH):
            logging.warning("EXONtools config file was not found!")
            logging.warning("Writing the default config settings to {0:s}".format(executor.DEFAULT_CONFIG_FILEPATH))
            DefaultConfig()
        config.read(executor.DEFAULT_CONFIG_FILEPATH)
        logging.info("Loading and testing configuration settings for all dependency programs")

        opts = config.options("PATHS")
        vers = config.options("VERSIONS")

        for sec in config.sections():
            for opt in opts:
                if opt not in config.options(sec):
                    logging.error("All program identifiers must be identical in all sections")
                    raise EXONtoolsError("Configuration file error")

        for opt in opts:
            if config.get("PATHS", opt).strip() == "" or config.get("PATHS", opt) is None:
                logging.warning("The path for '{0:s}' is not provided in the configuration file".format(opt))
                for prog in getlist(config.get("PROGRAMS", opt)):
                    executor.program_paths[prog.lower()] = None
            elif config.get("FORMATS", opt).lower() == "file":
                if len(getlist(config.get("PROGRAMS", opt))) != 1:
                    logging.error("Configuration file error. Only one program is allowed when program path leads to a file")
                    raise EXONtoolsError("Configuration file error")
                else:
                    prog = getlist(config.get("PROGRAMS", opt))[0].lower()
                    executor.prognames[prog] = opt
                    if config.get("PATHS", opt).lower() == "default":
                        p = subprocess.Popen('which ' + prog, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
                        default_path, _ = p.communicate()
                        if default_path and os.path.exists(default_path.strip()):
                            executor.program_paths[prog.lower()] = default_path.strip()
                            logging.debug("Setting the path for program '{0:s}' to {1:s}".format(prog, default_path.strip()))
                        else:
                            logging.error("The command '{0:s}' is not found in your system PATH".format(prog))
                            logging.error("Please install the program and provide the correct program path in the 'dependencies.ini' file")
                            raise EXONtoolsError("Configuration file error")
                    else:
                        prog_path = os.path.abspath(os.path.expanduser(config.get("PATHS", opt)))
                        if os.path.exists(prog_path):
                            executor.program_paths[prog.lower()] = prog_path
                            logging.debug("Setting the path for program '{0:s}' to {1:s}".format(opt, prog_path))
                        else:
                            logging.error("The path {1:s} for '{0:s}' program does not exist".format(prog, prog_path))
                            logging.error("Please install the program and provide the correct program path in the dependencies.ini")
                            raise EXONtoolsError("Configuration file error")
            elif config.get("FORMATS", opt).lower() == "directory":
                if len(getlist(config.get("PROGRAMS", opt))) < 1:
                    logging.error("At least one program must be listed for '{0:s}' PATH".format(opt))
                    raise EXONtoolsError("Configuration file error")
                else:
                    for prog in getlist(config.get("PROGRAMS", opt)):
                        executor.prognames[prog] = opt
                    if config.get("PATHS", opt).lower() == "default":
                        for prog in getlist(config.get("PROGRAMS", opt)):
                            p = subprocess.Popen('which ' + prog, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
                            default_path, _ = p.communicate()
                            if default_path and os.path.exists(default_path.strip()):
                                executor.program_paths[prog.lower()] = default_path.strip()
                                logging.debug("Setting the path for program '{0:s}' to {1:s}".format(prog, default_path.strip()))
                            else:
                                logging.error("The command '{0:s}' is not found in your system PATH".format(prog))
                                logging.error("Please install the program and provide the correct program path in the 'dependencies.ini' file")
                                raise EXONtoolsError("Configuration file error")
                    else:
                        for prog in getlist(config.get("PROGRAMS", opt)):
                            prog_path = os.path.join(os.path.abspath(os.path.expanduser(config.get("PATHS", opt))), prog)
                            if os.path.exists(prog_path):
                                executor.program_paths[prog.lower()] = prog_path
                                logging.debug("Setting the path for program '{0:s}' to {1:s}".format(opt, prog_path))
                            else:
                                logging.error("The path {1:s} for '{0:s}' program does not exist".format(prog, prog_path))
                                logging.error("Please install the program and provide the correct program path in the 'dependencies.ini' file")
                                raise EXONtoolsError("Configuration file error")
            else:
                logging.error("Unknown format is provided for '{0:s}'. Please correct the configuration file 'dependencies.ini'".format(opt))
                raise EXONtoolsError("Configuration file error. Unknown format")

        logging.debug("Fixing PBLAT path")
        try:
            if executor.program_paths["blat"].endswith('pblat'):
                newpath = os.path.abspath(os.path.expanduser(os.path.join(os.path.dirname(executor.program_paths["blat"]), 'blat')))
                if not os.path.exists(newpath):
                    shutil.copy2(executor.program_paths["blat"], newpath)
                executor.program_paths["blat"] = newpath
        except (KeyError, AttributeError):
            pass
        except (OSError, IOError):
            logging.error("The pblat path cannot be fixed")
            raise EXONtoolsError

        logging.debug("Uploading versions of the dependency programs")
        for prog in vers:
            ver = config.get("VERSIONS", prog).strip()
            executor.versions[prog.lower()] = ver

        executor.config_set = True
        logging.debug("Dependency PATHs were all tested: OK")

    @staticmethod
    def validate_conditions(conditions):
        """Validates this program's conditions for execution"""

        dryexceptions = ["pathexists"]

        for condition in conditions.items():
            if executor.dryrun and condition[0] in dryexceptions:
                pass
            else:
                getattr(EXT_validator, condition[0])(condition[1])

    @classmethod
    def setdry(cls):
        executor.dryrun = True

    @classmethod
    def setdebug(cls):
        executor.debug = True

    @classmethod
    def setconfig(cls, *args):
        executor.loadConfigs()
        if args:
            executor.progcheck(*args)
        executor.initalize_commands()

    @staticmethod
    def progcheck(*progs):
        if executor.program_paths:
            logging.debug("Program execution check")
            for prog in progs:
                try:
                    if executor.program_paths[prog.lower()] is None:
                        logging.error("The path for the program '{0:s}' is not provided in dependencies.ini".format(prog))
                        raise EXONtoolsError("Program execution check error")
                    elif not os.path.exists(executor.program_paths[prog.lower()]):
                        logging.error("The path for program '{0:s}' does not exist".format(prog))
                        raise EXONtoolsError("Program execution check error")
                    else:
                        executor.check_version(executor.prognames[prog.lower()], executor.versions[executor.prognames[prog.lower()]])
                        logging.debug("Program '{0:s}' check: OK".format(prog))
                except KeyError:
                    logging.error("The program '{0:s}' is not listed in the dependencies.ini file".format(prog))
                    raise EXONtoolsError("Program execution check error")
        else:
            logging.error("The configuration file with dependencies is not initialized. Please run executor.loadConfigs before checking the program path")
            raise EXONtoolsError("Program execution check error")

    @staticmethod
    def check_version(program, version):
        """Checks the dependency program version"""

        arg_templates = {
            "abyss": (" --version", lambda x: x.strip().split("\n")[0].split()[-1]),
            "bbtools": (" -v", lambda x: x.strip().split("Version: ")[1].split('\n')[0].strip()),
            "bowtie2": (" --version", lambda x: x.strip().split("\n")[0].split()[-1]),
            "blast": (" -version", lambda x: x.strip().split("\n")[0].split()[-1]),
            "bwa": ("", lambda x: x.strip().split("Version: ")[1].split('\n')[0].strip()),
            "cdhit": ("", lambda x: x.strip().split()[3]),
            "cutadapt": (" --version", lambda x: x.strip()),
            "fastqc": (" -v", lambda x: x.strip().split()[1]),
            "multiqc": (" --version", lambda x: x.strip().split()[2]),
            "samtools": (" --version", lambda x: x.strip().split("\n")[0].split()[-1]),
            "spades.py": (" --version", lambda x: x.strip().split("\n")[0].split()[-1]),
            "transabyss": (" --version", lambda x: x.strip()),
            "trimmomatic": (" -version", lambda x: x.strip()),
            "trinity": (" --version", lambda x: x.strip().split("\n")[0].split()[-1]),
            "uname": (" -a", lambda x: x.strip().split()[2])
        }

        if version.strip() != "NA":
            try:
                if program == "abyss":
                    version_com = os.path.join(os.path.pardir(executor.program_paths["abyss-pe"]), "ABYSS") + arg_templates[program][0]
                elif program == "cdhit":
                    version_com = executor.program_paths["cd-hit-est"] + arg_templates[program][0]
                elif program == "blast":
                    version_com = executor.program_paths["blastn"] + arg_templates[program][0]
                elif program == "bbtools":
                    version_com = executor.program_paths["bbduk.sh"] + arg_templates[program][0]
                else:
                    version_com = executor.program_paths[program] + arg_templates[program][0]

                p = subprocess.Popen(version_com, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                stdoutdata, stderrdata = p.communicate()
                if stdoutdata.strip():
                    check = arg_templates[program][1](stdoutdata)
                elif stderrdata.strip():
                    check = arg_templates[program][1](stderrdata)
                else:
                    logging.error("Probably you found the bug: the error occurred during parsing of dependency versions".format(program))
                    raise EXONtoolsError("Dependency program version check error")

                if check == "__BLEEDING_EDGE__":
                    return True

                ver = version.strip().replace(".", "")
                checkint = "".join([x for x in check if x.isdigit()])
                if int(ver) <= int(checkint[0:len(ver)]):
                    logging.debug("Version test for '{}' program: OK".format(program))
                    return True
                else:
                    logging.error("Version of the '{0:s}' program is not supported by the EXONtools pipeline".format(program))
                raise EXONtoolsError("Dependency program version check error")
            except KeyError:
                logging.error("The program '{0:s}' is not in the list of dependency version checks\nPlease use NA as its version identifier in the dependencies.ini or just delete/comment that line\n".format(program))
                raise EXONtoolsError("Dependency program version check error")
            except TypeError:
                pass
            except subprocess.CalledProcessError:
                logging.error("Fails the version check for '{0:s}' program. Please verify its path in the configuration file".format(program))
                raise EXONtoolsError("Dependency program version check error")

    @classmethod
    def exportpath(cls, *args):
        """exports path for directory of the dependency program"""
        for arg in args:
            p = subprocess.Popen('which ' + arg.lower(), shell=True, stdout=subprocess.PIPE, universal_newlines=True)
            default_path, _ = p.communicate()
            if default_path and os.path.exists(default_path):
                logging.warning("The default path for '{0:s}' program exists. Nothing to export.".format(arg))
            elif not executor.program_paths[arg.lower()]:
                logging.error("The path for the program {0:s} is not provided. Please verify your config file".format(arg))
                raise EXONtoolsError("Dependency program export path error")
            else:
                try:
                    dirname = os.path.dirname(executor.program_paths[arg.lower()])
                    os.environ['PATH'] = dirname + ":" + os.environ['PATH']
                    logging.debug("Exporting path for '{}' parent directory".format(arg))
                except KeyError:
                    logging.error("The path for the program {0:s} is not loaded. Please verify your config file".format(arg))
                    raise EXONtoolsError("Dependency program export path error")

    @classmethod
    def initalize_commands(cls):
        """Initalizes the commandTemplates dictionary."""
        cls.commandTemplates = {}

        logging.debug("Initializing commands")
        arg_templates = {
            # ./abyss-pe --directory= mpirun= np= j= k= c= e= n= E=0 name= in= se=
            "abyss-pe": " --directory={0:s} {1:s} j={2:d} {3:s} E=0 name={4:s} {5:s}",
            # "bbduk"  MEMORY %s %s THREADS qin=PHRED
            "bbduk.sh": " -Xmx{0:d}m {1:s} {2:s} t={3:d} qin={4:s}",
            # "bbmerge": MEMORY THREADS in1=INFILE1 in2=INFILE2 out=MERGED out1=OUTFILE1 out2=OUTFILE2...
            "bbmerge.sh": "-Xmx{0:d}m -t={1:d} in1={2:s} in2={3:s} out={4:s} outu1={5:s} outu2={6:s} interleaved=f ordered=t minoverlap={7:d} mininsert={8:d}",
            # ./blastn -db DBDIR/NAME -query INFILE -out OUTFILE -outfmt 6 -evalue 1e-10 -num_threads INT OTHER PARS
            "blastn": " -db {0:s} -query {1:s} -out {2:s} -outfmt {3:d} -evalue {4:s} -num_threads {5:d} {6:s}",
            # ./blastp -db DBDIR/NAME -query INFILE -out OUTFILE -outfmt 6 -evalue 1e-10 -num_threads INT OTHER PARS
            "blastp": " -db {0:s} -query {1:s} -out {2:s} -outfmt {3:d} -evalue {4:s} -num_threads {5:d} {6:s}",
            # ./blastx -db DBDIR/NAME -query INFILE -out OUTFILE -outfmt 6 -evalue 1e-10 -num_threads INT OTHER PARS
            "blastx": " -db {0:s} -query {1:s} -out {2:s} -outfmt {3:d} -evalue {4:s} -num_threads {5:d} {6:s}",
            # ./pblat DATABASE QUERY -noHead -dots=100 -threads=INT OUTPUT
            "blat": " {0:s} {1:s} -noHead -dots=100 -threads={2:d} {3:s}",
            # ./bowtie2-build REFERENCE_PATH INDEXPATH/name --quiet
            "bowtie2-build": " --quiet {0:s} {1:s}",
            # ./bowtie2 --threads -x INDEX -S NAME INPUTFILES PARS --phred -q
            "bowtie2": " --threads {0:d} -x {1:s} -S {2:s} {3:s} {4:s} --phred{5:d} -q",
            # ./bwa CMD PARAMETERS INFILES
            "bwa": " {0:s} {1:s} {2:s}",
            # ./cap3 INFILE -p 85 -f 1500 -e 500 -z 1 -o 16 -s 251
            # "cap3": " {0:s} -p 75 -f 1500 -e 500 -z 1 -o 100 -s 500",
            "cap3": " {0:s} -p 85 -f 1500 -e 500 -z 1 -o 30 -s 251",
            # ./cd-hit-est  -i INFILE -o OUTFILE -c IDENTITY -l MINLEN -M MEMORY -T THREADS -r 1 -d 30
            "cd-hit-est": " -i {0:s} -o {1:s} -c {2:0.2f} -l {3:d} -M {4:d} -T {5:d} -r 1 -d 30",
            # ./cutadapt -f fastq -e 0.15 -O 4 -n 5 OPTIONS OUTPATH INPATH
            "cutadapt": " -f fastq -e 0.15 -O 4 -n 5 {0:s} {1:s} {2:s}",
            # ./fastqc -t THREADS --quiet -o OUTPATH --extract -f FASTQ INPATH
            "fastqc": " -t {0:d} --quiet -o {1:s} --extract -f {2:s} {3:s}",
            # ./makeblastdb -in INFILE -out DBDIR/NAME -dbtype TYPE
            "makeblastdb": " -in {0:s} -out {1:s} -dbtype {2:s}",
            # ./multiqc -q -i TITLE -b EXONtools -o OUTPUT
            "multiqc": " -q -i EXONtools -b '{0:s}' -o {1:s} --export --module {2:s} {3:s}",
            # ./samtools CMD OPTIONS
            "samtools": " {0:s} {1:s}",
            # ./spades.py --disable-gzip-output --only-assembler -k -cov-cutoff -o -t -m -1 -2 -s
            "spades.py": " --disable-gzip-output {0:s} -o {1:s} -t {2:d} -m {3:d} {4:s}",
            # ./tblastn -db DBDIR/NAME -query INFILE -out OUTFILE -outfmt 6 -evalue 1e-10 -num_threads INT OTHER PARS
            "tblastn": " -db {0:s} -query {1:s} -out {2:s} -outfmt {3:d} -evalue {4:s} -num_threads {5:d} {6:s}",            # ./tblastx -db DBDIR/NAME -query INFILE -out OUTFILE -outfmt 6 -evalue 1e-10 -num_threads INT OTHER PARS
            "tblastx": " -db {0:s} -query {1:s} -out {2:s} -outfmt {3:d} -evalue {4:s} -num_threads {5:d} {6:s}",
            # ./transabyss --se --pe --name --outdir --kmer --cov --seros 0 --eros --pairs --threads --mpi --noref
            "transabyss": " {0:s} --name {1:s} --outdir {2:s} --seros 0 {3:s} --threads {4:d} --mpi {5:d} --noref",
            # " " + self.program_paths["trimmomatic"] + " %s -threads %d -phred%s -trimlog %s %s %s"
            "trimmomatic": "java -jar",
            # ./Trinity --seqType fq --max_memory --CPU --no_normalize_reads --left --right --output
            "trinity": " --full_cleanup --seqType fq --max_memory {0:d}G {1:s} --CPU {2:d} {3:s} --output {4:s}"
        }

        for prog in executor.program_paths:
            try:
                cls.commandTemplates.update({prog.upper(): cls.program_paths[prog] + arg_templates[prog]})
            except KeyError:
                logging.critical("The command template for '{0:s}' program does not exist".format(prog))
                raise EXONtoolsError("Command template error")
            except (AttributeError, TypeError):
                pass
