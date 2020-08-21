# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import sys
import os
import platform
import importlib
import subprocess
import logging
from numbers import Number

from mains.EXT_errors import EXONtoolsError


def testpy():
    """Checks python version and returns False if <2.7 or <3.5"""
    pycheck = platform.python_version_tuple()
    pycheck = [int(x) for x in pycheck]
    if (pycheck[0] == 2 and pycheck[1] > 6) or (pycheck[0] == 3 and pycheck[1] > 4):
        return True
    else:
        return False


def testpip():
    """Checks pip version and prints the warning that pip is not working correctly in Python 3"""
    try:
        import pip
    except ImportError:
        sys.exit("The PIP is not installed on your system.\nPlease run:\n'sudo easy_install pip'\nor install Python from Anaconda source (recommended)")
    return True


def import_module_test(module, force_install=False):
    """Test if module is installed. Also it can run the installation with pip using force_install"""
    if testpip():
        try:
            if "=" in module:
                modulek = module.split("=")[0]
            else:
                modulek = module
            modulek = modulek.replace("python-", "")
            importlib.import_module(modulek)
        except (ImportError):
            if module == "igraph":
                module = "python-igraph"
            if force_install:
                print("Installing '{}' python module with pip to your local bin directory.\nThis installation will not be executed again in the following EXONtools runs".format(module))
                try:
                    subprocess.check_call(["pip", 'install', "--user", module])
                except subprocess.CalledProcessError:
                    sys.exit("EXONtools cannot install the '{0}' module with 'force_install'. Please try to install it manually:\n 'pip install [--user] {0}'".format(module))
                print("The module '{}' was succesfully installed. Now you can run the EXONtools pipeline normally.".format(module))
            else:
                sys.exit("'{0}' module is required by the EXONtools but currently it is not installed in your system.\nPlease run: 'pip install [--user] {0}' in your console.\nAlternatively you can run the command 'python EXONtools.py force_install'".format(module))


def memory_check(memvalue, minmem=1, maxpercmem=0.85):
    """Function checks that the requested memory is sufficient for chosen command"""

    from psutil import virtual_memory

    if isinstance(memvalue, Number):
        pass
    else:
        logging.error("Memory value must have integer or float type format")
        raise TypeError

    mem = memvalue * 1000
    maxmem = virtual_memory()[0] / 1000000 * 0.85

    if memvalue < minmem:
        logging.error("The current pipeline step requires at least 1 Gb of memory.\n{0:s} Gb was requested by user. Please set a lower value\n".format(str(memvalue)))
        raise EXONtoolsError
    elif mem > maxmem:
        logging.error("Memory used by the current pipeline step cannot exceed {0:d} %.\n{1:s} Gb was requested by user. Please set a lower value\n".format(int(maxpercmem * 100), str(memvalue)))
        raise EXONtoolsError
    else:
        logging.info("{0:s} Gb memory was allocated for the current pipeline step".format(str(memvalue)))
        logging.debug("Memory check: OK")
        return True


def positive(values, strict=True):
    """Validates if provided values are all positive"""

    if not isinstance(values, list):
        logging.error("Values for 'Positive' validation must be provided as a list")
        raise EXONtoolsError("Positive value check error")
    elif values:
        for val in values:
            if isinstance(val, Number):
                if val <= 0 and strict:
                    logging.error("The provided value '{0:s}' is not positive. Aborting the operation.".format(str(val)))
                    raise EXONtoolsError("Positive value check error")
                if val < 0 and not strict:
                    logging.error("The provided value '{0:s}' is not zero or positive. Aborting the operation.".format(str(val)))
                    raise EXONtoolsError("Positive value check error")
            else:
                logging.error("The provided value '{0:s}' is not of a numeric type. Aborting the operation.".format(str(val)))
                raise EXONtoolsError("Positive value check error")
            logging.debug("Positive value check for '{0:s}': OK".format(str(val)))
        return True
    else:
        logging.error("Nothing to validate. Please provide the correct arguments. Aborting the operation.")
        raise EXONtoolsError("Positive value check error")


def pathexists(paths):
    """Returns True if all paths exist"""

    if not isinstance(paths, list):
        logging.error("Values for 'PathExists' validation must be provided as a list")
        raise EXONtoolsError("Non-existing path error")
    elif paths:
        for path_ in paths:
            if not os.path.exists(path_):
                logging.error("The path {0:s} does not exist".format(path_))
                raise EXONtoolsError("Non-existing path error")
        return True
    else:
        logging.error("Nothing to validate. Please provide the correct arguments. Aborting the operation.")
        raise EXONtoolsError("Non-existing path error")
