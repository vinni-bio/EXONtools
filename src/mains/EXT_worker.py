# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function
import logging


class worker(object):
    """This class will initiate the task performed by a single ('worker') in parallel analysis"""

    tasks_performed = 0
    dryrun = False

    def __init__(self, func, params):
        self.function = func
        self.params = params
        self.command = "EXONtools function <{0}>: {1}".format(self.function.__name__, self.params)
        worker.tasks_performed += 1

    def run(self):
        logging.debug("Worker is running the function '{0:s}'".format(self.function.__name__))
        logging.debug(self.command)
        if not worker.dryrun:
            return self.function(*self.params)

    @classmethod
    def reset(cls):
        cls.tasks_performed = 0
        logging.debug("Reset workers: OK")

    @classmethod
    def setdry(cls):
        worker.dryrun = True
