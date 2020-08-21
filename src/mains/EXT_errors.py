# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.


class Error(Exception):
    """Base class for other exceptions"""
    pass


class EXONtoolsError(Error):
    """Raise EXONtools exception"""

    def __init__(self, message="The EXONtools run is failed. Read the last log message."):
        self.message = message

    def __str__(self):
        return repr(self.message)
